// measure_batched_ks.cpp
#include "openfhe.h"

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include "utils/serial.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>

#if defined(__GLIBC__)
  #include <malloc.h>
  static inline void TrimMalloc() { malloc_trim(0); }
#else
  static inline void TrimMalloc() {}
#endif

using namespace lbcrypto;

// -------- Serialization helpers (size in bytes) --------
template <class T>
static size_t SerializedSizeBytesBinary(const T& obj) {
    std::stringstream ss;
    try {
        lbcrypto::Serial::Serialize(obj, ss, lbcrypto::SerType::BINARY); 
    } catch (const std::exception&) {
        return 0; 
    }
    std::streampos pos = ss.tellp();
    return (pos >= 0) ? static_cast<size_t>(pos) : 0;
}

static size_t CiphertextSizeBytes(const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ct) {
    if (!ct) return 0;
    return SerializedSizeBytesBinary(ct);
}

static inline double ToMB(size_t bytes) {
    return static_cast<double>(bytes) / 1'000'000.0; // SI MB
}

// Rotation keys size (all automorphism keys currently loaded in cc)
// static size_t RotationKeysSizeBytes(const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& cc) {
//     std::stringstream ss;
//     bool ok = cc->SerializeEvalAutomorphismKey(ss, lbcrypto::SerType::BINARY); 
//     if (!ok) return 0;
//     std::streampos pos = ss.tellp();
//     return (pos >= 0) ? static_cast<size_t>(pos) : 0;
// }


// -------- Helpers: bit-length from NativeInteger --------
static inline double BitsOf(const NativeInteger& x){
#if 1
    return static_cast<double>(x.GetMSB());
#else
    return std::floor(std::log2(x.ConvertToDouble()) + 1e-12);
#endif
}

static std::string JoinBits(const std::vector<double>& v){
    std::ostringstream os;
    for (size_t i=0;i<v.size();++i){
        if (i) os << ' ';
        os << static_cast<int>(std::round(v[i]));
    }
    return os.str();
}

// -------- Vector helpers --------
static std::vector<double> MakeMsg128(){
    std::vector<double> v(128);
    for (size_t i=0;i<128;++i) v[i] = 0.001 * static_cast<double>(i+1); // 0.001..0.128
    return v;
}
static std::vector<double> RollLeft(const std::vector<double>& v, int k){
    const int n = static_cast<int>(v.size());
    std::vector<double> out(n);
    int s = ((k % n) + n) % n;
    for (int i=0;i<n;++i) out[i] = v[(i + s) % n];
    return out;
}
static double MaxAbsErr(const std::vector<double>& a, const std::vector<double>& b){
    const size_t n = std::min(a.size(), b.size());
    double m = 0.0;
    for (size_t i=0;i<n;++i) m = std::max(m, std::abs(a[i]-b[i]));
    return m;
}

// -------- Encrypt zero --------
static Ciphertext<DCRTPoly> EncryptZero(const CryptoContext<DCRTPoly>& cc,
                                        const PublicKey<DCRTPoly>& pk,
                                        size_t len = 128) {
    // create a vector of zeros
    std::vector<double> zeros(len, 0.0);
    auto ptx_zero = cc->MakeCKKSPackedPlaintext(zeros);
    return cc->Encrypt(pk, ptx_zero);
}

// -------- P info (k_P, per-prime bits, sum logP) --------
struct PInfo {
    size_t kP = 0;
    std::vector<double> bits; // per p_j
    double sumLogP = 0.0;
};
static PInfo GetPInfo(const CryptoContext<DCRTPoly>& cc) {
    PInfo out;
    auto base = cc->GetCryptoParameters();
    auto params = std::dynamic_pointer_cast<CryptoParametersRNS>(base);
    if (!params) return out;
    auto pParams = params->GetParamsP(); // {p_1,...,p_k} in HYBRID/BATCHED
    if (!pParams) return out;
    const auto& towers = pParams->GetParams();
    out.kP = towers.size();
    out.bits.reserve(out.kP);
    for (size_t j = 0; j < out.kP; ++j) {
        auto pj = towers[j]->GetModulus();
        double b = BitsOf(pj);
        out.bits.push_back(b);
        out.sumLogP += b;
    }
    return out;
}

// -------- Q info (l_Q, per-prime bits, sum logQ) --------
struct QInfo {
    size_t lQ = 0;
    std::vector<double> bits; // per q_i
    double sumLogQ = 0.0;
};
static QInfo GetQInfo(const CryptoContext<DCRTPoly>& cc) {
    QInfo out;
    auto base = cc->GetCryptoParameters();
    auto qParams = base->GetElementParams(); // Q-only basis
    if (!qParams) return out;
    const auto& towers = qParams->GetParams();
    out.lQ = towers.size();
    out.bits.reserve(out.lQ);
    for (size_t i = 0; i < out.lQ; ++i) {
        auto qi = towers[i]->GetModulus();
        double b = BitsOf(qi);
        out.bits.push_back(b);
        out.sumLogQ += b;
    }
    return out;
}

// -------- Parameter presets --------
struct Preset {
    std::string name;       // e.g., "N=2^14"
    uint32_t ringDim;       // Ring dimension
    uint32_t multDepth;     // L
    uint32_t scalingBits;   // ScalingModSize
    uint32_t firstModBits;  // FirstModSize (0 = unused)
    SecurityLevel sec;
    uint32_t batchSize;
};
static std::vector<Preset> MakePresets(){
    return {
        {"N=2^14", 1<<14, 7, 34, 46, HEStd_128_classic, 128},
        {"N=2^15", 1<<15, 13, 40, 51, HEStd_128_classic, 128},
        {"N=2^16", 1<<16, 24, 45, 56, HEStd_128_classic, 128},
        // {"N=2^17", 1<<17, 38, 50, 61, HEStd_128_classic, 128}
    };
}

// -------- Build context --------
static CryptoContext<DCRTPoly> BuildContext(const Preset& ps, KeySwitchTechnique ksTech){
    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(ps.sec);
    p.SetRingDim(ps.ringDim);
    p.SetMultiplicativeDepth(ps.multDepth);
    p.SetScalingModSize(ps.scalingBits);
    if (ps.firstModBits) p.SetFirstModSize(ps.firstModBits);
    p.SetBatchSize(ps.batchSize);
    p.SetKeySwitchTechnique(ksTech);    // HYBRID for baseline, BATCHED for lazy+batchedKS
    p.SetNumLargeDigits(3);             // dnum=3
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);
    cc->Enable(KEYSWITCH);
    return cc;
}

// -------- Rotation keys --------
static void GenRotateKeys_Baseline(const CryptoContext<DCRTPoly>& cc,
                                   const PrivateKey<DCRTPoly>& sk,
                                   int k) {
    if (k <= 0) return;
    std::vector<int32_t> rots;
    rots.reserve(static_cast<size_t>(k));
    for (int i = 1; i <= k; ++i) rots.push_back(static_cast<int32_t>(i));
    cc->EvalRotateKeyGen(sk, rots);
}

static void GenRotateKeys_Lazy(const CryptoContext<DCRTPoly>& cc,
                               const PrivateKey<DCRTPoly>& sk,
                               int k) {
    if (k <= 0) return;
    std::vector<int32_t> rots;
    rots.reserve(static_cast<size_t>(k));
    for (int i = 1; i <= k; ++i) rots.push_back(static_cast<int32_t>(i));
    cc->EvalLazyRotateKeyGen(sk, rots);
}

// -------- Warm-ups to reduce timing noise --------
static void WarmUpBaseline(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk) {
    auto msg = MakeMsg128();
    auto ptx = cc->MakeCKKSPackedPlaintext(msg);
    auto ct  = cc->Encrypt(pk, ptx);
    for (int i=0;i<8;++i) ct = cc->EvalRotate(ct, 1);
}

static void WarmUpBatched(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk) {
    auto msg = MakeMsg128();
    auto ptx = cc->MakeCKKSPackedPlaintext(msg);
    auto ct  = cc->Encrypt(pk, ptx);
    for (int i=0;i<8;++i) ct = cc->EvalLazyRotate(ct, 1);
}

// -------- reference helpers (sum of k rotations) --------
static std::vector<double> MakeRefSum(const std::vector<double>& msg, int k) {
    std::vector<double> ref(msg.size(), 0.0);
    for (int i = 1; i <= k; ++i) {
        auto ri = RollLeft(msg, i);
        for (size_t j = 0; j < ref.size(); ++j) ref[j] += ri[j];
    }
    return ref;
}

static std::vector<double> DecryptToVec(const CryptoContext<DCRTPoly>& cc,
                                        const PrivateKey<DCRTPoly>& sk,
                                        const Ciphertext<DCRTPoly>& ct,
                                        size_t len = 128) {
    Plaintext out;
    cc->Decrypt(sk, ct, &out);
    out->SetLength(len);
    return out->GetRealPackedValue();
}

struct MeasOut {
    double ms = 0.0;
    Ciphertext<DCRTPoly> ct;
    ksprofile::Phases phases; 
    size_t ctBytes = 0;
};

// -------- Baseline: time k individual rotations (EvalAdd outside timer) --------
static MeasOut MeasureBaseline(const CryptoContext<DCRTPoly>& cc,
                               const PublicKey<DCRTPoly>& pk,
                               int k) {
    auto msg = MakeMsg128();
    auto ptx = cc->MakeCKKSPackedPlaintext(msg);
    auto ct0 = cc->Encrypt(pk, ptx);

    ksprofile::Reset();

    auto sum = EncryptZero(cc, pk);
    double ms_rot = 0.0;

    for (int i = 1; i <= k; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        auto r  = cc->EvalRotate(ct0, i);
        auto t1 = std::chrono::steady_clock::now();
        ms_rot +=  std::chrono::duration<double, std::milli>(t1 - t0).count();

        sum = cc->EvalAdd(sum, r);
    }
    auto phases = ksprofile::Snapshot();

    size_t bytes = CiphertextSizeBytes(sum);

    return {ms_rot, sum, phases, bytes};
}


// -------- Batched: time k lazy rotations + 1 EvalBatchedKS (EvalLazyAdd outside rotation timer) --------
static MeasOut MeasureBatched(const CryptoContext<DCRTPoly>& cc,
                              const PublicKey<DCRTPoly>& pk,
                              int k) {
    auto msg = MakeMsg128();
    auto ptx = cc->MakeCKKSPackedPlaintext(msg);
    auto ct0 = cc->Encrypt(pk, ptx);

    auto sum = EncryptZero(cc, pk);

    double ms_lazy_rot = 0.0;
    for (int i = 1; i <= k; ++i) {
        auto t0 = std::chrono::steady_clock::now();
        auto ri = cc->EvalLazyRotate(ct0, i);           // automorphism only
        auto t1 = std::chrono::steady_clock::now();
        ms_lazy_rot += std::chrono::duration<double, std::milli>(t1 - t0).count();

        sum = cc->EvalLazyAdd(sum, ri);
    }

    size_t bytesBeforeKS = CiphertextSizeBytes(sum);

    ksprofile::Reset();
    auto k0 = std::chrono::steady_clock::now();
    sum = cc->EvalBatchedKS(sum);                       // KS timed
    auto k1 = std::chrono::steady_clock::now();
    auto phases = ksprofile::Snapshot();

    double ms_total = ms_lazy_rot
                    + std::chrono::duration<double, std::milli>(k1 - k0).count();

    return {ms_total, sum, phases, bytesBeforeKS};
}

int main(int argc, char** argv){
    // args: [trials]
    const int trials = (argc>=2 ? std::max(1, std::atoi(argv[1])) : 10);

    const std::vector<int> K = {2,4,8,16,32,64};
    auto presets = MakePresets();

    std::ofstream csv("timings.csv");
    csv << "preset,N,depth_L,scalingBits,firstModBits,dnum,ksTech,"
        "k,trial,"
        "ms_baseline,ms_batched,"
        "base_modup_ms,base_inner_ms,base_moddown_ms,base_modup_calls,base_inner_calls,base_moddown_calls,"
        "bat_modup_ms,bat_inner_ms,bat_moddown_ms,bat_modup_calls,bat_inner_calls,bat_moddown_calls,"
        "err_baseline,err_batched,err_between,"
        "P_k,sum_logP_bits,P_bits,"
        "Q_l,sum_logQ_bits,Q_bits,"
        "ct_MB_baseline,ct_MB_batched_preKS,rotkeys_MB_base,rotkeys_MB_batched\n";

    for (const auto& ps : presets){
        std::cout << "=== " << ps.name << " | RingDim(set)="<<ps.ringDim
                  << " | L="<<ps.multDepth
                  << " | scale="<<ps.scalingBits
                  << " | first="<<ps.firstModBits
                  << " | dnum=3 | KS(Base=HYBRID, Batch=BATCHED) ===\n";

        for (int k : K){
            {
                lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
            }
            auto ccBatch = BuildContext(ps, BATCHED);
            auto kpBatch = ccBatch->KeyGen();
            GenRotateKeys_Lazy(ccBatch, kpBatch.secretKey, k);
            size_t rotBytesBatch = 0; // TODO

            auto ccBase  = BuildContext(ps, HYBRID);
            auto kpBase  = ccBase->KeyGen();
            GenRotateKeys_Baseline(ccBase, kpBase.secretKey, k);
            size_t rotBytesBase  = 0; // TODO

            // Verify ring dimension
            uint32_t actualRingDim = ccBatch->GetRingDimension();
            if (actualRingDim != ps.ringDim) {
                std::cerr << "WARNING: Ring dimension mismatch! Set=" << ps.ringDim
                          << ", Actual=" << actualRingDim << "\n";
            }

            auto pinfo = GetPInfo(ccBatch);
            auto qinfo = GetQInfo(ccBatch);
            std::cout << "[k="<<k<<"] N(actual)="<<actualRingDim
                      << " | P_k="<<pinfo.kP
                      << ", sum log2(P)≈"<<std::fixed<<std::setprecision(1)<<pinfo.sumLogP
                      << " bits, P bits: ["<<JoinBits(pinfo.bits)<<"]\n";
            std::cout << "         Q_l="<<qinfo.lQ
                      << ", sum log2(Q)≈"<<std::fixed<<std::setprecision(1)<<qinfo.sumLogQ
                      << " bits, Q bits: ["<<JoinBits(qinfo.bits)<<"]\n";

            // Collect timing statistics across trials
            std::vector<double> base_times, batched_times;
            base_times.reserve(trials);
            batched_times.reserve(trials);
            double max_err_base = 0.0, max_err_batched = 0.0, max_err_between = 0.0;

            for (int t=0;t<trials;++t){
                // Warm-ups (per context)
                WarmUpBaseline(ccBase,  kpBase.publicKey);
                WarmUpBatched(ccBatch,  kpBatch.publicKey);

                // Measure
                auto base    = MeasureBaseline(ccBase,  kpBase.publicKey,  k);
                auto batched = MeasureBatched(ccBatch,  kpBatch.publicKey, k);

                // Reference and errors
                const auto msg = MakeMsg128();
                const auto ref = MakeRefSum(msg, k);

                auto base_vec  = DecryptToVec(ccBase,  kpBase.secretKey,  base.ct,   128);
                auto batch_vec = DecryptToVec(ccBatch, kpBatch.secretKey, batched.ct,128);

                double err_baseline = MaxAbsErr(base_vec,  ref);
                double err_batched  = MaxAbsErr(batch_vec, ref);
                double err_between  = MaxAbsErr(base_vec,  batch_vec);

                // Collect statistics
                base_times.push_back(base.ms);
                batched_times.push_back(batched.ms);
                max_err_base = std::max(max_err_base, err_baseline);
                max_err_batched = std::max(max_err_batched, err_batched);
                max_err_between = std::max(max_err_between, err_between);

                // CSV
                csv << ps.name << "," << ccBatch->GetRingDimension() << ","
                    << ps.multDepth << "," << ps.scalingBits << "," << ps.firstModBits << ","
                    << 3 << "," << "HYBRID_vs_BATCHED" << ","
                    << k << "," << t << ","
                    << base.ms << "," << batched.ms << ","
                    << base.phases.modup_ms   << "," << base.phases.inner_ms   << "," << base.phases.moddown_ms   << ","
                    << base.phases.modup_calls<< "," << base.phases.inner_calls<< "," << base.phases.moddown_calls<< ","
                    << batched.phases.modup_ms   << "," << batched.phases.inner_ms   << "," << batched.phases.moddown_ms   << ","
                    << batched.phases.modup_calls<< "," << batched.phases.inner_calls<< "," << batched.phases.moddown_calls<< ","
                    << err_baseline << "," << err_batched << "," << err_between << ","
                    << pinfo.kP << "," << pinfo.sumLogP << ",\""
                    << JoinBits(pinfo.bits) << "\","
                    << qinfo.lQ << "," << qinfo.sumLogQ << ",\""
                    << JoinBits(qinfo.bits) << "\","
                    << std::fixed << std::setprecision(2)
                    << ToMB(base.ctBytes) << "," << ToMB(batched.ctBytes) << ","
                    << ToMB(rotBytesBase) << "," << ToMB(rotBytesBatch) << "\n";

                if ((t % 5) == 4) {
                    csv.flush();
                    TrimMalloc();
                }
            }

            // Compute and print summary statistics for this k
            auto compute_stats = [](const std::vector<double>& v) {
                double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
                double sq_sum = 0.0;
                for (double x : v) sq_sum += (x - mean) * (x - mean);
                double stddev = std::sqrt(sq_sum / v.size());
                return std::make_pair(mean, stddev);
            };

            auto [base_mean, base_std] = compute_stats(base_times);
            auto [bat_mean, bat_std] = compute_stats(batched_times);
            double speedup = base_mean / std::max(1e-9, bat_mean);

            std::cout << "k="<<std::setw(3)<<k<<" | "<<trials<<" trials completed\n"
                      << "  baseline: "<<std::fixed<<std::setprecision(3)<<base_mean<<" ± "<<base_std<<" ms\n"
                      << "  batched:  "<<bat_mean<<" ± "<<bat_std<<" ms\n"
                      << "  speedup:  "<<std::setprecision(2)<<speedup<<"x\n"
                      << "  max errors: base="<<std::setprecision(6)<<max_err_base
                      << ", batched="<<max_err_batched
                      << ", between="<<max_err_between<<"\n";

            std::cout << std::endl;

            ccBase.reset();
            ccBatch.reset();
            lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
            TrimMalloc();
        }
    }

    std::cout << "\nResults written to timings.csv\n";
    return 0;
}
