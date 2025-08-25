#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <fstream>

#include "openfhe.h"
#include "jkls18.h"

using namespace lbcrypto;

// -------- Toggle experiments (uncomment to enable) --------
// #define RUN_BASELINE 1
// #define RUN_LAZY 1
#define RUN_HOIST 1
// ---------------------------------------------------------

#if defined(__GLIBC__)
  #include <malloc.h>
  static inline void TrimMalloc() { malloc_trim(0); }
#else
  static inline void TrimMalloc() {}
#endif

static inline void HardReset() {
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    TrimMalloc();
}

// ================= User-tunable =================
static const int      kTrials = 100;
static const unsigned kSeed   = 1337u;
// static const int      kDims[] = {8};
static const int      kDims[] = {8, 16, 32, 64};
static const char*    kCSV    = "jkls18_bench_results.csv";
// ================================================

// ---------------- Timing/accuracy utils ----------------
struct Stats { double mean_ms{0.0}; double std_ms{0.0}; };
static Stats meanStd(const std::vector<double>& v) {
    if (v.empty()) return {};
    double m = std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    double var = 0.0; for (double x : v) var += (x - m)*(x - m);
    var /= (double)v.size();
    return {m, std::sqrt(var)};
}

struct Acc { double max_abs_err{0.0}; double mse{0.0}; };
static Acc accuracy(const std::vector<double>& a, const std::vector<double>& b) {
    const size_t n = std::min(a.size(), b.size());
    double maxe = 0.0; long double sse = 0.0L;
    for (size_t i = 0; i < n; ++i) {
        double e = std::abs(a[i] - b[i]);
        if (e > maxe) maxe = e;
        sse += (long double)e * (long double)e;
    }
    return {maxe, (double)(sse / (long double)std::max<size_t>(1, n))};
}

static std::vector<double> randomMatrix(int d, std::mt19937& gen) {
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    std::vector<double> M(d * d);
    for (double& x : M) x = dis(gen);
    return M;
}

static std::vector<double> matmulPlain(const std::vector<double>& A,
                                       const std::vector<double>& B, int d) {
    std::vector<double> C(d * d, 0.0);
    for (int i = 0; i < d; ++i)
        for (int k = 0; k < d; ++k) {
            double aik = A[i * d + k];
            for (int j = 0; j < d; ++j)
                C[i * d + j] += aik * B[k * d + j];
        }
    return C;
}

// ---------------- Context builder ----------------
struct ContextPack {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
    int num_slots{0}; // = d*d
};

static ContextPack makeContextForD(int d, KeySwitchTechnique ksTech, bool useLazy) {
    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(3);
    P.SetRingDim(1<<14);
    P.SetScalingModSize(50);
    P.SetSecurityLevel(HEStd_128_classic);
    P.SetBatchSize(d * d); // exactly d*d

    P.SetKeySwitchTechnique(ksTech); // HYBRID (baseline/hoist) or BATCHED (lazy)

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    // Rotation keys: include 0; cover [-d*d+1, ..., d*d-1]
    std::vector<int32_t> rotIndices;
    rotIndices.reserve(2 * d * d - 1);
    for (int i = -d*d + 1; i < d*d; ++i) rotIndices.push_back(i);

    if (useLazy) cc->EvalLazyRotateKeyGen(kp.secretKey, rotIndices);
    else         cc->EvalRotateKeyGen(kp.secretKey,     rotIndices);

    return {cc, kp.publicKey, kp.secretKey, d * d};
}

// ---------------- Encryption helper ----------------
static Ciphertext<DCRTPoly>
encryptMatrix(const CryptoContext<DCRTPoly>& cc,
              const PublicKey<DCRTPoly>& pk,
              const std::vector<double>& M) {
    auto pt = cc->MakeCKKSPackedPlaintext(M);
    return cc->Encrypt(pk, pt);
}

// ---------------- Adapters ----------------
struct IAlgo {
    virtual ~IAlgo() = default;
    virtual std::string getName() const = 0;
    virtual Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                           const Ciphertext<DCRTPoly>& B) = 0;
};

#ifdef RUN_BASELINE
struct BaselineAdapter : IAlgo {
    MATMULT_JKLS18 impl;
    BaselineAdapter(const CryptoContext<DCRTPoly>& cc,
                    const PublicKey<DCRTPoly>& pk, int d) : impl(cc, pk, d) {}
    std::string getName() const override { return "Baseline"; }
    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                   const Ciphertext<DCRTPoly>& B) override {
        return impl.eval_mult(A, B);
    }
};
#endif

#ifdef RUN_LAZY
struct LazyAdapter : IAlgo {
    MATMULT_JKLS18 impl;
    LazyAdapter(const CryptoContext<DCRTPoly>& cc,
                const PublicKey<DCRTPoly>& pk, int d) : impl(cc, pk, d) {}
    std::string getName() const override { return "Lazy"; }
    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                   const Ciphertext<DCRTPoly>& B) override {
        return impl.eval_mult_lazy(A, B); // lazy rotates + batched KS inside
    }
};
#endif

#ifdef RUN_HOIST
struct HoistAdapter : IAlgo {
    MATMULT_JKLS18 impl;
    HoistAdapter(const CryptoContext<DCRTPoly>& cc,
                 const PublicKey<DCRTPoly>& pk, int d) : impl(cc, pk, d) {}
    std::string getName() const override { return "Hoist"; }
    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                   const Ciphertext<DCRTPoly>& B) override {
        return impl.eval_mult_hoist(A, B);
    }
};
#endif

// ---------------- Runner ----------------
template <typename MakeAlgo>
static std::tuple<Stats, Acc>
run_one_method(const std::string& label, int d, int trials, unsigned seed,
               MakeAlgo makeAlgo,
               KeySwitchTechnique ksTech, bool useLazy) {
    ContextPack ctx = makeContextForD(d, ksTech, useLazy);
    std::unique_ptr<IAlgo> algo = makeAlgo(ctx.cc, ctx.pk, d);

    std::mt19937 gen(seed);
    std::vector<std::vector<double>> A_set, B_set, C_ref_set;
    A_set.reserve(trials); B_set.reserve(trials); C_ref_set.reserve(trials);
    for (int t = 0; t < trials; ++t) {
        auto A = randomMatrix(d, gen);
        auto B = randomMatrix(d, gen);
        auto C = matmulPlain(A, B, d);
        A_set.emplace_back(std::move(A));
        B_set.emplace_back(std::move(B));
        C_ref_set.emplace_back(std::move(C));
    }

    // warm-up
    {
        auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[0]);
        auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[0]);
        auto ctC = algo->eval_mult(ctA, ctB);
        Plaintext pt;
        ctx.cc->Decrypt(ctx.sk, ctC, &pt);
        pt->SetLength(d * d);
        (void)pt->GetRealPackedValue();
    }

    std::vector<double> times_ms; times_ms.reserve(trials);
    double max_abs_err_overall = 0.0; long double mse_sum = 0.0L;

    for (int t = 0; t < trials; ++t) {
        auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[t]);
        auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[t]);

        auto t0 = std::chrono::steady_clock::now();
        auto ctC = algo->eval_mult(ctA, ctB);
        auto t1 = std::chrono::steady_clock::now();

        Plaintext pt;
        ctx.cc->Decrypt(ctx.sk, ctC, &pt);
        pt->SetLength(d * d);
        std::vector<double> dec = pt->GetRealPackedValue();

        auto acc = accuracy(dec, C_ref_set[t]);
        if (acc.max_abs_err > max_abs_err_overall) max_abs_err_overall = acc.max_abs_err;
        mse_sum += acc.mse;

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        times_ms.push_back(ms);
    }

    Stats st = meanStd(times_ms);
    Acc ac{max_abs_err_overall, (double)(mse_sum / trials)};

    std::cout << std::fixed << std::setprecision(3);
    std::cout << label << " : mean " << st.mean_ms << " ms, std " << st.std_ms
              << " ms | max|err| " << std::setprecision(6) << ac.max_abs_err
              << ", mse " << ac.mse << std::setprecision(3) << "\n";

    ctx.cc->ClearEvalAutomorphismKeys();
    ctx.cc->ClearEvalMultKeys();
    return {st, ac};
}

// ---------------- Suites ----------------
struct Row { int d; Stats st; Acc ac; };

#ifdef RUN_BASELINE
static std::vector<Row> run_baseline_suite(const std::vector<int>& dims, int trials, unsigned seed) {
    std::vector<Row> rows;
    std::cout << "==== Baseline suite ====\n";
    for (int d : dims) {
        auto [st, ac] = run_one_method(
            "Baseline   ", d, trials, seed,
            [](const CryptoContext<DCRTPoly>& cc,
               const PublicKey<DCRTPoly>& pk, int d_) {
                return std::make_unique<BaselineAdapter>(cc, pk, d_);
            },
            HYBRID, /*useLazy=*/false
        );
        rows.push_back({d, st, ac});
        HardReset();
    }
    return rows;
}
#endif

#ifdef RUN_LAZY
static std::vector<Row> run_lazy_suite(const std::vector<int>& dims, int trials, unsigned seed) {
    std::vector<Row> rows;
    std::cout << "==== Lazy suite ====\n";
    for (int d : dims) {
        auto [st, ac] = run_one_method(
            "Lazy       ", d, trials, seed,
            [](const CryptoContext<DCRTPoly>& cc,
               const PublicKey<DCRTPoly>& pk, int d_) {
                return std::make_unique<LazyAdapter>(cc, pk, d_);
            },
            BATCHED, /*useLazy=*/true
        );
        rows.push_back({d, st, ac});
        HardReset();
    }
    return rows;
}
#endif

#ifdef RUN_HOIST
static std::vector<Row> run_hoist_suite(const std::vector<int>& dims, int trials, unsigned seed) {
    std::vector<Row> rows;
    std::cout << "==== Hoist suite ====\n";
    for (int d : dims) {
        auto [st, ac] = run_one_method(
            "Hoist      ", d, trials, seed,
            [](const CryptoContext<DCRTPoly>& cc,
               const PublicKey<DCRTPoly>& pk, int d_) {
                return std::make_unique<HoistAdapter>(cc, pk, d_);
            },
            HYBRID, /*useLazy=*/false // same keys/rotations as baseline
        );
        rows.push_back({d, st, ac});
        HardReset();
    }
    return rows;
}
#endif

// ---------------- main ----------------
int main() {
    HardReset();
    std::vector<int> dims(std::begin(kDims), std::end(kDims));

#ifdef RUN_LAZY
    auto lazyRows  = run_lazy_suite(dims, kTrials, kSeed);
    HardReset();
#endif
#ifdef RUN_BASELINE
    auto baseRows  = run_baseline_suite(dims, kTrials, kSeed);
    HardReset();
#endif
#ifdef RUN_HOIST
    auto hoistRows = run_hoist_suite(dims, kTrials, kSeed);
    HardReset();
#endif

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);

    // CSV header
    csv << "d";
#ifdef RUN_BASELINE
    csv << ",base_mean_ms,base_std_ms,base_max_abs_err,base_mse";
#endif
#ifdef RUN_LAZY
    csv << ",lazy_mean_ms,lazy_std_ms,lazy_max_abs_err,lazy_mse";
#endif
#ifdef RUN_HOIST
    csv << ",hoist_mean_ms,hoist_std_ms,hoist_max_abs_err,hoist_mse";
#endif
#if defined(RUN_BASELINE) && defined(RUN_LAZY)
    csv << ",speedup_lazy_x";
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
    csv << ",speedup_hoist_x";
#endif
    csv << "\n";

    std::cout << "==== Summary ====\n";
    for (size_t i = 0; i < dims.size(); ++i) {
        const int d = dims[i];
        csv << d;

#ifdef RUN_BASELINE
        const auto& b = baseRows[i];
        csv << "," << b.st.mean_ms << "," << b.st.std_ms
            << "," << b.ac.max_abs_err << "," << b.ac.mse;
#endif
#ifdef RUN_LAZY
        const auto& l = lazyRows[i];
        csv << "," << l.st.mean_ms << "," << l.st.std_ms
            << "," << l.ac.max_abs_err << "," << l.ac.mse;
#endif
#ifdef RUN_HOIST
        const auto& h = hoistRows[i];
        csv << "," << h.st.mean_ms << "," << h.st.std_ms
            << "," << h.ac.max_abs_err << "," << h.ac.mse;
#endif

        // Console summary + optional speedups
        std::cout << "d=" << d;
#ifdef RUN_BASELINE
        std::cout << " | base "  << std::setprecision(3) << b.st.mean_ms << " ms";
#endif
#ifdef RUN_LAZY
        std::cout << " | lazy "  << l.st.mean_ms << " ms";
#endif
#ifdef RUN_HOIST
        std::cout << " | hoist " << h.st.mean_ms << " ms";
#endif

#if defined(RUN_BASELINE) && defined(RUN_LAZY)
        {
            double speedup_lazy = b.st.mean_ms / std::max(1e-12, l.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_lazy << std::setprecision(6);
            std::cout << " | speedup(lazy) " << speedup_lazy << "x";
        }
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
        {
            double speedup_hoist = b.st.mean_ms / std::max(1e-12, h.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_hoist << std::setprecision(6);
            std::cout << " | speedup(hoist) " << speedup_hoist << "x";
        }
#endif
        std::cout << "\n";
        csv << "\n";
    }
    csv.close();

    std::cout << "Results written to " << kCSV << "\n";
    return 0;
}
