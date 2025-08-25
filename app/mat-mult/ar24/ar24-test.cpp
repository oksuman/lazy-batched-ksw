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
#include "ar24.h"
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

using namespace lbcrypto;

#if defined(__GLIBC__)
  #include <malloc.h>
  static inline void TrimMalloc() { malloc_trim(0); }
#else
  static inline void TrimMalloc() {}
#endif

static const int   kTrials = 100;
static const unsigned kSeed = 1337u;
static const int   kDims[] = {8, 16, 32, 64};
static const char* kCSV   = "ar24_bench_results.csv";

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

struct ContextPack {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
    int s{0};
    int num_slots{0};
};

static ContextPack makeContextForD(int d, KeySwitchTechnique ksTech) {
    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(2);
    P.SetRingDim(1<<14);
    P.SetScalingModSize(50);
    P.SetSecurityLevel(HEStd_128_classic);
    P.SetKeySwitchTechnique(ksTech);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    int max_batch = (int)cc->GetRingDimension() / 2;
    int s = std::min(max_batch / (d * d), d);
    if (s <= 0) s = 1;
    P.SetBatchSize(d * d * s);

    cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    return {cc, kp.publicKey, kp.secretKey, s, d * d * s};
}

static Ciphertext<DCRTPoly>
encryptMatrix(const CryptoContext<DCRTPoly>& cc,
                   const PublicKey<DCRTPoly>& pk,
                   const std::vector<double>& M, int d, int s) {
    const int num_slots = d * d * s;
    std::vector<double> slots(num_slots, 0.0);
    std::copy(M.begin(), M.end(), slots.begin());
    auto pt = cc->MakeCKKSPackedPlaintext(slots);
    return cc->Encrypt(pk, pt);
}

struct IAlgo {
    virtual ~IAlgo() = default;
    virtual bool isLazy() const = 0;
    virtual void plan_and_keygen(const CryptoContext<DCRTPoly>& cc,
                                 const PrivateKey<DCRTPoly>& sk,
                                 int num_slots) = 0;
    virtual Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                           const Ciphertext<DCRTPoly>& B) = 0;
};

struct BaselineAdapter : IAlgo {
    MATMULT_AR24 impl;
    BaselineAdapter(const CryptoContext<DCRTPoly>& cc,
                    const PublicKey<DCRTPoly>& pk, int d) : impl(cc, pk, d) {}
    bool isLazy() const override { return false; }
    void plan_and_keygen(const CryptoContext<DCRTPoly>& cc,
                         const PrivateKey<DCRTPoly>& sk,
                         int num_slots) override {
        RotationKeyCollector rk;
        rk.begin(num_slots, false);
        impl.eval_mult_plan(rk);
        rk.generate(cc, sk);
    }
    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                   const Ciphertext<DCRTPoly>& B) override {
        return impl.eval_mult(A, B);
    }
};

struct LazyAdapter : IAlgo {
    MATMULT_AR24 impl;
    LazyAdapter(const CryptoContext<DCRTPoly>& cc,
                const PublicKey<DCRTPoly>& pk, int d) : impl(cc, pk, d) {}
    bool isLazy() const override { return true; }
    void plan_and_keygen(const CryptoContext<DCRTPoly>& cc,
                        const PrivateKey<DCRTPoly>& sk,
                        int num_slots) override {
        RotationKeyCollectorLazy rk;
        rk.begin(num_slots);
        impl.eval_mult_lazy_plan(rk);
        rk.generate(cc, sk);
    }
    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& A,
                                   const Ciphertext<DCRTPoly>& B) override {
        return impl.eval_mult_lazy(A, B);
    }
};

template <typename MakeAlgo>
static std::tuple<Stats, Acc>
run_one_method(const std::string& label, int d, int trials, unsigned seed, KeySwitchTechnique ksTech,
               MakeAlgo makeAlgo) {
    ContextPack ctx = makeContextForD(d, ksTech);
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

    {
        algo->plan_and_keygen(ctx.cc, ctx.sk, ctx.num_slots);
    }

    {
        auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[0], d, ctx.s);
        auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[0], d, ctx.s);
        auto ctC = algo->eval_mult(ctA, ctB);
        Plaintext pt;
        ctx.cc->Decrypt(ctx.sk, ctC, &pt);
        pt->SetLength(d * d);
        (void)pt->GetRealPackedValue();
    }

    std::vector<double> times_ms; times_ms.reserve(trials);
    double max_abs_err_overall = 0.0; long double mse_sum = 0.0L;

    for (int t = 0; t < trials; ++t) {
        auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[t], d, ctx.s);
        auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[t], d, ctx.s);

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

    return {st, ac};
}

struct Row { int d; Stats st; Acc ac; };

static std::vector<Row> run_baseline_suite(const std::vector<int>& dims, int trials, unsigned seed) {
    std::vector<Row> rows;
    std::cout << "==== Baseline suite ====\n";
    for (int d : dims) {
        auto [st, ac] = run_one_method(
            "Baseline   ", d, trials, seed, HYBRID,
            [](const CryptoContext<DCRTPoly>& cc,
               const PublicKey<DCRTPoly>& pk, int d_) {
                return std::make_unique<BaselineAdapter>(cc, pk, d_);
            });
        rows.push_back({d, st, ac});
        lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
        TrimMalloc();
    }
    return rows;
}

static std::vector<Row> run_batched_suite(const std::vector<int>& dims, int trials, unsigned seed) {
    std::vector<Row> rows;
    std::cout << "==== Lazy suite ====\n";
    for (int d : dims) {
        auto [st, ac] = run_one_method(
            "Batched KS ", d, trials, seed, BATCHED,
            [](const CryptoContext<DCRTPoly>& cc,
               const PublicKey<DCRTPoly>& pk, int d_) {
                return std::make_unique<LazyAdapter>(cc, pk, d_);
            });
        rows.push_back({d, st, ac});
        lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
        TrimMalloc();
    }
    return rows;
}

int main() {
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    TrimMalloc();

    std::vector<int> dims(std::begin(kDims), std::end(kDims));

    auto batRows = run_batched_suite(dims, kTrials, kSeed);
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    TrimMalloc();
    auto baseRows = run_baseline_suite(dims, kTrials, kSeed);


    std::ofstream csv(kCSV);
    csv << "d,base_mean_ms,base_std_ms,bat_mean_ms,bat_std_ms,"
           "base_max_abs_err,base_mse,bat_max_abs_err,bat_mse,speedup_x\n";
    csv << std::fixed << std::setprecision(6);

    std::cout << "==== Summary ====\n";
    for (size_t i = 0; i < dims.size(); ++i) {
        const int d = dims[i];
        const auto& b = baseRows[i];
        const auto& t = batRows[i];
        const double speedup = b.st.mean_ms / std::max(1e-12, t.st.mean_ms);

        std::cout << "d=" << d
                  << " | base " << std::setprecision(3) << b.st.mean_ms << " ms"
                  << " vs bat " << t.st.mean_ms << " ms"
                  << " | speedup " << speedup << "x\n";

        csv << d << ","
            << b.st.mean_ms << "," << b.st.std_ms << ","
            << t.st.mean_ms << "," << t.st.std_ms << ","
            << b.ac.max_abs_err << "," << b.ac.mse << ","
            << t.ac.max_abs_err << "," << t.ac.mse << ","
            << std::setprecision(3) << speedup << std::setprecision(6) << "\n";
    }
    csv.close();

    std::cout << "Results written to " << kCSV << "\n";
    return 0;
}
