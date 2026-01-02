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
#include "kcl25.h"

using namespace lbcrypto;

// -------- Toggle experiments (uncomment to enable) --------
#define RUN_BASELINE 1
#define RUN_LAZY 1
#define DEBUG_MODE 1  
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
static const int      kTrials = 1;
static const unsigned kSeed   = 1337u;
static const int      kDims[] = {4};
static const char*    kCSV    = "kcl25_bench_results.csv";
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

struct Acc { double max_abs_err{0.0}; double mse{0.0}; double identity_max_err{0.0}; double identity_mse{0.0}; };
static Acc accuracy(const std::vector<double>& a, const std::vector<double>& b) {
    const size_t n = std::min(a.size(), b.size());
    double maxe = 0.0; long double sse = 0.0L;
    for (size_t i = 0; i < n; ++i) {
        double e = std::abs(a[i] - b[i]);
        if (e > maxe) maxe = e;
        sse += (long double)e * (long double)e;
    }
    return {maxe, (double)(sse / (long double)std::max<size_t>(1, n)), 0.0, 0.0};
}

// Matrix multiplication for verification
static std::vector<double> matrixMultiply(const std::vector<double>& A, const std::vector<double>& B, int d) {
    std::vector<double> C(d * d, 0.0);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            for (int k = 0; k < d; ++k) {
                C[i * d + j] += A[i * d + k] * B[k * d + j];
            }
        }
    }
    return C;
}

// Check how close M * M^(-1) is to identity matrix
static Acc verifyIdentity(const std::vector<double>& M, const std::vector<double>& Minv, int d) {
    auto product = matrixMultiply(M, Minv, d);

    // Compare with identity matrix
    std::vector<double> identity(d * d, 0.0);
    for (int i = 0; i < d; ++i) identity[i * d + i] = 1.0;

    double maxe = 0.0; long double sse = 0.0L;
    for (int i = 0; i < d * d; ++i) {
        double e = std::abs(product[i] - identity[i]);
        if (e > maxe) maxe = e;
        sse += (long double)e * (long double)e;
    }

    return {0.0, 0.0, maxe, (double)(sse / (long double)(d * d))};
}

// Generate random invertible matrix with elements in [-1, 1]
static std::vector<double> randomInvertibleMatrix(int d, std::mt19937& gen) {
    std::uniform_real_distribution<double> off_diag_dis(-0.15, 0.15);
    std::uniform_real_distribution<double> diag_dis(0.5, 0.9);
    std::vector<double> M(d * d, 0.0);

    // Generate with diagonal dominance to ensure invertibility
    for (int i = 0; i < d; ++i) {
        // First, generate off-diagonal elements
        double row_sum = 0.0;
        for (int j = 0; j < d; ++j) {
            if (i != j) {
                M[i * d + j] = off_diag_dis(gen); // Small off-diagonal: [-0.15, 0.15]
                row_sum += std::abs(M[i * d + j]);
            }
        }
        // Diagonal must be larger than sum of absolute values of off-diagonals
        // to ensure diagonal dominance (sufficient condition for invertibility)
        double diag_min = row_sum + 0.1;
        M[i * d + i] = std::max(diag_dis(gen), diag_min); // Ensure diagonal dominance
    }
    return M;
}

// Compute matrix inverse using plain arithmetic (Gauss-Jordan)
static std::vector<double> matrixInversePlain(const std::vector<double>& A, int d) {
    std::vector<double> M = A; // Copy input
    std::vector<double> inv(d * d, 0.0);

    // Initialize inv as identity
    for (int i = 0; i < d; ++i) {
        inv[i * d + i] = 1.0;
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < d; ++i) {
        // Find pivot
        double pivot = M[i * d + i];
        if (std::abs(pivot) < 1e-10) {
            // Singular matrix or numerical issues
            std::cerr << "Warning: Near-singular matrix at row " << i << "\n";
            return inv; // Return partial result
        }

        // Scale row i
        for (int j = 0; j < d; ++j) {
            M[i * d + j] /= pivot;
            inv[i * d + j] /= pivot;
        }

        // Eliminate column i in other rows
        for (int k = 0; k < d; ++k) {
            if (k != i) {
                double factor = M[k * d + i];
                for (int j = 0; j < d; ++j) {
                    M[k * d + j] -= factor * M[i * d + j];
                    inv[k * d + j] -= factor * inv[i * d + j];
                }
            }
        }
    }

    return inv;
}

// ---------------- Context builder ----------------
struct ContextPack {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
    int num_slots{0};
};

static ContextPack makeContextForD(int d, KeySwitchTechnique ksTech, bool useLazy) {
    CCParams<CryptoContextCKKSRNS> P;

    // Dimension-specific setup
    if (d == 4) {
        // d=4: r=18, depth=48, no bootstrapping
        P.SetMultiplicativeDepth(48);
        P.SetScalingModSize(50);
        // FirstModSize not set (use default)
    } else if (d == 8) {
        // d=8: r=21, depth=34, bootstrapping enabled
        P.SetMultiplicativeDepth(34);
        P.SetScalingModSize(59);
        P.SetFirstModSize(60);
    } else if (d == 16) {
        // d=16: r=25, depth=34, bootstrapping enabled
        P.SetMultiplicativeDepth(34);
        P.SetScalingModSize(59);
        P.SetFirstModSize(60);
    } else if (d == 32) {
        // d=32: r=28, depth=29, bootstrapping enabled
        P.SetMultiplicativeDepth(29);
        P.SetScalingModSize(59);
        P.SetFirstModSize(60);
    } else {
        // d=64 or other: r=31, depth=29, bootstrapping enabled
        P.SetMultiplicativeDepth(29);
        P.SetScalingModSize(59);
        P.SetFirstModSize(60);
    }

    P.SetRingDim(1<<17);  // 2^17 = 131,072
    P.SetSecurityLevel(HEStd_128_classic);

    int max_batch = (1<<17) / 2;  // N/2
    int s = std::min(max_batch / d / d, d);
    if (s <= 0) s = 1;
    P.SetBatchSize(s * d * d);  

    P.SetKeySwitchTechnique(ksTech);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    // Bootstrapping for d >= 8
    if (d >= 8) {
        cc->EvalBootstrapSetup({4, 5}); // level budget
    }

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    if (d >= 8) {
        cc->EvalBootstrapKeyGen(kp.secretKey, max_batch);
    }

    // Rotation keys: need full range for batched operations
    int max_rot = s * d * d;
    std::vector<int32_t> rotIndices;
    for (int i = -max_rot + 1; i < max_rot; ++i) rotIndices.push_back(i);

    // if (useLazy) cc->EvalLazyRotateKeyGen(kp.secretKey, rotIndices);
    if (useLazy) cc->EvalRotateKeyGen(kp.secretKey, rotIndices);
    else         cc->EvalRotateKeyGen(kp.secretKey,     rotIndices);

    return {cc, kp.publicKey, kp.secretKey, s * d * d};
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
    virtual Ciphertext<DCRTPoly> eval_inverse(const Ciphertext<DCRTPoly>& M) = 0;
    virtual void setSecretKey(const PrivateKey<DCRTPoly>& sk) = 0;
};

// Get iteration count r based on dimension d (from inverse_newCol_test.cpp)
static int getIterations(int d) {
    switch (d) {
        case 4:  return 14;
        case 8:  return 21;
        case 16: return 25;
        case 32: return 28;
        case 64: return 31;
        default: return 20; 
    }
}

static int getMultDepth(int d) {
    switch (d) {
        case 4:  return 2 * 18 + 12; // r=18
        case 8:  return 34;
        case 16: return 34;
        case 32: return 29;
        case 64: return 29;
        default: return 34;
    }
}

#ifdef RUN_BASELINE
struct BaselineAdapter : IAlgo {
    MATINV_KCL25 impl;
    PrivateKey<DCRTPoly> sk;
    BaselineAdapter(const CryptoContext<DCRTPoly>& cc,
                    const PublicKey<DCRTPoly>& pk, int d)
        : impl(cc, pk, d, getIterations(d), getMultDepth(d)) {}
    std::string getName() const override { return "Baseline"; }
    void setSecretKey(const PrivateKey<DCRTPoly>& secret_key) override {
        sk = secret_key;
    }
    Ciphertext<DCRTPoly> eval_inverse(const Ciphertext<DCRTPoly>& M) override {
#ifdef DEBUG_MODE
        return impl.eval_inverse_debug(M, sk);
#else
        return impl.eval_inverse(M);
#endif
    }
};
#endif

#ifdef RUN_LAZY
struct LazyAdapter : IAlgo {
    MATINV_KCL25 impl;
    PrivateKey<DCRTPoly> sk;
    LazyAdapter(const CryptoContext<DCRTPoly>& cc,
                const PublicKey<DCRTPoly>& pk, int d)
        : impl(cc, pk, d, getIterations(d), getMultDepth(d)) {}
    std::string getName() const override { return "Lazy"; }
    void setSecretKey(const PrivateKey<DCRTPoly>& secret_key) override {
        sk = secret_key;
    }
    Ciphertext<DCRTPoly> eval_inverse(const Ciphertext<DCRTPoly>& M) override {
#ifdef DEBUG_MODE
        return impl.eval_inverse_lazy_debug(M, sk);
#else
        return impl.eval_inverse_lazy(M);
#endif
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
    algo->setSecretKey(ctx.sk);

    std::mt19937 gen(seed);
    std::vector<std::vector<double>> M_set, Minv_ref_set;
    M_set.reserve(trials); Minv_ref_set.reserve(trials);

    for (int t = 0; t < trials; ++t) {
        auto M = randomInvertibleMatrix(d, gen);
        auto Minv = matrixInversePlain(M, d);
        M_set.emplace_back(std::move(M));
        Minv_ref_set.emplace_back(std::move(Minv));
    }

    // warm-up
    // {
    //     auto ctM = encryptMatrix(ctx.cc, ctx.pk, M_set[0]);
    //     auto ctMinv = algo->eval_inverse(ctM);
    //     Plaintext pt;
    //     ctx.cc->Decrypt(ctx.sk, ctMinv, &pt);
    //     pt->SetLength(d * d);
    //     (void)pt->GetRealPackedValue();
    // }

    std::vector<double> times_ms; times_ms.reserve(trials);
    double max_abs_err_overall = 0.0; long double mse_sum = 0.0L;
    double identity_max_err_overall = 0.0; long double identity_mse_sum = 0.0L;

    for (int t = 0; t < trials; ++t) {
        auto ctM = encryptMatrix(ctx.cc, ctx.pk, M_set[t]);

        auto t0 = std::chrono::steady_clock::now();
        auto ctMinv = algo->eval_inverse(ctM);
        auto t1 = std::chrono::steady_clock::now();

        Plaintext pt;
        ctx.cc->Decrypt(ctx.sk, ctMinv, &pt);
        pt->SetLength(d * d);
        std::vector<double> dec = pt->GetRealPackedValue();

        auto acc = accuracy(dec, Minv_ref_set[t]);
        if (acc.max_abs_err > max_abs_err_overall) max_abs_err_overall = acc.max_abs_err;
        mse_sum += acc.mse;

        // Verify M * M^(-1) = I
        auto id_acc = verifyIdentity(M_set[t], dec, d);
        if (id_acc.identity_max_err > identity_max_err_overall) identity_max_err_overall = id_acc.identity_max_err;
        identity_mse_sum += id_acc.identity_mse;

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        times_ms.push_back(ms);
    }

    Stats st = meanStd(times_ms);
    Acc ac{max_abs_err_overall, (double)(mse_sum / trials), identity_max_err_overall, (double)(identity_mse_sum / trials)};

    std::cout << std::fixed << std::setprecision(3);
    std::cout << label << " : mean " << st.mean_ms << " ms, std " << st.std_ms
              << " ms | max|err| " << std::setprecision(6) << ac.max_abs_err
              << ", mse " << ac.mse
              << " | M*M^-1: max|err| " << ac.identity_max_err
              << ", mse " << ac.identity_mse << std::setprecision(3) << "\n";

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

// ---------------- main ----------------
int main() {
    HardReset();
    std::vector<int> dims(std::begin(kDims), std::end(kDims));

#ifdef RUN_BASELINE
    auto baseRows  = run_baseline_suite(dims, kTrials, kSeed);
    HardReset();
#endif
#ifdef RUN_LAZY
    auto lazyRows  = run_lazy_suite(dims, kTrials, kSeed);
    HardReset();
#endif

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);

    // CSV header
    csv << "d";
#ifdef RUN_BASELINE
    csv << ",base_mean_ms,base_std_ms,base_max_abs_err,base_mse,base_id_max_err,base_id_mse";
#endif
#ifdef RUN_LAZY
    csv << ",lazy_mean_ms,lazy_std_ms,lazy_max_abs_err,lazy_mse,lazy_id_max_err,lazy_id_mse";
#endif
#if defined(RUN_BASELINE) && defined(RUN_LAZY)
    csv << ",speedup_lazy_x";
#endif
    csv << "\n";

    std::cout << "==== Summary ====\n";
    for (size_t i = 0; i < dims.size(); ++i) {
        const int d = dims[i];
        csv << d;

#ifdef RUN_BASELINE
        const auto& b = baseRows[i];
        csv << "," << b.st.mean_ms << "," << b.st.std_ms
            << "," << b.ac.max_abs_err << "," << b.ac.mse
            << "," << b.ac.identity_max_err << "," << b.ac.identity_mse;
#endif
#ifdef RUN_LAZY
        const auto& l = lazyRows[i];
        csv << "," << l.st.mean_ms << "," << l.st.std_ms
            << "," << l.ac.max_abs_err << "," << l.ac.mse
            << "," << l.ac.identity_max_err << "," << l.ac.identity_mse;
#endif

        std::cout << "d=" << d;
#ifdef RUN_BASELINE
        std::cout << " | base "  << std::setprecision(3) << b.st.mean_ms << " ms";
#endif
#ifdef RUN_LAZY
        std::cout << " | lazy "  << l.st.mean_ms << " ms";
#endif

#if defined(RUN_BASELINE) && defined(RUN_LAZY)
        {
            double speedup_lazy = b.st.mean_ms / std::max(1e-12, l.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_lazy << std::setprecision(6);
            std::cout << " | speedup(lazy) " << speedup_lazy << "x";
        }
#endif
        std::cout << "\n";
        csv << "\n";
    }
    csv.close();

    std::cout << "Results written to " << kCSV << "\n";
    return 0;
}
