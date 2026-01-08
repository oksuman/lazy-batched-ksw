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
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

using namespace lbcrypto;

// -------- Toggle experiments (uncomment to enable) --------
#define RUN_BASELINE 1
#define RUN_LAZY 1
// #define DEBUG_MODE 1  
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
static const int      kDims[] = {4,8,16,32};
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

// Check invertibility (from matrix-mult-fhe/utils/matrix_utils.h)
static bool isInvertible(const std::vector<double>& matrix, int d) {
    double max_elem = 0;
    double min_elem = std::numeric_limits<double>::max();
    for (size_t i = 0; i < static_cast<size_t>(d * d); i++) {
        max_elem = std::max(max_elem, std::abs(matrix[i]));
        min_elem = std::min(min_elem, std::abs(matrix[i]));
    }
    if (min_elem < 1e-6 || max_elem / min_elem > 1e6) {
        return false;
    }
    return true;
}

// Generate random invertible matrix with elements in [-1, 1] (matrix-mult-fhe style)
static std::vector<double> randomInvertibleMatrix(int d, std::mt19937& gen) {
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    std::vector<double> M(static_cast<size_t>(d * d));

    do {
        for (int i = 0; i < d * d; i++) {
            M[static_cast<size_t>(i)] = dis(gen);
        }
    } while (!isInvertible(M, d));

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

// ---------------- Helper functions ----------------
// Get iteration count r based on dimension d (from inverse_newCol_test.cpp)
static int getIterations(int d) {
    switch (d) {
        case 4:  return 16;
        case 8:  return 20;
        case 16: return 24;
        case 32: return 28;
        case 64: return 32;
        default: return 20;
    }
}

static int getMultDepth(int d) {
    switch (d) {
        case 4:  return 41; // 16*2+9
        case 8:  return 29;
        case 16: return 29;
        case 32: return 29;
        case 64: return 29;
        default: return 34;
    }
}

// ---------------- Context builder ----------------
struct ContextPack {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
    int num_slots{0};
    int s{0};  // batch count
};

static ContextPack makeContextForD(int d, KeySwitchTechnique ksTech, bool useLazy) {
    CCParams<CryptoContextCKKSRNS> P;

    // Dimension-specific setup
    if (d == 4) {
        // d=4: r=18, depth=48, no bootstrapping
        P.SetMultiplicativeDepth(41);
        P.SetScalingModSize(50);
        // P.SetFirstModSize(60);
    } else if (d == 8) {
        // d=8: r=21, depth=34, bootstrapping enabled
        P.SetMultiplicativeDepth(29);
        P.SetScalingModSize(59);
        P.SetFirstModSize(60);
    } else if (d == 16) {
        // d=16: r=25, depth=34, bootstrapping enabled
        P.SetMultiplicativeDepth(29);
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

 
    P.SetSecurityLevel(HEStd_128_classic);

    // int max_batch = 1 << 15;
    // int s = std::min(max_batch / d / d, d);
    int s = d;
    int batchSize = d * d;
    P.SetBatchSize(batchSize);  
    P.SetKeySwitchTechnique(ksTech);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    // Bootstrapping for d >= 8
    if (d >= 8) {
        cc->EvalBootstrapSetup({3, 3}, {0, 0}, batchSize);
    }

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    if (d >= 8) {
        cc->EvalBootstrapKeyGen(kp.secretKey, batchSize);
    }

    // Use rotation collectors to generate only the needed rotation keys
    // Create a temporary MATINV_KCL25 instance to run the plan function
    MATINV_KCL25 temp_impl(cc, kp.publicKey, d, getIterations(d), getMultDepth(d));

    if (useLazy) {
        RotationKeyCollectorLazy rk;
        rk.begin(s * d * d);
        temp_impl.eval_inverse_lazy_plan(rk);
        auto collected = rk.getCollectedAutoIndices();
        std::cout << "  [Lazy] Generating optimized rotation keys\n";
        rk.generate(cc, kp.secretKey);
    } else {
        RotationKeyCollector rk;
        rk.begin(s * d * d, false);
        temp_impl.eval_inverse_plan(rk);
        std::cout << "  [Baseline] Generating optimized rotation keys\n";
        rk.generate(cc, kp.secretKey);
    }
    std::cout << "CKKS scheme is using ring dimension "
              << cc->GetRingDimension() << std::endl
              << std::endl;
    return {cc, kp.publicKey, kp.secretKey, s * d * d, s};
}

// ---------------- Encryption helper ----------------
static Ciphertext<DCRTPoly>
encryptMatrix(const CryptoContext<DCRTPoly>& cc,
              const PublicKey<DCRTPoly>& pk,
              const std::vector<double>& M,
              int num_slots) {
    auto pt = cc->MakeCKKSPackedPlaintext(M);
    auto ct = cc->Encrypt(pk, pt);
    // SetSlots to s*d*d automatically replicates the data s times
    ct->SetSlots(num_slots);
    return ct;
}

// ---------------- Adapters ----------------
struct IAlgo {
    virtual ~IAlgo() = default;
    virtual std::string getName() const = 0;
    virtual Ciphertext<DCRTPoly> eval_inverse(const Ciphertext<DCRTPoly>& M) = 0;
    virtual void setSecretKey(const PrivateKey<DCRTPoly>& sk) = 0;
};

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
    //     auto ctM = encryptMatrix(ctx.cc, ctx.pk, M_set[0], ctx.num_slots);
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
        auto ctM = encryptMatrix(ctx.cc, ctx.pk, M_set[t], ctx.num_slots);

        auto t0 = std::chrono::steady_clock::now();
        auto ctMinv = algo->eval_inverse(ctM);
        auto t1 = std::chrono::steady_clock::now();

        Plaintext pt;
        ctx.cc->Decrypt(ctx.sk, ctMinv, &pt);
        pt->SetLength(d * d);
        std::vector<double> dec = pt->GetRealPackedValue();

        auto acc = accuracy(dec, Minv_ref_set[t]);

        // Debug: print expected vs computed
        std::cout << "\n=== Comparison ===\n";
        std::cout << "Expected inverse (first 16): ";
        for (int i = 0; i < std::min(16, d*d); ++i) std::cout << Minv_ref_set[t][i] << " ";
        std::cout << "\nComputed inverse (first 16): ";
        for (int i = 0; i < std::min(16, (int)dec.size()); ++i) std::cout << dec[i] << " ";
        std::cout << "\nmax_abs_err: " << acc.max_abs_err << ", mse: " << acc.mse << "\n";

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
struct Row {
    int d;
#ifdef RUN_BASELINE
    Stats st_base;
    Acc ac_base;
#endif
#ifdef RUN_LAZY
    Stats st_lazy;
    Acc ac_lazy;
#endif
};

// Run baseline and lazy for a single dimension, with cleanup between
static Row run_dimension(int d, int trials, unsigned seed) {
    Row row;
    row.d = d;

#ifdef RUN_BASELINE
    std::cout << "==== d=" << d << " Baseline ====\n";
    auto [st_base, ac_base] = run_one_method(
        "Baseline   ", d, trials, seed,
        [](const CryptoContext<DCRTPoly>& cc,
           const PublicKey<DCRTPoly>& pk, int d_) {
            return std::make_unique<BaselineAdapter>(cc, pk, d_);
        },
        HYBRID, /*useLazy=*/false
    );
    row.st_base = st_base;
    row.ac_base = ac_base;
    HardReset();  // Clean up after baseline
#endif

#ifdef RUN_LAZY
    std::cout << "==== d=" << d << " Lazy ====\n";
    auto [st_lazy, ac_lazy] = run_one_method(
        "Lazy       ", d, trials, seed,
        [](const CryptoContext<DCRTPoly>& cc,
           const PublicKey<DCRTPoly>& pk, int d_) {
            return std::make_unique<LazyAdapter>(cc, pk, d_);
        },
        BATCHED, /*useLazy=*/true
    );
    row.st_lazy = st_lazy;
    row.ac_lazy = ac_lazy;
    HardReset();  // Clean up after lazy
#endif

    return row;
}

// ---------------- main ----------------
int main() {
    HardReset();
    std::vector<int> dims(std::begin(kDims), std::end(kDims));

    // Run dimension-by-dimension: baseline vs lazy for each d before moving to next
    std::vector<Row> rows;
    rows.reserve(dims.size());

    for (int d : dims) {
        rows.push_back(run_dimension(d, kTrials, kSeed));
    }

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

    std::cout << "\n==== Summary ====\n";
    for (const auto& row : rows) {
        csv << row.d;

#ifdef RUN_BASELINE
        csv << "," << row.st_base.mean_ms << "," << row.st_base.std_ms
            << "," << row.ac_base.max_abs_err << "," << row.ac_base.mse
            << "," << row.ac_base.identity_max_err << "," << row.ac_base.identity_mse;
#endif
#ifdef RUN_LAZY
        csv << "," << row.st_lazy.mean_ms << "," << row.st_lazy.std_ms
            << "," << row.ac_lazy.max_abs_err << "," << row.ac_lazy.mse
            << "," << row.ac_lazy.identity_max_err << "," << row.ac_lazy.identity_mse;
#endif

        std::cout << "d=" << row.d;
#ifdef RUN_BASELINE
        std::cout << " | base "  << std::setprecision(3) << row.st_base.mean_ms << " ms";
#endif
#ifdef RUN_LAZY
        std::cout << " | lazy "  << row.st_lazy.mean_ms << " ms";
#endif

#if defined(RUN_BASELINE) && defined(RUN_LAZY)
        {
            double speedup_lazy = row.st_base.mean_ms / std::max(1e-12, row.st_lazy.mean_ms);
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
