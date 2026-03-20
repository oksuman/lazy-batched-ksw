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
#include <cstdlib>
#include <unistd.h>
#include <sys/wait.h>

#include "openfhe.h"
#include "jkls18.h"
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"
#include "memory_tracker.h"

using namespace lbcrypto;

// -------- Toggle experiments --------
#define RUN_BASELINE 1
#define RUN_HOIST 1
#define RUN_DOUBLE_HOIST 1
#define RUN_LAZY 1
// ------------------------------------

// ================= User-tunable =================
static const int      kTrials = 1;
static const unsigned kSeed   = 1337u;
static const int      kDims[] = {8};
static const char*    kCSV    = "jkls18_bench_results.csv";
// ================================================

// -------- Parameter presets --------
struct Preset {
    std::string name;
    uint32_t ringDim;
    uint32_t multDepth;
    uint32_t scalingBits;
    uint32_t firstModBits;
    SecurityLevel sec;
};

static std::vector<Preset> MakePresets() {
    return {
        {"N=2^14", 1<<14, 7, 34, 46, HEStd_128_classic},
        // {"N=2^15", 1<<15, 13, 40, 51, HEStd_128_classic},
        // {"N=2^16", 1<<16, 24, 45, 56, HEStd_128_classic}
    };
}

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
};

static ContextPack makeContextForD(const Preset& ps, int d, KeySwitchTechnique ksTech) {
    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(ps.multDepth);
    P.SetRingDim(ps.ringDim);
    P.SetScalingModSize(ps.scalingBits);
    P.SetFirstModSize(ps.firstModBits);
    P.SetNumLargeDigits(3);
    P.SetSecurityLevel(ps.sec);
    P.SetBatchSize(d * d);
    P.SetKeySwitchTechnique(ksTech);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    return {cc, kp.publicKey, kp.secretKey};
}

// ---------------- Encryption helper ----------------
static Ciphertext<DCRTPoly>
encryptMatrix(const CryptoContext<DCRTPoly>& cc,
              const PublicKey<DCRTPoly>& pk,
              const std::vector<double>& M) {
    auto pt = cc->MakeCKKSPackedPlaintext(M);
    return cc->Encrypt(pk, pt);
}

// ---------------- Method enum ----------------
enum class Method { Baseline, Lazy, Hoist, DoubleHoist };

// Result communicated from child to parent via pipe
struct BenchResult {
    double mean_ms;
    double std_ms;
    double max_abs_err;
    double mse;
    double peak_rss_mb;
};

// Run a single method benchmark in a forked child process for fair peak RSS measurement.
static BenchResult run_in_subprocess(
    const Preset& ps, int d, int trials, unsigned seed,
    Method method)
{
    int pipefd[2];
    if (pipe(pipefd) != 0) {
        perror("pipe");
        exit(1);
    }

    pid_t pid = fork();
    if (pid < 0) {
        perror("fork");
        exit(1);
    }

    if (pid == 0) {
        // ---- Child process ----
        close(pipefd[0]);

        KeySwitchTechnique ksTech = (method == Method::Lazy) ? BATCHED : HYBRID;
        ContextPack ctx = makeContextForD(ps, d, ksTech);

        // Generate rotation keys (same index set for all methods)
        MATMULT_JKLS18 planner(ctx.cc, ctx.pk, d);
        {
            int actualSlots = static_cast<int>(ctx.cc->GetRingDimension() / 2);
            RotationKeyCollector rk;
            rk.begin(actualSlots, method == Method::Lazy);
            planner.eval_mult_hoist_plan(rk);
            rk.generate(ctx.cc, ctx.sk);
        }

        // Create algo
        MATMULT_JKLS18 impl(ctx.cc, ctx.pk, d);

        // Prepare test data
        std::mt19937 gen(seed);
        std::vector<std::vector<double>> A_set, B_set, C_ref_set;
        for (int t = 0; t < trials; ++t) {
            auto A = randomMatrix(d, gen);
            auto B = randomMatrix(d, gen);
            auto C = matmulPlain(A, B, d);
            A_set.emplace_back(std::move(A));
            B_set.emplace_back(std::move(B));
            C_ref_set.emplace_back(std::move(C));
        }

        // Warm-up
        {
            auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[0]);
            auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[0]);
            Ciphertext<DCRTPoly> ctC;
            switch (method) {
                case Method::Baseline:    ctC = impl.eval_mult(ctA, ctB); break;
                case Method::Hoist:       ctC = impl.eval_mult_hoist(ctA, ctB); break;
                case Method::DoubleHoist: ctC = impl.eval_mult_double_hoist(ctA, ctB); break;
                case Method::Lazy:        ctC = impl.eval_mult_lazy(ctA, ctB); break;
            }
            Plaintext pt;
            ctx.cc->Decrypt(ctx.sk, ctC, &pt);
        }

        // Timed runs
        std::vector<double> times_ms;
        double max_abs_err_overall = 0.0;
        long double mse_sum = 0.0L;

        for (int t = 0; t < trials; ++t) {
            auto ctA = encryptMatrix(ctx.cc, ctx.pk, A_set[t]);
            auto ctB = encryptMatrix(ctx.cc, ctx.pk, B_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            Ciphertext<DCRTPoly> ctC;
            switch (method) {
                case Method::Baseline:    ctC = impl.eval_mult(ctA, ctB); break;
                case Method::Hoist:       ctC = impl.eval_mult_hoist(ctA, ctB); break;
                case Method::DoubleHoist: ctC = impl.eval_mult_double_hoist(ctA, ctB); break;
                case Method::Lazy:        ctC = impl.eval_mult_lazy(ctA, ctB); break;
            }
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
        BenchResult res;
        res.mean_ms = st.mean_ms;
        res.std_ms = st.std_ms;
        res.max_abs_err = max_abs_err_overall;
        res.mse = (double)(mse_sum / trials);
        res.peak_rss_mb = getPeakRSSMB();

        write(pipefd[1], &res, sizeof(res));
        close(pipefd[1]);
        _exit(0);
    }

    // ---- Parent process ----
    close(pipefd[1]);

    BenchResult res{};
    ssize_t n = read(pipefd[0], &res, sizeof(res));
    close(pipefd[0]);

    int status;
    waitpid(pid, &status, 0);

    if (n != sizeof(res) || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        std::cerr << "  ERROR: child process failed\n";
        return {};
    }

    return res;
}

static const char* methodName(Method m) {
    switch (m) {
        case Method::Baseline:    return "Baseline";
        case Method::Hoist:       return "Hoist";
        case Method::DoubleHoist: return "DoubleHoist";
        case Method::Lazy:        return "Lazy";
    }
    return "Unknown";
}

// ---------------- Result row ----------------
struct Row { int d; Stats st; Acc ac; double peak_rss_mb{0.0}; };

static Row run_method_for_dim(const Preset& ps, int d, int trials, unsigned seed, Method method) {
    std::cout << "    d=" << d << " [" << methodName(method) << "] ... " << std::flush;
    auto res = run_in_subprocess(ps, d, trials, seed, method);
    std::cout << std::fixed << std::setprecision(3) << res.mean_ms << " ms"
              << " | err " << std::setprecision(6) << res.max_abs_err
              << " | rss " << std::setprecision(1) << res.peak_rss_mb << " MB\n";
    return {d, {res.mean_ms, res.std_ms}, {res.max_abs_err, res.mse}, res.peak_rss_mb};
}

// ---------------- main ----------------
int main() {
    const char* omp_threads = std::getenv("OMP_NUM_THREADS");
    std::cout << "OMP_NUM_THREADS = " << (omp_threads ? omp_threads : "not set (using default)") << "\n\n";

    std::vector<int> dims(std::begin(kDims), std::end(kDims));
    auto presets = MakePresets();

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);

    // CSV header
    csv << "preset,N,depth_L,scalingBits,firstModBits,dnum,d";
#ifdef RUN_BASELINE
    csv << ",base_mean_ms,base_std_ms,base_max_abs_err,base_mse,base_peak_rss_mb";
#endif
#ifdef RUN_HOIST
    csv << ",hoist_mean_ms,hoist_std_ms,hoist_max_abs_err,hoist_mse,hoist_peak_rss_mb";
#endif
#ifdef RUN_DOUBLE_HOIST
    csv << ",dhoist_mean_ms,dhoist_std_ms,dhoist_max_abs_err,dhoist_mse,dhoist_peak_rss_mb";
#endif
#ifdef RUN_LAZY
    csv << ",lazy_mean_ms,lazy_std_ms,lazy_max_abs_err,lazy_mse,lazy_peak_rss_mb";
#endif
#if defined(RUN_BASELINE) && defined(RUN_LAZY)
    csv << ",speedup_lazy_x";
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
    csv << ",speedup_hoist_x";
#endif
#if defined(RUN_BASELINE) && defined(RUN_DOUBLE_HOIST)
    csv << ",speedup_dhoist_x";
#endif
    csv << "\n";

    for (const auto& ps : presets) {
        std::cout << "\n=== " << ps.name << " | RingDim=" << ps.ringDim
                  << " | L=" << ps.multDepth
                  << " | scale=" << ps.scalingBits
                  << " | first=" << ps.firstModBits
                  << " | dnum=3 ===\n";

        for (size_t i = 0; i < dims.size(); ++i) {
            const int d = dims[i];

#ifdef RUN_BASELINE
            auto base = run_method_for_dim(ps, d, kTrials, kSeed, Method::Baseline);
#endif
#ifdef RUN_HOIST
            auto hoist = run_method_for_dim(ps, d, kTrials, kSeed, Method::Hoist);
#endif
#ifdef RUN_DOUBLE_HOIST
            auto dhoist = run_method_for_dim(ps, d, kTrials, kSeed, Method::DoubleHoist);
#endif
#ifdef RUN_LAZY
            auto lazy = run_method_for_dim(ps, d, kTrials, kSeed, Method::Lazy);
#endif

            // CSV row
            csv << ps.name << "," << ps.ringDim << ","
                << ps.multDepth << "," << ps.scalingBits << "," << ps.firstModBits << ","
                << 3 << "," << d;

#ifdef RUN_BASELINE
            csv << "," << base.st.mean_ms << "," << base.st.std_ms
                << "," << base.ac.max_abs_err << "," << base.ac.mse
                << "," << std::setprecision(1) << base.peak_rss_mb << std::setprecision(6);
#endif
#ifdef RUN_HOIST
            csv << "," << hoist.st.mean_ms << "," << hoist.st.std_ms
                << "," << hoist.ac.max_abs_err << "," << hoist.ac.mse
                << "," << std::setprecision(1) << hoist.peak_rss_mb << std::setprecision(6);
#endif
#ifdef RUN_DOUBLE_HOIST
            csv << "," << dhoist.st.mean_ms << "," << dhoist.st.std_ms
                << "," << dhoist.ac.max_abs_err << "," << dhoist.ac.mse
                << "," << std::setprecision(1) << dhoist.peak_rss_mb << std::setprecision(6);
#endif
#ifdef RUN_LAZY
            csv << "," << lazy.st.mean_ms << "," << lazy.st.std_ms
                << "," << lazy.ac.max_abs_err << "," << lazy.ac.mse
                << "," << std::setprecision(1) << lazy.peak_rss_mb << std::setprecision(6);
#endif

            // Console summary
            std::cout << "  --- d=" << d << " summary ---\n";
#ifdef RUN_BASELINE
            std::cout << "    Baseline    : " << std::setprecision(3) << base.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << base.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_HOIST
            std::cout << "    Hoist       : " << std::setprecision(3) << hoist.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << hoist.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_DOUBLE_HOIST
            std::cout << "    DoubleHoist : " << std::setprecision(3) << dhoist.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << dhoist.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_LAZY
            std::cout << "    Lazy        : " << std::setprecision(3) << lazy.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << lazy.peak_rss_mb << " MB\n";
#endif

#if defined(RUN_BASELINE) && defined(RUN_LAZY)
            double speedup_lazy = base.st.mean_ms / std::max(1e-12, lazy.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_lazy << std::setprecision(6);
            std::cout << "    speedup(lazy)   = " << std::setprecision(3) << speedup_lazy << "x\n";
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
            double speedup_hoist = base.st.mean_ms / std::max(1e-12, hoist.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_hoist << std::setprecision(6);
            std::cout << "    speedup(hoist)  = " << std::setprecision(3) << speedup_hoist << "x\n";
#endif
#if defined(RUN_BASELINE) && defined(RUN_DOUBLE_HOIST)
            double speedup_dhoist = base.st.mean_ms / std::max(1e-12, dhoist.st.mean_ms);
            csv << "," << std::setprecision(3) << speedup_dhoist << std::setprecision(6);
            std::cout << "    speedup(dhoist) = " << std::setprecision(3) << speedup_dhoist << "x\n";
#endif
            csv << "\n";
        }
        csv.flush();
    }
    csv.close();

    std::cout << "\nResults written to " << kCSV << "\n";
    return 0;
}
