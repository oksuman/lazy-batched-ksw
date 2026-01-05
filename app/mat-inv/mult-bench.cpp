#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "openfhe.h"
#include "kcl25.h"
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

using namespace lbcrypto;

// ================= Memory management =================
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

// ================= Configuration =================
static const int kDims[] = {8};
static const int kThreads[] = {8};
static const int kMultDepths[] = {24, 32, 36, 40};
static const int kScalingFactors[] = {30, 34, 40, 45};
static const int kTrials = 1;
static const unsigned kSeed = 1337u;
static const char* kCSV = "mult_bench_results.csv";
// =================================================

struct Stats {
    double mean_ms{0.0};
    double std_ms{0.0};
};

static Stats meanStd(const std::vector<double>& v) {
    if (v.empty()) return {};
    double m = std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    double var = 0.0;
    for (double x : v) var += (x - m) * (x - m);
    var /= (double)v.size();
    return {m, std::sqrt(var)};
}

struct Accuracy {
    double max_abs_err{0.0};
    double mse{0.0};
};

static Accuracy computeAccuracy(const std::vector<double>& computed, const std::vector<double>& expected) {
    const size_t n = std::min(computed.size(), expected.size());
    double maxe = 0.0;
    long double sse = 0.0L;
    for (size_t i = 0; i < n; ++i) {
        double e = std::abs(computed[i] - expected[i]);
        if (e > maxe) maxe = e;
        sse += (long double)e * (long double)e;
    }
    return {maxe, (double)(sse / (long double)std::max<size_t>(1, n))};
}

// Plain matrix multiplication
static std::vector<double> matrixMultiplyPlain(const std::vector<double>& A,
                                                const std::vector<double>& B, int d) {
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

// Plain matrix transpose
static std::vector<double> matrixTransposePlain(const std::vector<double>& A, int d) {
    std::vector<double> T(d * d);
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            T[j * d + i] = A[i * d + j];
        }
    }
    return T;
}

// Generate random matrix
static std::vector<double> randomMatrix(int d, std::mt19937& gen) {
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    std::vector<double> M(d * d);
    for (int i = 0; i < d * d; ++i) {
        M[i] = dis(gen);
    }
    return M;
}

// Plain trace computation (sum of diagonal elements, replicated)
static double matrixTracePlain(const std::vector<double>& A, int d) {
    double trace = 0.0;
    for (int i = 0; i < d; ++i) {
        trace += A[i * d + i];
    }
    return trace;
}

struct BenchResult {
    int d;
    int threads;
    int multDepth;
    int scalingFactor;
    int ringDim;

    // Multiplication
    double mult_base_ms;
    double mult_lazy_ms;
    double mult_base_err;
    double mult_lazy_err;

    // Transpose
    double trans_base_ms;
    double trans_lazy_ms;
    double trans_base_err;
    double trans_lazy_err;

    // Trace
    double trace_base_ms;
    double trace_lazy_ms;
    double trace_base_err;
    double trace_lazy_err;
};

static BenchResult runBenchmark(int d, int numThreads, int multDepth, int scalingModSize) {
    BenchResult result;
    result.d = d;
    result.threads = numThreads;
    result.multDepth = multDepth;
    result.scalingFactor = scalingModSize;

    // Set OpenMP threads
    #ifdef _OPENMP
    omp_set_num_threads(numThreads);
    #endif

    int batchSize = d * d;

    std::mt19937 gen(kSeed);

    // Generate test matrices
    std::vector<std::vector<double>> A_set, B_set;
    for (int t = 0; t < kTrials; ++t) {
        A_set.push_back(randomMatrix(d, gen));
        B_set.push_back(randomMatrix(d, gen));
    }

    // ============== Baseline ==============
    {
        CCParams<CryptoContextCKKSRNS> P;
        P.SetMultiplicativeDepth(multDepth);
        P.SetScalingModSize(scalingModSize);
        P.SetFirstModSize(scalingModSize + 10);
        P.SetSecurityLevel(HEStd_128_classic);
        P.SetBatchSize(batchSize);
        P.SetKeySwitchTechnique(HYBRID);

        auto cc = GenCryptoContext(P);

        // Ensure ring dimension >= 2^16
        if (cc->GetRingDimension() < (1 << 16)) {
            P.SetRingDim(1 << 16);
            cc = GenCryptoContext(P);
        }

        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);

        // Record and calculate batching params from actual ring dimension
        result.ringDim = static_cast<int>(cc->GetRingDimension());
        int max_batch = result.ringDim / 2;
        int s = std::min(max_batch / d / d, d);
        int num_slots = s * d * d;

        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);

        // Generate rotation keys
        MATINV_KCL25 impl(cc, kp.publicKey, d, 1, multDepth);
        RotationKeyCollector rk;
        rk.begin(num_slots, false);
        impl.eval_inverse_plan(rk);
        rk.generate(cc, kp.secretKey);

        // Encrypt test data
        auto encryptMatrix = [&](const std::vector<double>& M) {
            auto pt = cc->MakeCKKSPackedPlaintext(M);
            auto ct = cc->Encrypt(kp.publicKey, pt);
            ct->SetSlots(num_slots);
            return ct;
        };

        // Warm-up
        {
            auto ctA = encryptMatrix(A_set[0]);
            auto ctB = encryptMatrix(B_set[0]);
            auto ctC = impl.eval_mult(ctA, ctB);
            auto ctT = impl.eval_transpose(ctA);
        }

        // Benchmark multiplication
        std::vector<double> mult_times;
        double mult_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);
            auto ctB = encryptMatrix(B_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctC = impl.eval_mult(ctA, ctB);
            auto t1 = std::chrono::steady_clock::now();

            mult_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctC, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            auto expected = matrixMultiplyPlain(A_set[t], B_set[t], d);
            auto acc = computeAccuracy(computed, expected);
            if (acc.max_abs_err > mult_max_err) mult_max_err = acc.max_abs_err;
        }
        result.mult_base_ms = meanStd(mult_times).mean_ms;
        result.mult_base_err = mult_max_err;

        // Benchmark transpose
        std::vector<double> trans_times;
        double trans_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctT = impl.eval_transpose(ctA);
            auto t1 = std::chrono::steady_clock::now();

            trans_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctT, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            auto expected = matrixTransposePlain(A_set[t], d);
            auto acc = computeAccuracy(computed, expected);
            if (acc.max_abs_err > trans_max_err) trans_max_err = acc.max_abs_err;
        }
        result.trans_base_ms = meanStd(trans_times).mean_ms;
        result.trans_base_err = trans_max_err;

        // Benchmark trace
        std::vector<double> trace_times;
        double trace_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctTrace = impl.eval_trace(ctA, batchSize);
            auto t1 = std::chrono::steady_clock::now();

            trace_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctTrace, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            double expected = matrixTracePlain(A_set[t], d);
            double err = std::abs(computed[0] - expected);
            if (err > trace_max_err) trace_max_err = err;
        }
        result.trace_base_ms = meanStd(trace_times).mean_ms;
        result.trace_base_err = trace_max_err;

        cc->ClearEvalAutomorphismKeys();
        cc->ClearEvalMultKeys();
    }
    HardReset();

    // ============== Lazy ==============
    {
        CCParams<CryptoContextCKKSRNS> P;
        P.SetMultiplicativeDepth(multDepth);
        P.SetScalingModSize(scalingModSize);
        P.SetFirstModSize(scalingModSize + 10);
        P.SetSecurityLevel(HEStd_128_classic);
        P.SetBatchSize(batchSize);
        P.SetKeySwitchTechnique(BATCHED);

        auto cc = GenCryptoContext(P);

        // Ensure ring dimension >= 2^16
        if (cc->GetRingDimension() < (1 << 16)) {
            P.SetRingDim(1 << 16);
            cc = GenCryptoContext(P);
        }

        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);

        // Calculate batching params from actual ring dimension
        int max_batch = static_cast<int>(cc->GetRingDimension() / 2);
        int s = std::min(max_batch / d / d, d);
        int num_slots = s * d * d;

        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);

        // Generate lazy rotation keys
        MATINV_KCL25 impl(cc, kp.publicKey, d, 1, multDepth);
        RotationKeyCollectorLazy rk;
        rk.begin(num_slots);
        impl.eval_inverse_lazy_plan(rk);
        rk.generate(cc, kp.secretKey);

        // Encrypt test data
        auto encryptMatrix = [&](const std::vector<double>& M) {
            auto pt = cc->MakeCKKSPackedPlaintext(M);
            auto ct = cc->Encrypt(kp.publicKey, pt);
            ct->SetSlots(num_slots);
            return ct;
        };

        // Warm-up
        {
            auto ctA = encryptMatrix(A_set[0]);
            auto ctB = encryptMatrix(B_set[0]);
            auto ctC = impl.eval_mult_lazy(ctA, ctB);
            auto ctT = impl.eval_transpose_lazy(ctA);
        }

        // Benchmark multiplication
        std::vector<double> mult_times;
        double mult_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);
            auto ctB = encryptMatrix(B_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctC = impl.eval_mult_lazy(ctA, ctB);
            auto t1 = std::chrono::steady_clock::now();

            mult_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctC, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            auto expected = matrixMultiplyPlain(A_set[t], B_set[t], d);
            auto acc = computeAccuracy(computed, expected);
            if (acc.max_abs_err > mult_max_err) mult_max_err = acc.max_abs_err;
        }
        result.mult_lazy_ms = meanStd(mult_times).mean_ms;
        result.mult_lazy_err = mult_max_err;

        // Benchmark transpose
        std::vector<double> trans_times;
        double trans_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctT = impl.eval_transpose_lazy(ctA);
            auto t1 = std::chrono::steady_clock::now();

            trans_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctT, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            auto expected = matrixTransposePlain(A_set[t], d);
            auto acc = computeAccuracy(computed, expected);
            if (acc.max_abs_err > trans_max_err) trans_max_err = acc.max_abs_err;
        }
        result.trans_lazy_ms = meanStd(trans_times).mean_ms;
        result.trans_lazy_err = trans_max_err;

        // Benchmark trace
        std::vector<double> trace_times;
        double trace_max_err = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            auto ctA = encryptMatrix(A_set[t]);

            auto t0 = std::chrono::steady_clock::now();
            auto ctTrace = impl.eval_trace_lazy(ctA, batchSize);
            auto t1 = std::chrono::steady_clock::now();

            trace_times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());

            Plaintext pt;
            cc->Decrypt(kp.secretKey, ctTrace, &pt);
            pt->SetLength(d * d);
            auto computed = pt->GetRealPackedValue();
            double expected = matrixTracePlain(A_set[t], d);
            double err = std::abs(computed[0] - expected);
            if (err > trace_max_err) trace_max_err = err;
        }
        result.trace_lazy_ms = meanStd(trace_times).mean_ms;
        result.trace_lazy_err = trace_max_err;

        cc->ClearEvalAutomorphismKeys();
        cc->ClearEvalMultKeys();
    }
    HardReset();

    return result;
}

int main() {
    HardReset();  // Initial cleanup

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);
    csv << "d,threads,multDepth,scalingFactor,ringDim,"
        << "mult_base_ms,mult_lazy_ms,mult_speedup,"
        << "mult_base_err,mult_lazy_err,"
        << "trans_base_ms,trans_lazy_ms,trans_speedup,"
        << "trans_base_err,trans_lazy_err,"
        << "trace_base_ms,trace_lazy_ms,trace_speedup,"
        << "trace_base_err,trace_lazy_err\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "=== Matrix Mult / Transpose / Trace Benchmark ===\n\n";

    int total = sizeof(kDims)/sizeof(kDims[0]) *
                sizeof(kThreads)/sizeof(kThreads[0]) *
                sizeof(kMultDepths)/sizeof(kMultDepths[0]) *
                sizeof(kScalingFactors)/sizeof(kScalingFactors[0]);
    int current = 0;

    for (int d : kDims) {
        for (int threads : kThreads) {
            for (int depth : kMultDepths) {
                for (int sf : kScalingFactors) {
                    current++;
                    std::cout << "[" << current << "/" << total << "] "
                              << "d=" << d << ", threads=" << threads
                              << ", depth=" << depth << ", sf=" << sf << " ... ";
                    std::cout.flush();

                    try {
                        auto r = runBenchmark(d, threads, depth, sf);

                        double mult_speedup = r.mult_base_ms / std::max(1e-9, r.mult_lazy_ms);
                        double trans_speedup = r.trans_base_ms / std::max(1e-9, r.trans_lazy_ms);
                        double trace_speedup = r.trace_base_ms / std::max(1e-9, r.trace_lazy_ms);

                        std::cout << "ringDim=" << r.ringDim
                                  << ", mult: " << r.mult_base_ms << "/" << r.mult_lazy_ms
                                  << " ms (" << mult_speedup << "x), "
                                  << "trans: " << r.trans_base_ms << "/" << r.trans_lazy_ms
                                  << " ms (" << trans_speedup << "x), "
                                  << "trace: " << r.trace_base_ms << "/" << r.trace_lazy_ms
                                  << " ms (" << trace_speedup << "x)\n";

                        csv << r.d << "," << r.threads << "," << r.multDepth << "," << r.scalingFactor << "," << r.ringDim << ","
                            << r.mult_base_ms << "," << r.mult_lazy_ms << "," << mult_speedup << ","
                            << r.mult_base_err << "," << r.mult_lazy_err << ","
                            << r.trans_base_ms << "," << r.trans_lazy_ms << "," << trans_speedup << ","
                            << r.trans_base_err << "," << r.trans_lazy_err << ","
                            << r.trace_base_ms << "," << r.trace_lazy_ms << "," << trace_speedup << ","
                            << r.trace_base_err << "," << r.trace_lazy_err << "\n";
                        csv.flush();
                    } catch (const std::exception& e) {
                        std::cout << "FAILED: " << e.what() << "\n";
                        csv << d << "," << threads << "," << depth << "," << sf << ",0,"
                            << "ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,"
                            << "ERROR,ERROR,ERROR,ERROR,ERROR\n";
                    }
                }
            }
        }
    }

    csv.close();
    std::cout << "\nResults written to " << kCSV << "\n";
    return 0;
}
