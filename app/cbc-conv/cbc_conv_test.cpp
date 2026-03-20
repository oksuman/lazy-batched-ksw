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
#include "cbc_conv.h"
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"
#include "memory_tracker.h"

using namespace lbcrypto;

// -------- Toggle experiments --------
#define RUN_BASELINE  1
#define RUN_HOIST     1
#define RUN_LAZY      1
#define RUN_TWOSTAGE  1
// ------------------------------------

// ================= User-tunable =================
static const int      kTrials = 1;
static const unsigned kSeed   = 1337u;
static const char*    kCSV    = "cbc_conv_bench_results.csv";
// ================================================

// -------- Convolution configurations --------
struct ConvConfig {
    std::string name;
    int c_in;
    int c_out;
    int W;
    int gap;
};

static std::vector<ConvConfig> MakeConfigs() {
    return {
        {"A", 16, 16, 32, 1},
        // {"B", 32, 32, 16, 1},
        // {"C", 64, 64,  8, 1},
    };
}

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

// -------- Timing/accuracy utils --------
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

// -------- Random data generation --------
struct ConvData {
    std::vector<std::vector<double>> input_channels;
    std::vector<std::vector<std::vector<std::vector<double>>>> weights;
    std::vector<double> bias;
};

static ConvData generateRandomData(const ConvConfig& cfg, std::mt19937& gen) {
    ConvData data;
    int numSlots = cfg.W * cfg.W;

    std::uniform_real_distribution<double> input_dist(0.0, 1.0);
    std::uniform_real_distribution<double> weight_dist(-1.0, 1.0);

    data.input_channels.resize(cfg.c_in);
    for (int ch = 0; ch < cfg.c_in; ch++) {
        data.input_channels[ch].resize(numSlots);
        for (double& x : data.input_channels[ch]) x = input_dist(gen);
    }

    data.weights.resize(cfg.c_out);
    for (int i = 0; i < cfg.c_out; i++) {
        data.weights[i].resize(cfg.c_in);
        for (int ch = 0; ch < cfg.c_in; ch++) {
            data.weights[i][ch].resize(3);
            for (int j = 0; j < 3; j++) {
                data.weights[i][ch][j].resize(3);
                for (int k = 0; k < 3; k++)
                    data.weights[i][ch][j][k] = weight_dist(gen);
            }
        }
    }

    data.bias.resize(cfg.c_out, 0.0);
    return data;
}

// -------- Method enum --------
enum class Method { Baseline, Hoist, Lazy, TwoStage };

// -------- Result communicated via pipe --------
struct BenchResult {
    double mean_ms;
    double std_ms;
    double max_abs_err;
    double mse;
    double peak_rss_mb;
};

static BenchResult run_in_subprocess(
    const Preset& ps, const ConvConfig& cfg, int trials, unsigned seed,
    Method method)
{
    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); exit(1); }

    pid_t pid = fork();
    if (pid < 0) { perror("fork"); exit(1); }

    if (pid == 0) {
        // ---- Child process ----
        close(pipefd[0]);

        int batchSize = cfg.W * cfg.W;
        KeySwitchTechnique ksTech = (method == Method::Lazy || method == Method::TwoStage) ? BATCHED : HYBRID;

        CCParams<CryptoContextCKKSRNS> P;
        P.SetMultiplicativeDepth(ps.multDepth);
        P.SetRingDim(ps.ringDim);
        P.SetScalingModSize(ps.scalingBits);
        P.SetFirstModSize(ps.firstModBits);
        P.SetNumLargeDigits(3);
        P.SetSecurityLevel(ps.sec);
        P.SetBatchSize(batchSize);
        P.SetKeySwitchTechnique(ksTech);

        auto cc = GenCryptoContext(P);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);

        auto kp = cc->KeyGen();
        cc->EvalMultKeyGen(kp.secretKey);

        CBC_Conv conv(cc, kp.publicKey, cfg.W, cfg.W, cfg.c_in, cfg.c_out, cfg.gap);

        // Generate rotation keys
        {
            int actualSlots = static_cast<int>(cc->GetRingDimension() / 2);
            bool useLazy = (method == Method::Lazy || method == Method::TwoStage);
            RotationKeyCollector rk;
            rk.begin(actualSlots, useLazy);
            if (method == Method::TwoStage)
                conv.eval_twostage_plan(rk);
            else
                conv.eval_hoisted_plan(rk);
            rk.generate(cc, kp.secretKey);
        }

        // Generate random data
        std::mt19937 gen(seed);
        ConvData data = generateRandomData(cfg, gen);
        conv.encodeKernels(data.weights, data.bias);

        // Encrypt input channels
        std::vector<Ciphertext<DCRTPoly>> input_cts(cfg.c_in);
        for (int ch = 0; ch < cfg.c_in; ch++) {
            auto pt = cc->MakeCKKSPackedPlaintext(data.input_channels[ch]);
            input_cts[ch] = cc->Encrypt(kp.publicKey, pt);
        }

        // Plaintext reference
        auto ref = conv.eval_plain(data.input_channels, data.weights, data.bias);

        // Warm-up
        {
            std::vector<Ciphertext<DCRTPoly>> warmup;
            if (method == Method::TwoStage)
                warmup = conv.eval_twostage(input_cts);
            else if (method == Method::Lazy)
                warmup = conv.eval_lazy(input_cts);
            else if (method == Method::Hoist)
                warmup = conv.eval_hoisted(input_cts);
            else
                warmup = conv.eval_baseline(input_cts);
            Plaintext pt;
            cc->Decrypt(kp.secretKey, warmup[0], &pt);
        }

        // Timed trials
        std::vector<double> times_ms; times_ms.reserve(trials);
        double max_abs_err_overall = 0.0; long double mse_sum = 0.0L;

        for (int t = 0; t < trials; ++t) {
            auto t0 = std::chrono::steady_clock::now();
            std::vector<Ciphertext<DCRTPoly>> results;
            if (method == Method::TwoStage)
                results = conv.eval_twostage(input_cts);
            else if (method == Method::Lazy)
                results = conv.eval_lazy(input_cts);
            else if (method == Method::Hoist)
                results = conv.eval_hoisted(input_cts);
            else
                results = conv.eval_baseline(input_cts);
            auto t1 = std::chrono::steady_clock::now();

            // Accuracy: check all output channels
            for (int out = 0; out < cfg.c_out; out++) {
                Plaintext pt;
                cc->Decrypt(kp.secretKey, results[out], &pt);
                pt->SetLength(batchSize);
                std::vector<double> dec = pt->GetRealPackedValue();
                auto acc = accuracy(dec, ref[out]);
                if (acc.max_abs_err > max_abs_err_overall) max_abs_err_overall = acc.max_abs_err;
                mse_sum += acc.mse;
            }

            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            times_ms.push_back(ms);
        }

        Stats st = meanStd(times_ms);
        BenchResult res;
        res.mean_ms = st.mean_ms;
        res.std_ms = st.std_ms;
        res.max_abs_err = max_abs_err_overall;
        res.mse = (double)(mse_sum / (trials * cfg.c_out));
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
        case Method::Baseline:  return "Baseline";
        case Method::Hoist:     return "Hoist";
        case Method::Lazy:      return "Lazy";
        case Method::TwoStage:  return "TwoStage";
    }
    return "Unknown";
}

// -------- Result row --------
struct Row { std::string config; Stats st; Acc ac; double peak_rss_mb{0.0}; };

static Row run_method_for_config(const Preset& ps, const ConvConfig& cfg,
                                  int trials, unsigned seed, Method method) {
    std::cout << "    Config " << cfg.name
              << " (c_in=" << cfg.c_in << ", c_out=" << cfg.c_out
              << ", W=" << cfg.W << ", gap=" << cfg.gap
              << ") [" << methodName(method) << "] ... " << std::flush;
    auto res = run_in_subprocess(ps, cfg, trials, seed, method);
    std::cout << std::fixed << std::setprecision(3) << res.mean_ms << " ms"
              << " | err " << std::setprecision(6) << res.max_abs_err
              << " | rss " << std::setprecision(1) << res.peak_rss_mb << " MB\n";
    return {cfg.name, {res.mean_ms, res.std_ms}, {res.max_abs_err, res.mse}, res.peak_rss_mb};
}

// -------- main --------
int main() {
    const char* omp_threads = std::getenv("OMP_NUM_THREADS");
    std::cout << "OMP_NUM_THREADS = "
              << (omp_threads ? omp_threads : "not set (using default)") << "\n\n";

    auto configs = MakeConfigs();
    auto presets = MakePresets();

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);

    // CSV header
    csv << "preset,N,depth_L,scalingBits,firstModBits,dnum"
        << ",config,c_in,c_out,W,gap";
#ifdef RUN_BASELINE
    csv << ",base_mean_ms,base_std_ms,base_max_abs_err,base_mse,base_peak_rss_mb";
#endif
#ifdef RUN_HOIST
    csv << ",hoist_mean_ms,hoist_std_ms,hoist_max_abs_err,hoist_mse,hoist_peak_rss_mb";
#endif
#ifdef RUN_LAZY
    csv << ",lazy_mean_ms,lazy_std_ms,lazy_max_abs_err,lazy_mse,lazy_peak_rss_mb";
#endif
#ifdef RUN_TWOSTAGE
    csv << ",twostage_mean_ms,twostage_std_ms,twostage_max_abs_err,twostage_mse,twostage_peak_rss_mb";
#endif
#if defined(RUN_BASELINE) && defined(RUN_LAZY)
    csv << ",speedup_lazy_x";
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
    csv << ",speedup_hoist_x";
#endif
#if defined(RUN_BASELINE) && defined(RUN_TWOSTAGE)
    csv << ",speedup_twostage_x";
#endif
    csv << "\n";

    for (const auto& ps : presets) {
        std::cout << "\n=== " << ps.name << " | RingDim=" << ps.ringDim
                  << " | L=" << ps.multDepth
                  << " | scale=" << ps.scalingBits
                  << " | first=" << ps.firstModBits
                  << " | dnum=3 ===\n";

        for (size_t i = 0; i < configs.size(); ++i) {
            const auto& cfg = configs[i];

#ifdef RUN_BASELINE
            auto base = run_method_for_config(ps, cfg, kTrials, kSeed, Method::Baseline);
#endif
#ifdef RUN_HOIST
            auto hoist = run_method_for_config(ps, cfg, kTrials, kSeed, Method::Hoist);
#endif
#ifdef RUN_LAZY
            auto lazy = run_method_for_config(ps, cfg, kTrials, kSeed, Method::Lazy);
#endif
#ifdef RUN_TWOSTAGE
            auto twostage = run_method_for_config(ps, cfg, kTrials, kSeed, Method::TwoStage);
#endif

            // CSV row
            csv << ps.name << "," << ps.ringDim << ","
                << ps.multDepth << "," << ps.scalingBits << ","
                << ps.firstModBits << "," << 3
                << "," << cfg.name << "," << cfg.c_in << ","
                << cfg.c_out << "," << cfg.W << "," << cfg.gap;

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
#ifdef RUN_LAZY
            csv << "," << lazy.st.mean_ms << "," << lazy.st.std_ms
                << "," << lazy.ac.max_abs_err << "," << lazy.ac.mse
                << "," << std::setprecision(1) << lazy.peak_rss_mb << std::setprecision(6);
#endif
#ifdef RUN_TWOSTAGE
            csv << "," << twostage.st.mean_ms << "," << twostage.st.std_ms
                << "," << twostage.ac.max_abs_err << "," << twostage.ac.mse
                << "," << std::setprecision(1) << twostage.peak_rss_mb << std::setprecision(6);
#endif

            // Console summary
            std::cout << "  --- Config " << cfg.name << " summary ---\n";
#ifdef RUN_BASELINE
            std::cout << "    Baseline : " << std::setprecision(3) << base.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << base.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_HOIST
            std::cout << "    Hoist    : " << std::setprecision(3) << hoist.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << hoist.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_LAZY
            std::cout << "    Lazy     : " << std::setprecision(3) << lazy.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << lazy.peak_rss_mb << " MB\n";
#endif
#ifdef RUN_TWOSTAGE
            std::cout << "    TwoStage : " << std::setprecision(3) << twostage.st.mean_ms
                      << " ms | rss " << std::setprecision(1) << twostage.peak_rss_mb << " MB\n";
#endif

#if defined(RUN_BASELINE) && defined(RUN_LAZY)
            double sp_lazy = base.st.mean_ms / std::max(1e-12, lazy.st.mean_ms);
            csv << "," << std::setprecision(3) << sp_lazy << std::setprecision(6);
            std::cout << "    speedup(lazy)  = " << std::setprecision(3) << sp_lazy << "x\n";
#endif
#if defined(RUN_BASELINE) && defined(RUN_HOIST)
            double sp_hoist = base.st.mean_ms / std::max(1e-12, hoist.st.mean_ms);
            csv << "," << std::setprecision(3) << sp_hoist << std::setprecision(6);
            std::cout << "    speedup(hoist) = " << std::setprecision(3) << sp_hoist << "x\n";
#endif
#if defined(RUN_BASELINE) && defined(RUN_TWOSTAGE)
            double sp_twostage = base.st.mean_ms / std::max(1e-12, twostage.st.mean_ms);
            csv << "," << std::setprecision(3) << sp_twostage << std::setprecision(6);
            std::cout << "    speedup(twostage) = " << std::setprecision(3) << sp_twostage << "x\n";
#endif
            csv << "\n";
        }
        csv.flush();
    }
    csv.close();

    std::cout << "\nResults written to " << kCSV << "\n";
    return 0;
}
