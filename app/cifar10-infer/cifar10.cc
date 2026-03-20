/*
 * Portions of this file were adapted from the
 * "polycircuit" project (Apache License 2.0)
 * originally published by the Fair Math organization.
 *
 * Original source repository:
 *     https://github.com/fairmath/polycircuit
 *
 * The original code is licensed under the Apache License, Version 2.0.
 * You must comply with the Apache License 2.0 when using/modifying
 * or redistributing the original or derived portions.
 *
 * See the original LICENSE in that repository for details:
 *     https://github.com/fairmath/polycircuit/blob/main/LICENSE
 *
 * This file is included in this project under the terms of the MIT License.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <unistd.h>
#include <sys/wait.h>
#include <omp.h>

#include "memory_tracker.h"
#include "polycircuit/component/CIFAR10ImageClassification/CIFAR10ImageClassification.hpp"

#include <openfhe/pke/constants-defs.h>
#include <openfhe/pke/cryptocontext.h>
#include <openfhe/pke/openfhe.h>
#include <openfhe/pke/scheme/ckksrns/gen-cryptocontext-ckksrns-params.h>

struct ParamSet {
    std::string name;
    uint32_t ringDim;
    uint32_t multDepth;
    uint32_t scalingBits;
    uint32_t firstModBits;
};

static std::vector<ParamSet> get_param_sets() {
    return {
        {"N=2^14", 1u << 14, 7, 34, 46},
        {"N=2^15", 1u << 15, 13, 40, 51},
        {"N=2^16", 1u << 16, 24, 45, 56}
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

enum class Method { Baseline, Lazy };

struct ChildResult {
    double mean_ms;
    double std_ms;
    int correct_count;
    int total_tests;
    double max_abs_err;
    double mse;
    double peak_rss_mb;
};

static void configure_openmp(int num_threads) {
    setenv("OMP_PROC_BIND", "true", 1);
    setenv("OMP_PLACES", "cores", 1);
    setenv("OMP_DYNAMIC", "false", 1);
    setenv("OMP_NUM_THREADS", std::to_string(num_threads).c_str(), 1);
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
}

static ChildResult run_in_subprocess(
    const ParamSet& ps,
    const std::vector<std::pair<std::string, int>>& test_images,
    int repetitions, int num_threads,
    Method method,
    const std::vector<std::vector<double>>& baseline_logits)
{
    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); exit(1); }

    pid_t pid = fork();
    if (pid < 0) { perror("fork"); exit(1); }

    if (pid == 0) {
        // ---- Child process ----
        close(pipefd[0]);
        configure_openmp(num_threads);

        lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
        params.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        params.SetNumLargeDigits(3);
        params.SetMultiplicativeDepth(ps.multDepth);
        params.SetScalingModSize(ps.scalingBits);
        params.SetFirstModSize(ps.firstModBits);
        params.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
        params.SetKeySwitchTechnique(
            method == Method::Lazy ? lbcrypto::BATCHED : lbcrypto::HYBRID);
        params.SetBatchSize(4096);
        params.SetRingDim(ps.ringDim);

        auto cc = lbcrypto::GenCryptoContext(params);
        cc->Enable(lbcrypto::PKESchemeFeature::PKE);
        cc->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
        cc->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
        cc->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
        cc->Enable(lbcrypto::PKESchemeFeature::FHE);

        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);

        std::vector<int> rot_indices = {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600};
        if (method == Method::Lazy)
            cc->EvalLazyRotateKeyGen(keys.secretKey, rot_indices);
        else
            cc->EvalRotateKeyGen(keys.secretKey, rot_indices);

        // Warm-up
        {
            std::ifstream ifs(test_images[0].first);
            std::vector<double> img{std::istream_iterator<double>{ifs},
                                    std::istream_iterator<double>{}};
            auto ptxt = cc->MakeCKKSPackedPlaintext(img);
            auto c = cc->Encrypt(ptxt, keys.publicKey);
            if (method == Method::Lazy) {
                auto r = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                    polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc, std::move(c)).evaluate_lazy());
                lbcrypto::Plaintext rptx;
                cc->Decrypt(r, keys.secretKey, &rptx);
            } else {
                auto r = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                    polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc, std::move(c)).evaluate());
                lbcrypto::Plaintext rptx;
                cc->Decrypt(r, keys.secretKey, &rptx);
            }
        }

        std::vector<double> times_ms;
        int correct_count = 0;
        int total_tests = 0;
        double max_abs_err_overall = 0.0;
        long double mse_sum = 0.0L;
        int logit_idx = 0;

        for (const auto& [filename, ground_truth] : test_images) {
            std::ifstream ifs(filename);
            if (!ifs.is_open()) {
                std::cerr << "Unable to read: " << filename << "\n";
                _exit(1);
            }
            std::vector<double> image_data{std::istream_iterator<double>{ifs},
                                           std::istream_iterator<double>{}};

            for (int rep = 0; rep < repetitions; rep++) {
                auto ptxt = cc->MakeCKKSPackedPlaintext(image_data);
                auto c = cc->Encrypt(ptxt, keys.publicKey);

                auto t0 = std::chrono::steady_clock::now();
                lbcrypto::Ciphertext<lbcrypto::DCRTPoly> result;
                if (method == Method::Lazy) {
                    result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                        polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(
                            cc, std::move(c)).evaluate_lazy());
                } else {
                    result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                        polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(
                            cc, std::move(c)).evaluate());
                }
                auto t1 = std::chrono::steady_clock::now();

                lbcrypto::Plaintext rptx;
                cc->Decrypt(result, keys.secretKey, &rptx);
                rptx->SetLength(10);
                auto classes = rptx->GetRealPackedValue();

                int predicted = static_cast<int>(std::distance(
                    classes.begin(), std::max_element(classes.begin(), classes.end())));
                if (predicted == ground_truth) correct_count++;
                total_tests++;

                // Logits comparison with baseline (only if baseline logits available)
                if (!baseline_logits.empty() && logit_idx < (int)baseline_logits.size()) {
                    const auto& ref = baseline_logits[logit_idx];
                    for (int k = 0; k < 10; k++) {
                        double e = std::abs(classes[k] - ref[k]);
                        if (e > max_abs_err_overall) max_abs_err_overall = e;
                        mse_sum += (long double)e * (long double)e;
                    }
                }
                logit_idx++;

                times_ms.push_back(
                    std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
        }

        Stats st = meanStd(times_ms);
        ChildResult res;
        res.mean_ms = st.mean_ms;
        res.std_ms = st.std_ms;
        res.correct_count = correct_count;
        res.total_tests = total_tests;
        res.max_abs_err = max_abs_err_overall;
        res.mse = total_tests > 0
            ? (double)(mse_sum / (long double)(total_tests * 10))
            : 0.0;
        res.peak_rss_mb = getPeakRSSMB();

        write(pipefd[1], &res, sizeof(res));
        close(pipefd[1]);
        _exit(0);
    }

    // ---- Parent process ----
    close(pipefd[1]);

    ChildResult res{};
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

// Collect baseline logits for cross-method comparison
static std::vector<std::vector<double>> collect_baseline_logits(
    const ParamSet& ps,
    const std::vector<std::pair<std::string, int>>& test_images,
    int repetitions, int num_threads)
{
    int pipefd[2];
    if (pipe(pipefd) != 0) { perror("pipe"); exit(1); }

    pid_t pid = fork();
    if (pid < 0) { perror("fork"); exit(1); }

    int total = static_cast<int>(test_images.size()) * repetitions;

    if (pid == 0) {
        close(pipefd[0]);
        configure_openmp(num_threads);

        lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params;
        params.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        params.SetNumLargeDigits(3);
        params.SetMultiplicativeDepth(ps.multDepth);
        params.SetScalingModSize(ps.scalingBits);
        params.SetFirstModSize(ps.firstModBits);
        params.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
        params.SetKeySwitchTechnique(lbcrypto::HYBRID);
        params.SetBatchSize(4096);
        params.SetRingDim(ps.ringDim);

        auto cc = lbcrypto::GenCryptoContext(params);
        cc->Enable(lbcrypto::PKESchemeFeature::PKE);
        cc->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
        cc->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
        cc->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
        cc->Enable(lbcrypto::PKESchemeFeature::FHE);

        auto keys = cc->KeyGen();
        cc->EvalMultKeyGen(keys.secretKey);
        cc->EvalRotateKeyGen(keys.secretKey,
            {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600});

        // Write all logits as flat doubles: total * 10
        for (const auto& [filename, gt] : test_images) {
            std::ifstream ifs(filename);
            std::vector<double> img{std::istream_iterator<double>{ifs},
                                    std::istream_iterator<double>{}};
            for (int rep = 0; rep < repetitions; rep++) {
                auto ptxt = cc->MakeCKKSPackedPlaintext(img);
                auto c = cc->Encrypt(ptxt, keys.publicKey);
                auto result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                    polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(
                        cc, std::move(c)).evaluate());
                lbcrypto::Plaintext rptx;
                cc->Decrypt(result, keys.secretKey, &rptx);
                rptx->SetLength(10);
                auto classes = rptx->GetRealPackedValue();
                double logits[10];
                for (int k = 0; k < 10; k++) logits[k] = classes[k];
                write(pipefd[1], logits, sizeof(logits));
            }
        }
        close(pipefd[1]);
        _exit(0);
    }

    close(pipefd[1]);

    std::vector<std::vector<double>> all_logits;
    all_logits.reserve(total);
    for (int i = 0; i < total; i++) {
        double logits[10];
        ssize_t n = read(pipefd[0], logits, sizeof(logits));
        if (n != sizeof(logits)) break;
        all_logits.push_back(std::vector<double>(logits, logits + 10));
    }
    close(pipefd[0]);

    int status;
    waitpid(pid, &status, 0);

    return all_logits;
}

int main(int argc, char* argv[]) try {
    int REPETITIONS = 25;
    int num_threads = 1;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-r" || arg == "--repetitions") && i + 1 < argc)
            REPETITIONS = std::stoi(argv[++i]);
        else if ((arg == "-t" || arg == "--threads") && i + 1 < argc)
            num_threads = std::stoi(argv[++i]);
        else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [-r repetitions] [-t threads]\n";
            return EXIT_SUCCESS;
        }
    }
    const std::string data_dir = "data/";
    const std::string output_file = "cifar10_bench_results.csv";

    const std::vector<std::pair<std::string, int>> test_images = {
        {data_dir + "class-1.txt", 1},
        {data_dir + "class-6.txt", 6},
        {data_dir + "class-7.txt", 7},
        {data_dir + "class-9.txt", 9}
    };

    auto param_sets = get_param_sets();

    std::cout << "============================================\n"
              << "CIFAR-10 Inference Benchmark\n"
              << "============================================\n"
              << "Images: " << test_images.size() << "\n"
              << "Repetitions per image: " << REPETITIONS << "\n"
              << "Threads: " << num_threads << "\n"
              << "Parameter sets: " << param_sets.size() << "\n"
              << "============================================\n\n";

    std::ofstream csv(output_file);
    csv << std::fixed << std::setprecision(6);
    csv << "preset,N,depth_L,scalingBits,firstModBits,dnum"
        << ",base_mean_ms,base_std_ms,base_accuracy,base_peak_rss_mb"
        << ",lazy_mean_ms,lazy_std_ms,lazy_accuracy,lazy_max_abs_err,lazy_mse,lazy_peak_rss_mb"
        << ",speedup_lazy_x\n";

    for (const auto& ps : param_sets) {
        std::cout << "\n=== " << ps.name << " | RingDim=" << ps.ringDim
                  << " | L=" << ps.multDepth
                  << " | scale=" << ps.scalingBits
                  << " | first=" << ps.firstModBits
                  << " | dnum=3 ===\n";

        // Collect baseline logits for accuracy comparison
        std::cout << "  Collecting baseline logits...\n";
        auto baseline_logits = collect_baseline_logits(
            ps, test_images, REPETITIONS, num_threads);

        // Baseline
        std::cout << "  [Baseline] ... " << std::flush;
        std::vector<std::vector<double>> empty_logits;
        auto base = run_in_subprocess(ps, test_images, REPETITIONS, num_threads,
                                       Method::Baseline, empty_logits);
        double base_acc = 100.0 * base.correct_count / std::max(1, base.total_tests);
        std::cout << std::fixed << std::setprecision(3) << base.mean_ms << " ms"
                  << " | acc " << std::setprecision(1) << base_acc << "%"
                  << " | rss " << base.peak_rss_mb << " MB\n";

        // Lazy
        std::cout << "  [Lazy]     ... " << std::flush;
        auto lazy = run_in_subprocess(ps, test_images, REPETITIONS, num_threads,
                                       Method::Lazy, baseline_logits);
        double lazy_acc = 100.0 * lazy.correct_count / std::max(1, lazy.total_tests);
        std::cout << std::fixed << std::setprecision(3) << lazy.mean_ms << " ms"
                  << " | acc " << std::setprecision(1) << lazy_acc << "%"
                  << " | err " << std::setprecision(6) << lazy.max_abs_err
                  << " | rss " << std::setprecision(1) << lazy.peak_rss_mb << " MB\n";

        double speedup = base.mean_ms / std::max(1e-12, lazy.mean_ms);
        std::cout << "  speedup(lazy) = " << std::setprecision(3) << speedup << "x\n";

        // CSV row
        csv << ps.name << "," << ps.ringDim << ","
            << ps.multDepth << "," << ps.scalingBits << ","
            << ps.firstModBits << "," << 3
            << "," << std::setprecision(3) << base.mean_ms
            << "," << base.std_ms
            << "," << std::setprecision(1) << base_acc
            << "," << base.peak_rss_mb
            << "," << std::setprecision(3) << lazy.mean_ms
            << "," << lazy.std_ms
            << "," << std::setprecision(1) << lazy_acc
            << "," << std::setprecision(6) << lazy.max_abs_err
            << "," << lazy.mse
            << "," << std::setprecision(1) << lazy.peak_rss_mb
            << "," << std::setprecision(3) << speedup
            << "\n";
        csv.flush();
    }
    csv.close();

    std::cout << "\nResults written to " << output_file << "\n";
    return EXIT_SUCCESS;
}
catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
}
catch (...) {
    std::cerr << "An unknown exception was thrown." << std::endl;
    return EXIT_FAILURE;
}
