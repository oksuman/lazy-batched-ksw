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
#include <omp.h>

#include "polycircuit/component/CIFAR10ImageClassification/CIFAR10ImageClassification.hpp"

#include <openfhe/pke/constants-defs.h>
#include <openfhe/pke/cryptocontext.h>
#include <openfhe/pke/openfhe.h>
#include <openfhe/pke/scheme/ckksrns/gen-cryptocontext-ckksrns-params.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct ParamSet {
    std::string name;
    uint32_t ringDim;
    uint32_t multDepth;
    uint32_t scalingBits;
    uint32_t firstModBits;
};

std::vector<ParamSet> get_param_sets() {
    return {
        {"N=2^13", 1u << 13, 5, 20, 30},
        {"N=2^14", 1u << 14, 7, 34, 46},
        {"N=2^15", 1u << 15, 13, 40, 51},
        {"N=2^16", 1u << 16, 24, 45, 56}
    };
}

struct BenchmarkResult {
    std::string param_name;
    int threads;
    std::string method;
    std::string image_name;
    int ground_truth;
    int predicted;
    bool correct;
    double inference_time_ms;
};

struct ExperimentStats {
    std::string param_name;
    int threads;
    std::string method;
    int total_tests;
    int correct_count;
    double accuracy;
    double total_time_ms;
    double avg_time_ms;
    double min_time_ms;
    double max_time_ms;
    double std_dev_ms;
};

ExperimentStats calculate_stats(const std::vector<BenchmarkResult>& results) {
    ExperimentStats stats{};
    if (results.empty()) {
        return stats;
    }

    stats.param_name = results[0].param_name;
    stats.threads = results[0].threads;
    stats.method = results[0].method;
    stats.total_tests = static_cast<int>(results.size());
    stats.correct_count = static_cast<int>(std::count_if(results.begin(), results.end(),
                                                        [](const auto& r) { return r.correct; }));
    stats.accuracy = 100.0 * static_cast<double>(stats.correct_count) / static_cast<double>(stats.total_tests);

    std::vector<double> times;
    times.reserve(results.size());
    for (const auto& r : results) {
        times.push_back(r.inference_time_ms);
    }

    stats.total_time_ms = std::accumulate(times.begin(), times.end(), 0.0);
    stats.avg_time_ms = stats.total_time_ms / static_cast<double>(times.size());
    stats.min_time_ms = *std::min_element(times.begin(), times.end());
    stats.max_time_ms = *std::max_element(times.begin(), times.end());

    double variance = 0.0;
    for (double t : times) {
        double d = t - stats.avg_time_ms;
        variance += d * d;
    }
    stats.std_dev_ms = std::sqrt(variance / static_cast<double>(times.size()));

    return stats;
}

void save_results_to_csv(const std::string& filename,
                         const std::vector<ExperimentStats>& all_stats) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Unable to create output file: " + filename);
    }

    ofs << "param_set,threads,method,total_tests,correct,accuracy,"
        << "total_time_ms,avg_time_ms,min_time_ms,max_time_ms,std_dev_ms\n";

    for (const auto& stats : all_stats) {
        ofs << stats.param_name << ","
            << stats.threads << ","
            << stats.method << ","
            << stats.total_tests << ","
            << stats.correct_count << ","
            << std::fixed << std::setprecision(2) << stats.accuracy << ","
            << std::setprecision(3) << stats.total_time_ms << ","
            << stats.avg_time_ms << ","
            << stats.min_time_ms << ","
            << stats.max_time_ms << ","
            << stats.std_dev_ms << "\n";
    }

    ofs.close();
    std::cout << "\nResults saved to: " << filename << std::endl;
}

static void configure_openmp(int num_threads) {
    setenv("OMP_PROC_BIND", "true", 1);
    setenv("OMP_PLACES", "cores", 1);
    setenv("OMP_DYNAMIC", "false", 1);
    setenv("OMP_NUM_THREADS", std::to_string(num_threads).c_str(), 1);

    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);

    int actual = 0;
#pragma omp parallel
    {
#pragma omp single
        { actual = omp_get_num_threads(); }
    }

    std::cout << "OpenMP requested=" << num_threads
              << ", omp_get_max_threads()=" << omp_get_max_threads()
              << ", actual_in_parallel=" << actual
              << std::endl;
}

int main(int argc, char* argv[]) try {
    po::options_description desc("Allowed parameters");
    desc.add_options()
        ("help,h", "produce help message")
        ("repetitions,r", po::value<int>()->default_value(5), "number of repetitions per image (default: 5)")
        ("threads,t", po::value<std::string>()->default_value("1,32"), "comma-separated thread counts (default: 1,32)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << '\n';
        return EXIT_SUCCESS;
    }

    const int REPETITIONS = vm["repetitions"].as<int>();
    const std::string data_dir = "/home/user/shkim/lazy-batched-ksw/app/cifar10-infer/data/";
    const std::string output_file = "/home/user/shkim/lazy-batched-ksw/app/cifar10-infer/cifar10_benchmark_results.csv";

    std::vector<int> thread_counts;
    {
        std::string threads_str = vm["threads"].as<std::string>();
        std::stringstream ss(threads_str);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (!item.empty()) {
                thread_counts.push_back(std::stoi(item));
            }
        }
        if (thread_counts.empty()) {
            throw std::runtime_error("No valid thread counts were provided.");
        }
    }

    const std::vector<std::pair<std::string, int>> test_images = {
        {data_dir + "class-1.txt", 1},
        {data_dir + "class-6.txt", 6},
        {data_dir + "class-7.txt", 7},
        {data_dir + "class-9.txt", 9}
    };

    auto param_sets = get_param_sets();

    std::cout << "============================================" << std::endl;
    std::cout << "CIFAR-10 Inference Benchmark" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Images: " << test_images.size() << std::endl;
    std::cout << "Repetitions per image: " << REPETITIONS << std::endl;
    std::cout << "Total tests per experiment: " << (test_images.size() * REPETITIONS) << std::endl;
    std::cout << "Parameter sets: " << param_sets.size() << std::endl;
    std::cout << "Thread counts: ";
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        std::cout << thread_counts[i];
        if (i < thread_counts.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    std::cout << "============================================\n" << std::endl;

    lbcrypto::SecurityLevel securityLevel = lbcrypto::HEStd_128_classic;
    std::vector<ExperimentStats> all_stats;

    for (const auto& ps : param_sets) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Parameter Set: " << ps.name << std::endl;
        std::cout << "  RingDim: " << ps.ringDim << std::endl;
        std::cout << "  MultDepth: " << ps.multDepth << std::endl;
        std::cout << "  ScalingBits: " << ps.scalingBits << std::endl;
        std::cout << "  FirstModBits: " << ps.firstModBits << std::endl;
        std::cout << "========================================" << std::endl;

        for (int num_threads : thread_counts) {
            std::cout << "\n>>> Thread Count: " << num_threads << " <<<" << std::endl;

            configure_openmp(num_threads);

            std::cout << "\n  [Baseline - HYBRID]" << std::endl;

            lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params_baseline;
            params_baseline.SetSecurityLevel(securityLevel);
            params_baseline.SetNumLargeDigits(3);
            params_baseline.SetMultiplicativeDepth(ps.multDepth);
            params_baseline.SetScalingModSize(ps.scalingBits);
            params_baseline.SetFirstModSize(ps.firstModBits);
            params_baseline.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
            params_baseline.SetKeySwitchTechnique(lbcrypto::HYBRID);
            params_baseline.SetBatchSize(4096);
            params_baseline.SetRingDim(ps.ringDim);

            auto cc_baseline = lbcrypto::GenCryptoContext(params_baseline);
            cc_baseline->Enable(lbcrypto::PKESchemeFeature::PKE);
            cc_baseline->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
            cc_baseline->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
            cc_baseline->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
            cc_baseline->Enable(lbcrypto::PKESchemeFeature::FHE);

            auto keys_baseline = cc_baseline->KeyGen();
            cc_baseline->EvalMultKeyGen(keys_baseline.secretKey);
            cc_baseline->EvalRotateKeyGen(keys_baseline.secretKey,
                                          {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600});

            std::vector<BenchmarkResult> baseline_results;
            baseline_results.reserve(test_images.size() * static_cast<size_t>(REPETITIONS));

            for (const auto& [filename, ground_truth] : test_images) {
                std::ifstream ifs(filename);
                if (!ifs.is_open()) {
                    throw std::runtime_error("Unable to read: " + filename);
                }
                std::vector<double> image_data{std::istream_iterator<double>{ifs},
                                               std::istream_iterator<double>{}};
                ifs.close();

                for (int rep = 0; rep < REPETITIONS; rep++) {
                    lbcrypto::Plaintext ptxt = cc_baseline->MakeCKKSPackedPlaintext(image_data);
                    auto c = cc_baseline->Encrypt(ptxt, keys_baseline.publicKey);

                    auto start = std::chrono::high_resolution_clock::now();
                    auto result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                        polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc_baseline, std::move(c)).evaluate()
                    );
                    auto end = std::chrono::high_resolution_clock::now();

                    lbcrypto::Plaintext rptx;
                    cc_baseline->Decrypt(result, keys_baseline.secretKey, &rptx);
                    rptx->SetLength(10);
                    auto classes = rptx->GetRealPackedValue();
                    int predicted = static_cast<int>(std::distance(classes.begin(),
                                                                   std::max_element(classes.begin(), classes.end())));

                    double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

                    baseline_results.push_back({
                        ps.name,
                        num_threads,
                        "Baseline",
                        filename,
                        ground_truth,
                        predicted,
                        predicted == ground_truth,
                        time_ms
                    });
                }
            }

            auto baseline_stats = calculate_stats(baseline_results);
            all_stats.push_back(baseline_stats);
            std::cout << "  Baseline: " << std::fixed << std::setprecision(3)
                      << baseline_stats.avg_time_ms << " ms avg, "
                      << std::setprecision(2) << baseline_stats.accuracy << "% accuracy" << std::endl;

            std::cout << "\n  [Lazy - BATCHED]" << std::endl;

            lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> params_lazy;
            params_lazy.SetSecurityLevel(securityLevel);
            params_lazy.SetNumLargeDigits(3);
            params_lazy.SetMultiplicativeDepth(ps.multDepth);
            params_lazy.SetScalingModSize(ps.scalingBits);
            params_lazy.SetFirstModSize(ps.firstModBits);
            params_lazy.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
            params_lazy.SetKeySwitchTechnique(lbcrypto::BATCHED);
            params_lazy.SetBatchSize(4096);
            params_lazy.SetRingDim(ps.ringDim);

            auto cc_lazy = lbcrypto::GenCryptoContext(params_lazy);
            cc_lazy->Enable(lbcrypto::PKESchemeFeature::PKE);
            cc_lazy->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
            cc_lazy->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
            cc_lazy->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
            cc_lazy->Enable(lbcrypto::PKESchemeFeature::FHE);

            auto keys_lazy = cc_lazy->KeyGen();
            cc_lazy->EvalMultKeyGen(keys_lazy.secretKey);
            cc_lazy->EvalLazyRotateKeyGen(keys_lazy.secretKey,
                                          {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600});

            std::vector<BenchmarkResult> lazy_results;
            lazy_results.reserve(test_images.size() * static_cast<size_t>(REPETITIONS));

            for (const auto& [filename, ground_truth] : test_images) {
                std::ifstream ifs(filename);
                if (!ifs.is_open()) {
                    throw std::runtime_error("Unable to read: " + filename);
                }
                std::vector<double> image_data{std::istream_iterator<double>{ifs},
                                               std::istream_iterator<double>{}};
                ifs.close();

                for (int rep = 0; rep < REPETITIONS; rep++) {
                    lbcrypto::Plaintext ptxt = cc_lazy->MakeCKKSPackedPlaintext(image_data);
                    auto c = cc_lazy->Encrypt(ptxt, keys_lazy.publicKey);

                    auto start = std::chrono::high_resolution_clock::now();
                    auto result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                        polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc_lazy, std::move(c)).evaluate_lazy()
                    );
                    auto end = std::chrono::high_resolution_clock::now();

                    lbcrypto::Plaintext rptx;
                    cc_lazy->Decrypt(result, keys_lazy.secretKey, &rptx);
                    rptx->SetLength(10);
                    auto classes = rptx->GetRealPackedValue();
                    int predicted = static_cast<int>(std::distance(classes.begin(),
                                                                   std::max_element(classes.begin(), classes.end())));

                    double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

                    lazy_results.push_back({
                        ps.name,
                        num_threads,
                        "Lazy",
                        filename,
                        ground_truth,
                        predicted,
                        predicted == ground_truth,
                        time_ms
                    });
                }
            }

            auto lazy_stats = calculate_stats(lazy_results);
            all_stats.push_back(lazy_stats);
            std::cout << "  Lazy:     " << std::fixed << std::setprecision(3)
                      << lazy_stats.avg_time_ms << " ms avg, "
                      << std::setprecision(2) << lazy_stats.accuracy << "% accuracy" << std::endl;

            double speedup = baseline_stats.avg_time_ms / std::max(1e-9, lazy_stats.avg_time_ms);
            std::cout << "  Speedup:  " << std::setprecision(3) << speedup << "x" << std::endl;
        }
    }

    save_results_to_csv(output_file, all_stats);

    std::cout << "\n\n============================================" << std::endl;
    std::cout << "SUMMARY TABLE" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::left << std::setw(20) << "Param Set"
              << std::setw(10) << "Threads"
              << std::setw(12) << "Method"
              << std::setw(15) << "Avg Time (ms)"
              << std::setw(12) << "Accuracy (%)" << std::endl;
    std::cout << std::string(69, '-') << std::endl;

    for (const auto& stats : all_stats) {
        std::cout << std::left << std::setw(20) << stats.param_name
                  << std::setw(10) << stats.threads
                  << std::setw(12) << stats.method
                  << std::fixed << std::setprecision(3) << std::setw(15) << stats.avg_time_ms
                  << std::setprecision(2) << std::setw(12) << stats.accuracy << std::endl;
    }
    std::cout << "============================================" << std::endl;

    return EXIT_SUCCESS;
}
catch (const po::error& ex) {
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
}
catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
}
catch (...) {
    std::cerr << "An unknown exception was thrown." << std::endl;
    return EXIT_FAILURE;
}
