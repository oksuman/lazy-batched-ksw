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

//  OMP_NUM_THREADS=4 ./cifar -r 5
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

#include "polycircuit/component/CIFAR10ImageClassification/CIFAR10ImageClassification.hpp"

#include <openfhe/pke/constants-defs.h>
#include <openfhe/pke/cryptocontext.h>
#include <openfhe/pke/openfhe.h>
#include <openfhe/pke/scheme/ckksrns/gen-cryptocontext-ckksrns-params.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct BenchmarkResult {
    std::string image_name;
    int ground_truth;
    int predicted;
    bool correct;
    double inference_time_ms;
};

struct ExperimentStats {
    int total_tests;
    int correct_count;
    double accuracy;
    double total_time_ms;
    double avg_time_ms;
    double min_time_ms;
    double max_time_ms;
    double std_dev_ms;
    std::vector<BenchmarkResult> results;
};

/*
minimim required CKKS params
{
  "indexes_for_rotation_key": [-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600],
  "mult_depth": 5,
  "ring_dimension": 8192,
  "scale_mod_size": 25,
  "first_mod_size": 30,
  "batch_size": 4096,
  "enable_bootstrapping": false,
  "levels_available_after_bootstrap": 0,
  "level_budget": [0,0]
}
*/

ExperimentStats calculate_stats(const std::vector<BenchmarkResult>& results) {
    ExperimentStats stats;
    stats.results = results;
    stats.total_tests = results.size();
    stats.correct_count = std::count_if(results.begin(), results.end(),
                                        [](const auto& r) { return r.correct; });
    stats.accuracy = 100.0 * stats.correct_count / stats.total_tests;

    std::vector<double> times;
    for (const auto& r : results) {
        times.push_back(r.inference_time_ms);
    }

    stats.total_time_ms = std::accumulate(times.begin(), times.end(), 0.0);
    stats.avg_time_ms = stats.total_time_ms / times.size();
    stats.min_time_ms = *std::min_element(times.begin(), times.end());
    stats.max_time_ms = *std::max_element(times.begin(), times.end());

    double variance = 0.0;
    for (double t : times) {
        variance += (t - stats.avg_time_ms) * (t - stats.avg_time_ms);
    }
    stats.std_dev_ms = std::sqrt(variance / times.size());

    return stats;
}

void save_results_to_file(const std::string& filename,
                          const ExperimentStats& exp1_stats,
                          const ExperimentStats& exp2_stats,
                          int repetitions_per_image) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Unable to create output file: " + filename);
    }

    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    ofs << "=========================================\n";
    ofs << "CIFAR-10 Inference Benchmark Results\n";
    ofs << "=========================================\n";
    ofs << "Date: " << std::ctime(&now_time);
    ofs << "Repetitions per image: " << repetitions_per_image << "\n";
    ofs << "Total images tested: 4 (class-1, class-6, class-7, class-9)\n";
    ofs << "Total tests per experiment: " << exp1_stats.total_tests << "\n\n";

    ofs << "=========================================\n";
    ofs << "EXPERIMENT 1: Standard evaluate()\n";
    ofs << "=========================================\n";
    ofs << "Total tests:      " << exp1_stats.total_tests << "\n";
    ofs << "Correct:          " << exp1_stats.correct_count << "\n";
    ofs << "Accuracy:         " << std::fixed << std::setprecision(2) << exp1_stats.accuracy << "%\n";
    ofs << "Total time:       " << std::setprecision(2) << exp1_stats.total_time_ms << " ms\n";
    ofs << "Average time:     " << std::setprecision(3) << exp1_stats.avg_time_ms << " ms\n";
    ofs << "Min time:         " << std::setprecision(3) << exp1_stats.min_time_ms << " ms\n";
    ofs << "Max time:         " << std::setprecision(3) << exp1_stats.max_time_ms << " ms\n";
    ofs << "Std deviation:    " << std::setprecision(3) << exp1_stats.std_dev_ms << " ms\n\n";

    ofs << "=========================================\n";
    ofs << "EXPERIMENT 2: Lazy evaluate_lazy()\n";
    ofs << "=========================================\n";
    ofs << "Total tests:      " << exp2_stats.total_tests << "\n";
    ofs << "Correct:          " << exp2_stats.correct_count << "\n";
    ofs << "Accuracy:         " << std::fixed << std::setprecision(2) << exp2_stats.accuracy << "%\n";
    ofs << "Total time:       " << std::setprecision(2) << exp2_stats.total_time_ms << " ms\n";
    ofs << "Average time:     " << std::setprecision(3) << exp2_stats.avg_time_ms << " ms\n";
    ofs << "Min time:         " << std::setprecision(3) << exp2_stats.min_time_ms << " ms\n";
    ofs << "Max time:         " << std::setprecision(3) << exp2_stats.max_time_ms << " ms\n";
    ofs << "Std deviation:    " << std::setprecision(3) << exp2_stats.std_dev_ms << " ms\n\n";

    ofs << "=========================================\n";
    ofs << "COMPARISON\n";
    ofs << "=========================================\n";
    double speedup = exp1_stats.avg_time_ms / exp2_stats.avg_time_ms;
    double time_diff = exp1_stats.avg_time_ms - exp2_stats.avg_time_ms;
    double percent_diff = 100.0 * time_diff / exp1_stats.avg_time_ms;

    ofs << "Speedup (Exp1/Exp2):    " << std::setprecision(3) << speedup << "x\n";
    ofs << "Time difference:        " << std::setprecision(3) << time_diff << " ms\n";
    ofs << "Percentage difference:  " << std::setprecision(2) << percent_diff << "%\n";
    ofs << "Winner:                 " << (exp2_stats.avg_time_ms < exp1_stats.avg_time_ms ? "Experiment 2 (Lazy)" : "Experiment 1 (Standard)") << "\n\n";

    ofs.close();
    std::cout << "\nResults saved to: " << filename << std::endl;
}

int main(int argc, char *argv[]) try
{
    po::options_description desc("Allowed parameters");
    desc.add_options()
    ("help,h", "produce help message")
    ("repetitions,r", po::value<int>()->default_value(25), "number of repetitions per image (default: 25)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << '\n';
        return EXIT_SUCCESS;
    }

    const int REPETITIONS = vm["repetitions"].as<int>();
    const std::string data_dir = "/home/user/shkim/lazy-batched-ksw/app/cifar10-infer/data/";
    const std::string output_file = "/home/user/shkim/lazy-batched-ksw/app/cifar10-infer/benchmark_results.txt";

    std::vector<std::pair<std::string, int>> test_images = {
        {data_dir + "class-1.txt", 1},
        {data_dir + "class-6.txt", 6},
        {data_dir + "class-7.txt", 7},
        {data_dir + "class-9.txt", 9}
    };

    std::vector<std::string> class_names = {
        "airplane", "automobile", "bird", "cat", "deer",
        "dog", "frog", "horse", "ship", "truck"
    };

    std::cout << "============================================" << std::endl;
    std::cout << "CIFAR-10 Inference Benchmark" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Images: " << test_images.size() << std::endl;
    std::cout << "Repetitions per image: " << REPETITIONS << std::endl;
    std::cout << "Total tests per experiment: " << (test_images.size() * REPETITIONS) << std::endl;
    std::cout << "============================================\n" << std::endl;

    // Security level setting
    lbcrypto::SecurityLevel securityLevel = lbcrypto::HEStd_128_classic;

    // ========== EXPERIMENT 1: Standard evaluate() ==========
    std::cout << "\n>>> EXPERIMENT 1: Standard evaluate() <<<\n" << std::endl;

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters1;
    parameters1.SetSecurityLevel(securityLevel);
    parameters1.SetNumLargeDigits(3);
    parameters1.SetMultiplicativeDepth(5);
    parameters1.SetScalingModSize(34);
    parameters1.SetFirstModSize(46);
    parameters1.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
    parameters1.SetKeySwitchTechnique(lbcrypto::HYBRID);
    parameters1.SetBatchSize(4096);
    // parameters1.SetRingDim(8192);
    parameters1.SetRingDim(16384);

    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc1 = lbcrypto::GenCryptoContext(parameters1);
    cc1->Enable(lbcrypto::PKESchemeFeature::PKE);
    cc1->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
    cc1->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
    cc1->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
    cc1->Enable(lbcrypto::PKESchemeFeature::FHE);

    auto keys1 = cc1->KeyGen();
    cc1->EvalMultKeyGen(keys1.secretKey);
    cc1->EvalRotateKeyGen(keys1.secretKey, {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600});

    std::cout << "Context setup complete (Ring dimension: " << cc1->GetRingDimension() << ")" << std::endl;

    std::vector<BenchmarkResult> exp1_results;

    for (const auto& [filename, ground_truth] : test_images) {
        std::cout << "\nTesting " << filename << " (" << REPETITIONS << " repetitions)..." << std::endl;

        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            throw std::runtime_error("Unable to read: " + filename);
        }
        std::vector<double> image_data{std::istream_iterator<double>{ifs}, std::istream_iterator<double>{}};
        ifs.close();

        for (int rep = 0; rep < REPETITIONS; rep++) {
            lbcrypto::Plaintext ptxt = cc1->MakeCKKSPackedPlaintext(image_data);
            auto c = cc1->Encrypt(ptxt, keys1.publicKey);

            auto start = std::chrono::high_resolution_clock::now();
            auto result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc1, std::move(c)).evaluate()
            );
            auto end = std::chrono::high_resolution_clock::now();

            lbcrypto::Plaintext rptx;
            cc1->Decrypt(result, keys1.secretKey, &rptx);
            rptx->SetLength(10);
            auto classes = rptx->GetRealPackedValue();
            int predicted = std::distance(classes.begin(), std::max_element(classes.begin(), classes.end()));

            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

            exp1_results.push_back({
                filename,
                ground_truth,
                predicted,
                predicted == ground_truth,
                time_ms
            });

            if ((rep + 1) % 5 == 0) {
                std::cout << "  Progress: " << (rep + 1) << "/" << REPETITIONS << std::endl;
            }
        }
    }

    auto exp1_stats = calculate_stats(exp1_results);
    std::cout << "\n--- Experiment 1 Summary ---" << std::endl;
    std::cout << "Accuracy: " << std::fixed << std::setprecision(2) << exp1_stats.accuracy << "%" << std::endl;
    std::cout << "Avg time: " << std::setprecision(3) << exp1_stats.avg_time_ms << " ms" << std::endl;

    // ========== EXPERIMENT 2: Lazy evaluate_lazy() ==========
    std::cout << "\n\n>>> EXPERIMENT 2: Lazy evaluate_lazy() <<<\n" << std::endl;

    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters2;
    parameters2.SetSecurityLevel(securityLevel);
    parameters2.SetNumLargeDigits(3);
    parameters2.SetMultiplicativeDepth(5);
    parameters2.SetScalingModSize(34);
    parameters2.SetFirstModSize(46);
    parameters2.SetScalingTechnique(lbcrypto::ScalingTechnique::FLEXIBLEAUTO);
    parameters2.SetKeySwitchTechnique(lbcrypto::BATCHED);
    parameters2.SetBatchSize(4096);
    // parameters2.SetRingDim(8192);
    parameters2.SetRingDim(16384);

    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc2 = lbcrypto::GenCryptoContext(parameters2);
    cc2->Enable(lbcrypto::PKESchemeFeature::PKE);
    cc2->Enable(lbcrypto::PKESchemeFeature::KEYSWITCH);
    cc2->Enable(lbcrypto::PKESchemeFeature::LEVELEDSHE);
    cc2->Enable(lbcrypto::PKESchemeFeature::ADVANCEDSHE);
    cc2->Enable(lbcrypto::PKESchemeFeature::FHE);

    auto keys2 = cc2->KeyGen();
    cc2->EvalMultKeyGen(keys2.secretKey);
    cc2->EvalLazyRotateKeyGen(keys2.secretKey, {-3, -2, -1, 10, 20, 40, 50, 100, 200, 400, 800, 1600});

    std::cout << "Context setup complete (Ring dimension: " << cc2->GetRingDimension() << ")" << std::endl;

    std::vector<BenchmarkResult> exp2_results;

    for (const auto& [filename, ground_truth] : test_images) {
        std::cout << "\nTesting " << filename << " (" << REPETITIONS << " repetitions)..." << std::endl;

        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            throw std::runtime_error("Unable to read: " + filename);
        }
        std::vector<double> image_data{std::istream_iterator<double>{ifs}, std::istream_iterator<double>{}};
        ifs.close();

        for (int rep = 0; rep < REPETITIONS; rep++) {
            lbcrypto::Plaintext ptxt = cc2->MakeCKKSPackedPlaintext(image_data);
            auto c = cc2->Encrypt(ptxt, keys2.publicKey);

            auto start = std::chrono::high_resolution_clock::now();
            auto result = std::get<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>(
                polycircuit::CIFAR10ImageClassification<lbcrypto::DCRTPoly>(cc2, std::move(c)).evaluate_lazy()
            );
            auto end = std::chrono::high_resolution_clock::now();

            lbcrypto::Plaintext rptx;
            cc2->Decrypt(result, keys2.secretKey, &rptx);
            rptx->SetLength(10);
            auto classes = rptx->GetRealPackedValue();
            int predicted = std::distance(classes.begin(), std::max_element(classes.begin(), classes.end()));

            double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

            exp2_results.push_back({
                filename,
                ground_truth,
                predicted,
                predicted == ground_truth,
                time_ms
            });

            if ((rep + 1) % 5 == 0) {
                std::cout << "  Progress: " << (rep + 1) << "/" << REPETITIONS << std::endl;
            }
        }
    }

    auto exp2_stats = calculate_stats(exp2_results);
    std::cout << "\n--- Experiment 2 Summary ---" << std::endl;
    std::cout << "Accuracy: " << std::fixed << std::setprecision(2) << exp2_stats.accuracy << "%" << std::endl;
    std::cout << "Avg time: " << std::setprecision(3) << exp2_stats.avg_time_ms << " ms" << std::endl;

    // ========== Save Results ==========
    save_results_to_file(output_file, exp1_stats, exp2_stats, REPETITIONS);

    // ========== Final Comparison ==========
    std::cout << "\n============================================" << std::endl;
    std::cout << "FINAL COMPARISON" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "Experiment 1 (Standard):  " << std::setprecision(3) << exp1_stats.avg_time_ms << " ms avg" << std::endl;
    std::cout << "Experiment 2 (Lazy):      " << std::setprecision(3) << exp2_stats.avg_time_ms << " ms avg" << std::endl;
    std::cout << "Speedup:                  " << std::setprecision(3) << (exp1_stats.avg_time_ms / exp2_stats.avg_time_ms) << "x" << std::endl;
    std::cout << "============================================" << std::endl;

    return EXIT_SUCCESS;
}
catch (const po::error& ex)
{
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
}
catch (const std::exception& ex)
{
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
}
catch (...)
{
    std::cerr << "An unknown exception was thrown." << std::endl;
    return EXIT_FAILURE;
}

