// Task 1: Unit-level performance microbenchmarks
// 1a. EvalDirectRotate vs EvalRotate (single rotation comparison)
// 1b. Column shifting: baseline, hoisting, double-hoisting (placeholder), lazy-batched
// 1c. mod Q vs mod PQ arithmetic cost (ptMult, ctAdd)

#include "openfhe.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace lbcrypto;

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
        {"N=2^15", 1<<15, 13, 40, 51, HEStd_128_classic},
        // {"N=2^16", 1<<16, 24, 45, 56, HEStd_128_classic},
    };
}

// -------- Timing utils --------
struct Stats { double mean_ms{0.0}; double std_ms{0.0}; };
static Stats meanStd(const std::vector<double>& v) {
    if (v.empty()) return {};
    double m = std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    double var = 0.0; for (double x : v) var += (x - m)*(x - m);
    var /= (double)v.size();
    return {m, std::sqrt(var)};
}

static std::vector<double> makeMsg(int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; i++) v[i] = 0.001 * (i + 1);
    return v;
}

// ============================================================
// 1a. EvalDirectRotate vs EvalRotate (single rotation)
// ============================================================
static void bench_1a(const Preset& ps, int trials, std::ofstream& csv) {
    std::cout << "  [1a] EvalDirectRotate vs EvalRotate ...\n";

    std::vector<int32_t> indices = {1, 5, 13};

    for (int32_t idx : indices) {
        // EvalRotate (standard HYBRID)
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(HYBRID);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalRotateKeyGen(kp.secretKey, {idx});

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);

            // Warm-up
            for (int i = 0; i < 3; i++) { auto r = cc->EvalRotate(ct, idx); }

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                auto t0 = std::chrono::steady_clock::now();
                auto r = cc->EvalRotate(ct, idx);
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1a_single_rotation," << idx << ",EvalRotate,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    idx=" << idx << " EvalRotate: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }

        // EvalDirectRotate (BATCHED)
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(BATCHED);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalLazyRotateKeyGen(kp.secretKey, {idx});

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);

            // Warm-up
            for (int i = 0; i < 3; i++) { auto r = cc->EvalDirectRotate(ct, idx); }

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                auto t0 = std::chrono::steady_clock::now();
                auto r = cc->EvalDirectRotate(ct, idx);
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1a_single_rotation," << idx << ",EvalDirectRotate,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    idx=" << idx << " EvalDirectRotate: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }
    }
}

// ============================================================
// 1b. Column shifting (rotate-accumulate pattern, k rotations of same ct)
//     Methods: baseline, hoisting, hoisted-direct, lazy-batched
// ============================================================
static void bench_1b(const Preset& ps, int trials, std::ofstream& csv) {
    std::cout << "  [1b] Rotate-accumulate pattern ...\n";

    std::vector<int> ks = {4, 8, 16};

    for (int k : ks) {
        std::vector<int32_t> indices(k);
        for (int i = 0; i < k; i++) indices[i] = i + 1;

        // Method 1: Baseline (individual EvalRotate)
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(HYBRID);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalRotateKeyGen(kp.secretKey, indices);

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);
            auto zero_pt = cc->MakeCKKSPackedPlaintext(std::vector<double>(128, 0.0));
            auto acc = cc->Encrypt(kp.publicKey, zero_pt);

            // Warm-up
            for (int i = 0; i < k; i++) cc->EvalAddInPlace(acc, cc->EvalRotate(ct, i + 1));

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                acc = cc->Encrypt(kp.publicKey, zero_pt);
                auto t0 = std::chrono::steady_clock::now();
                for (int i = 0; i < k; i++)
                    cc->EvalAddInPlace(acc, cc->EvalRotate(ct, i + 1));
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1b_rotate_accum," << k << ",Baseline,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    k=" << k << " Baseline: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }

        // Method 2: Hoisting (EvalFastRotation)
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(HYBRID);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalRotateKeyGen(kp.secretKey, indices);

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);
            auto zero_pt = cc->MakeCKKSPackedPlaintext(std::vector<double>(128, 0.0));
            auto acc = cc->Encrypt(kp.publicKey, zero_pt);
            usint M = 2 * cc->GetRingDimension();

            // Warm-up
            {
                auto digits = cc->EvalFastRotationPrecompute(ct);
                for (int i = 0; i < k; i++)
                    cc->EvalAddInPlace(acc, cc->EvalFastRotation(ct, i + 1, M, digits));
            }

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                acc = cc->Encrypt(kp.publicKey, zero_pt);
                auto t0 = std::chrono::steady_clock::now();
                auto digits = cc->EvalFastRotationPrecompute(ct);
                for (int i = 0; i < k; i++)
                    cc->EvalAddInPlace(acc, cc->EvalFastRotation(ct, i + 1, M, digits));
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1b_rotate_accum," << k << ",Hoisting,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    k=" << k << " Hoisting: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }

        // Method 3: Hoisted DirectRotate
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(BATCHED);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalLazyRotateKeyGen(kp.secretKey, indices);

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);
            auto zero_pt = cc->MakeCKKSPackedPlaintext(std::vector<double>(128, 0.0));
            auto acc = cc->Encrypt(kp.publicKey, zero_pt);

            // Warm-up
            {
                auto digits = cc->EvalDirectRotatePrecompute(ct);
                for (int i = 0; i < k; i++)
                    cc->EvalAddInPlace(acc, cc->EvalDirectRotate(ct, i + 1, digits));
            }

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                acc = cc->Encrypt(kp.publicKey, zero_pt);
                auto t0 = std::chrono::steady_clock::now();
                auto digits = cc->EvalDirectRotatePrecompute(ct);
                for (int i = 0; i < k; i++)
                    cc->EvalAddInPlace(acc, cc->EvalDirectRotate(ct, i + 1, digits));
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1b_rotate_accum," << k << ",HoistedDirect,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    k=" << k << " HoistedDirect: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }

        // Method 4: Lazy-batched (EvalLazyRotate + EvalBatchedKS)
        {
            CCParams<CryptoContextCKKSRNS> P;
            P.SetMultiplicativeDepth(ps.multDepth);
            P.SetRingDim(ps.ringDim);
            P.SetScalingModSize(ps.scalingBits);
            P.SetFirstModSize(ps.firstModBits);
            P.SetNumLargeDigits(3);
            P.SetSecurityLevel(ps.sec);
            P.SetBatchSize(128);
            P.SetKeySwitchTechnique(BATCHED);

            auto cc = GenCryptoContext(P);
            cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
            auto kp = cc->KeyGen();
            cc->EvalLazyRotateKeyGen(kp.secretKey, indices);

            auto msg = makeMsg(128);
            auto pt = cc->MakeCKKSPackedPlaintext(msg);
            auto ct = cc->Encrypt(kp.publicKey, pt);

            // Warm-up
            {
                auto acc = ct->Clone();
                for (int i = 0; i < k; i++)
                    cc->EvalLazyAddInPlace(acc, cc->EvalLazyRotate(ct, i + 1));
                cc->EvalBatchedKS(acc);
            }

            std::vector<double> times;
            for (int t = 0; t < trials; t++) {
                auto t0 = std::chrono::steady_clock::now();
                auto acc = ct->Clone();
                for (int i = 0; i < k; i++)
                    cc->EvalLazyAddInPlace(acc, cc->EvalLazyRotate(ct, i + 1));
                auto result = cc->EvalBatchedKS(acc);
                auto t1 = std::chrono::steady_clock::now();
                times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            auto st = meanStd(times);
            csv << ps.name << ",1b_rotate_accum," << k << ",LazyBatched,"
                << st.mean_ms << "," << st.std_ms << "\n";
            std::cout << "    k=" << k << " LazyBatched: " << std::fixed << std::setprecision(3)
                      << st.mean_ms << " ms\n";
        }
    }
}

// ============================================================
// 1c. mod Q vs mod PQ arithmetic cost
// ============================================================
static void bench_1c(const Preset& ps, int trials, std::ofstream& csv) {
    std::cout << "  [1c] mod Q vs mod PQ arithmetic ...\n";

    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(ps.multDepth);
    P.SetRingDim(ps.ringDim);
    P.SetScalingModSize(ps.scalingBits);
    P.SetFirstModSize(ps.firstModBits);
    P.SetNumLargeDigits(3);
    P.SetSecurityLevel(ps.sec);
    P.SetBatchSize(128);
    P.SetKeySwitchTechnique(HYBRID);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);
    auto kp = cc->KeyGen();
    cc->EvalRotateKeyGen(kp.secretKey, {1});

    auto msg = makeMsg(128);
    auto pt = cc->MakeCKKSPackedPlaintext(msg);
    auto ct = cc->Encrypt(kp.publicKey, pt);
    auto ct2 = cc->Encrypt(kp.publicKey, pt);

    // Extend ct to PQ basis for comparison
    auto ctExt = cc->KeySwitchExt(ct, true);
    auto ct2Ext = cc->KeySwitchExt(ct2, true);

    // Get tower counts for reporting
    size_t towersQ = ct->GetElements()[0].GetNumOfElements();
    size_t towersPQ = ctExt->GetElements()[0].GetNumOfElements();
    std::cout << "    Q towers=" << towersQ << ", PQ towers=" << towersPQ << "\n";

    // mod Q: raw polynomial add (same method as PQ)
    {
        std::vector<double> times;
        for (int t = 0; t < trials; t++) {
            auto a = ct->Clone();
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < 100; i++) {
                auto& cv1 = a->GetElements();
                const auto& cv2 = ct2->GetElements();
                for (size_t j = 0; j < cv1.size(); j++) cv1[j] += cv2[j];
            }
            auto t1 = std::chrono::steady_clock::now();
            times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count() / 100.0);
        }
        auto st = meanStd(times);
        csv << ps.name << ",1c_arith,,polyAdd_modQ," << st.mean_ms << "," << st.std_ms << "\n";
        std::cout << "    polyAdd (mod Q): " << std::fixed << std::setprecision(4)
                  << st.mean_ms << " ms\n";
    }

    // mod Q: raw polynomial multiply
    {
        const auto& ptPoly = ct2->GetElements()[0];
        std::vector<double> times;
        for (int t = 0; t < trials; t++) {
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < 100; i++) {
                auto r = ct->Clone();
                auto& cv = r->GetElements();
                for (auto& c : cv) c *= ptPoly;
            }
            auto t1 = std::chrono::steady_clock::now();
            times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count() / 100.0);
        }
        auto st = meanStd(times);
        csv << ps.name << ",1c_arith,,polyMult_modQ," << st.mean_ms << "," << st.std_ms << "\n";
        std::cout << "    polyMult (mod Q): " << std::fixed << std::setprecision(4)
                  << st.mean_ms << " ms\n";
    }

    // mod PQ: raw polynomial add
    {
        std::vector<double> times;
        for (int t = 0; t < trials; t++) {
            auto a = ctExt->Clone();
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < 100; i++) {
                auto& cv1 = a->GetElements();
                const auto& cv2 = ct2Ext->GetElements();
                for (size_t j = 0; j < cv1.size(); j++) cv1[j] += cv2[j];
            }
            auto t1 = std::chrono::steady_clock::now();
            times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count() / 100.0);
        }
        auto st = meanStd(times);
        csv << ps.name << ",1c_arith,,polyAdd_modPQ," << st.mean_ms << "," << st.std_ms << "\n";
        std::cout << "    polyAdd (mod PQ): " << std::fixed << std::setprecision(4)
                  << st.mean_ms << " ms\n";
    }

    // mod PQ: raw polynomial multiply
    {
        const auto& ptPoly = ct2Ext->GetElements()[0];
        std::vector<double> times;
        for (int t = 0; t < trials; t++) {
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < 100; i++) {
                auto r = ctExt->Clone();
                auto& cv = r->GetElements();
                for (auto& c : cv) c *= ptPoly;
            }
            auto t1 = std::chrono::steady_clock::now();
            times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count() / 100.0);
        }
        auto st = meanStd(times);
        csv << ps.name << ",1c_arith,,polyMult_modPQ," << st.mean_ms << "," << st.std_ms << "\n";
        std::cout << "    ptMult (mod PQ): " << std::fixed << std::setprecision(4)
                  << st.mean_ms << " ms\n";
    }
}

// ============================================================
// Main
// ============================================================
int main(int argc, char** argv) {
    const int trials = (argc >= 2) ? std::max(1, std::atoi(argv[1])) : 10;

    std::cout << "=== Unit-Level Microbenchmarks ===\n";
    const char* omp = std::getenv("OMP_NUM_THREADS");
    std::cout << "OMP_NUM_THREADS = " << (omp ? omp : "default") << "\n";
    std::cout << "Trials = " << trials << "\n\n";

    auto presets = MakePresets();

    std::ofstream csv("unit_ops_bench.csv");
    csv << std::fixed << std::setprecision(6);
    csv << "preset,benchmark,k_or_idx,method,mean_ms,std_ms\n";

    for (const auto& ps : presets) {
        std::cout << "\n=== " << ps.name << " ===\n";
        bench_1a(ps, trials, csv);
        csv.flush();
        bench_1b(ps, trials, csv);
        csv.flush();
        bench_1c(ps, trials, csv);
        csv.flush();
    }

    csv.close();
    std::cout << "\nResults written to unit_ops_bench.csv\n";
    return 0;
}
