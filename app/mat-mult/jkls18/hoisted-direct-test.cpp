// Correctness and performance test for hoisted EvalDirectRotate.
// Compares: EvalDirectRotate(ct, i) vs EvalDirectRotate(ct, i, precomp)
// and measures timing for k rotations with/without hoisting.

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <numeric>

#include "openfhe.h"

using namespace lbcrypto;

// ========== Helpers ==========

struct ContextPack {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
};

static ContextPack makeContext() {
    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(7);
    P.SetRingDim(1 << 14);
    P.SetScalingModSize(34);
    P.SetFirstModSize(46);
    P.SetNumLargeDigits(3);
    P.SetSecurityLevel(HEStd_128_classic);
    P.SetBatchSize(64);
    P.SetKeySwitchTechnique(BATCHED);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    return {cc, kp.publicKey, kp.secretKey};
}

static std::vector<double> decrypt(const ContextPack& ctx, ConstCiphertext<DCRTPoly> ct, int n) {
    Plaintext pt;
    ctx.cc->Decrypt(ctx.sk, ct, &pt);
    pt->SetLength(n);
    return pt->GetRealPackedValue();
}

static double maxAbsErr(const std::vector<double>& a, const std::vector<double>& b, int n) {
    double mx = 0;
    for (int i = 0; i < n; i++) {
        double e = std::abs(a[i] - b[i]);
        if (e > mx) mx = e;
    }
    return mx;
}

// Rotate a plaintext vector left by `index` positions (with wraparound)
static std::vector<double> rotatePlain(const std::vector<double>& v, int32_t index) {
    int n = static_cast<int>(v.size());
    std::vector<double> out(n);
    for (int i = 0; i < n; i++) {
        out[i] = v[((i + index) % n + n) % n];
    }
    return out;
}

// ========== Main ==========

int main() {
    std::cout << "===== Hoisted EvalDirectRotate Test =====" << std::endl;

    auto ctx = makeContext();
    const int N = 64;

    // Rotation indices to test
    std::vector<int32_t> indices = {1, 2, 3, 5, 7, 8, 13, -1, -3};

    // Generate lazy rotation keys
    ctx.cc->EvalLazyAtIndexKeyGen(ctx.sk, indices);

    // Create test vector
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> vals(N);
    for (auto& v : vals) v = dist(rng);

    auto pt = ctx.cc->MakeCKKSPackedPlaintext(vals);
    auto ct = ctx.cc->Encrypt(ctx.pk, pt);

    // ====== Correctness Test ======
    std::cout << "\n--- Correctness: hoisted vs non-hoisted EvalDirectRotate ---" << std::endl;
    std::cout << std::setw(10) << "Index"
              << std::setw(22) << "Err(hoist vs plain)"
              << std::setw(22) << "Err(norm vs hoist)" << std::endl;

    auto precomp = ctx.cc->EvalDirectRotatePrecompute(ct);

    bool allPass = true;
    for (int32_t idx : indices) {
        auto res_normal  = ctx.cc->EvalDirectRotate(ct, idx);
        auto res_hoisted = ctx.cc->EvalDirectRotate(ct, idx, precomp);

        auto dec_normal  = decrypt(ctx, res_normal, N);
        auto dec_hoisted = decrypt(ctx, res_hoisted, N);
        auto expected    = rotatePlain(vals, idx);

        double err_vs_plain  = maxAbsErr(dec_hoisted, expected, N);
        double err_vs_normal = maxAbsErr(dec_normal, dec_hoisted, N);

        bool pass = (err_vs_plain < 1e-3) && (err_vs_normal < 1e-6);
        std::cout << std::setw(10) << idx
                  << std::setw(22) << std::scientific << std::setprecision(6) << err_vs_plain
                  << std::setw(22) << err_vs_normal
                  << (pass ? "  PASS" : "  FAIL") << std::endl;

        if (!pass) allPass = false;
    }

    std::cout << "\nOverall: " << (allPass ? "ALL PASSED" : "SOME FAILED") << std::endl;

    // ====== Performance Test ======
    std::cout << "\n--- Performance Test ---" << std::endl;
    const int kTrials = 5;
    const int kRotations = static_cast<int>(indices.size());

    // Non-hoisted: k separate EvalDirectRotate calls
    {
        std::vector<double> times;
        for (int t = 0; t < kTrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            for (int32_t idx : indices) {
                auto res = ctx.cc->EvalDirectRotate(ct, idx);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start).count();
            times.push_back(ms);
        }
        double mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        std::cout << "Non-hoisted (" << kRotations << " rotations): "
                  << std::fixed << std::setprecision(2) << mean << " ms" << std::endl;
    }

    // Hoisted: 1 precompute + k hoisted rotations
    {
        std::vector<double> times;
        for (int t = 0; t < kTrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            auto digits = ctx.cc->EvalDirectRotatePrecompute(ct);
            for (int32_t idx : indices) {
                auto res = ctx.cc->EvalDirectRotate(ct, idx, digits);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start).count();
            times.push_back(ms);
        }
        double mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        std::cout << "Hoisted    (" << kRotations << " rotations): "
                  << std::fixed << std::setprecision(2) << mean << " ms" << std::endl;
    }

    return allPass ? 0 : 1;
}
