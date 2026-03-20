// Incremental tests for double hoisting building blocks.
// Tests each stage of key switching independently before combining.

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>

#include "openfhe.h"
#include "math/dftransform.h"
#include "encoding/ckkspackedencoding.h"

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
    P.SetKeySwitchTechnique(HYBRID);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    // Generate rotation keys for indices we'll test
    std::vector<int> indices = {0, 1, 2, 3, 5, 7, 8, 13};
    cc->EvalRotateKeyGen(kp.secretKey, indices);

    return {cc, kp.publicKey, kp.secretKey};
}

static double maxAbsErr(const std::vector<double>& a, const std::vector<double>& b, int n) {
    double mx = 0;
    for (int i = 0; i < n; i++) {
        double e = std::abs(a[i] - b[i]);
        if (e > mx) mx = e;
    }
    return mx;
}

static std::vector<double> decrypt(const ContextPack& ctx, ConstCiphertext<DCRTPoly> ct, int n) {
    Plaintext pt;
    ctx.cc->Decrypt(ctx.sk, ct, &pt);
    pt->SetLength(n);
    return pt->GetRealPackedValue();
}

// Constants needed for 128-bit encoding (only defined when NATIVEINT==128 in OpenFHE)
#if NATIVEINT != 128
static constexpr int128_t My_Max128BitValue() {
    return static_cast<int128_t>(((uint128_t)1 << 127) - ((uint128_t)1 << 73) - (uint128_t)1);
}
[[maybe_unused]] static bool my_is128BitOverflow(double d) {
    return std::abs(d) > static_cast<double>(My_Max128BitValue());
}
enum { MY_MAX_DOUBLE_PRECISION = 52 };
#endif

// ========== QP-basis plaintext encoding ==========
// Encodes real values directly into QP basis, matching MakeCKKSPackedPlaintext's
// NATIVEINT=64 Encode logic but outputting to all QP towers (not just Q).
// For FLEXIBLEAUTOEXT level 0: uses GetScalingFactorRealBig(0) as scale, NSD=1→2.
static Plaintext MakeQPPlaintext(const CryptoContext<DCRTPoly>& cc,
                                  const std::shared_ptr<DCRTPoly::Params>& paramsQP,
                                  const std::vector<double>& realValues,
                                  usint slots) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());

    // Match MakeCKKSPackedPlaintext: FLEXIBLEAUTOEXT level 0 uses BigScaling, NSD=1
    double scFact = cryptoParams->GetScalingFactorRealBig(0);
    size_t noiseScaleDeg = 1;  // forced for FLEXIBLEAUTOEXT level 0

    // Convert to complex
    std::vector<std::complex<double>> value(realValues.size());
    for (size_t i = 0; i < realValues.size(); i++)
        value[i] = {realValues[i], 0.0};
    value.resize(slots);
    for (size_t i = 0; i < value.size(); i++) value[i].imag(0.0);

    Plaintext p = Plaintext(std::make_shared<CKKSPackedEncoding>(
        paramsQP, cc->GetEncodingParams(), value, noiseScaleDeg, 0, scFact, slots));

    DCRTPoly& plainElement = p->GetElement<DCRTPoly>();
    usint N = cc->GetRingDimension();

    // IFFT (same as CKKSPackedEncoding::Encode)
    std::vector<std::complex<double>> inverse = value;
    inverse.resize(slots);
    DiscreteFourierTransform::FFTSpecialInv(inverse, N * 2);

    // Scale by scalingFactor (same as NATIVEINT=64 Encode path)
    double powP = scFact;

    constexpr int32_t MAX_BITS_IN_WORD = 61;
    int32_t logc = 0;
    for (size_t i = 0; i < slots; ++i) {
        inverse[i] *= powP;
        if (inverse[i].real() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].real()))));
            if (logc < logci) logc = logci;
        }
        if (inverse[i].imag() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].imag()))));
            if (logc < logci) logc = logci;
        }
    }

    int32_t logValid  = (logc <= MAX_BITS_IN_WORD) ? logc : MAX_BITS_IN_WORD;
    int32_t logApprox = logc - logValid;
    double approxFactor = pow(2, logApprox);

    // Round to int64 (same as NATIVEINT=64 Encode)
    std::vector<int64_t> temp(2 * slots);
    for (size_t i = 0; i < slots; ++i) {
        double dre = inverse[i].real() / approxFactor;
        double dim = inverse[i].imag() / approxFactor;
        int64_t re = std::llround(dre);
        int64_t im = std::llround(dim);
        temp[i]         = (re < 0) ? Max64BitValue() + re : re;
        temp[i + slots] = (im < 0) ? Max64BitValue() + im : im;
    }

    // FitToNativeVector for ALL QP towers (Q + P)
    const auto& nativeParams = plainElement.GetParams()->GetParams();
    int64_t bigBound = Max64BitValue();
    NativeInteger bigValueHf(static_cast<uint64_t>(bigBound >> 1));
    uint32_t dslots = temp.size();
    uint32_t gap = N / dslots;

    for (size_t i = 0; i < nativeParams.size(); i++) {
        NativeVector nativeVec(N, nativeParams[i]->GetModulus());
        NativeInteger modulus = nativeVec.GetModulus();
        NativeInteger diff = NativeInteger(static_cast<uint64_t>(bigBound)) - modulus;
        for (usint k = 0; k < dslots; k++) {
            NativeInteger n(static_cast<uint64_t>(temp[k]));
            if (n > bigValueHf)
                nativeVec[gap * k] = n.ModSub(diff, modulus);
            else
                nativeVec[gap * k] = n.Mod(modulus);
        }
        NativePoly element = plainElement.GetElementAtIndex(i);
        element.SetValues(std::move(nativeVec), Format::COEFFICIENT);
        plainElement.SetElementAtIndex(i, std::move(element));
    }

    // Handle approxFactor: multiply back by 2^logApprox in CRT (all QP towers)
    if (logApprox > 0) {
        usint numTowers = nativeParams.size();
        std::vector<DCRTPoly::Integer> moduli(numTowers);
        for (usint i = 0; i < numTowers; i++)
            moduli[i] = nativeParams[i]->GetModulus();

        constexpr int32_t MAX_LOG_STEP = 60;
        int32_t remLog = logApprox;
        int32_t logStep = (remLog <= MAX_LOG_STEP) ? remLog : MAX_LOG_STEP;
        DCRTPoly::Integer intStep = uint64_t(1) << logStep;
        std::vector<DCRTPoly::Integer> crtApprox(numTowers, intStep);
        remLog -= logStep;
        while (remLog > 0) {
            logStep = (remLog <= MAX_LOG_STEP) ? remLog : MAX_LOG_STEP;
            intStep = uint64_t(1) << logStep;
            std::vector<DCRTPoly::Integer> crtSF(numTowers, intStep);
            crtApprox = CKKSPackedEncoding::CRTMult(crtApprox, crtSF, moduli);
            remLog -= logStep;
        }
        plainElement = plainElement.Times(crtApprox);
    }

    // noiseScaleDeg=1 for FLEXIBLEAUTOEXT level 0, no CRT powP multiplication needed

    p->SetFormat(Format::EVALUATION);
    // Match MakeCKKSPackedPlaintext post-encode: NSD set to 2, SF stays as-is
    p->SetNoiseScaleDeg(2);
    // scalingFactor = pow(scFact, 1) = scFact (already set by constructor)

    return p;
}

// ========== Test 1: KeySwitchExt → KeySwitchDown == identity ==========
// Extend ct to QP basis, then bring back to Q. Should recover original.
static void test1_ext_down(const ContextPack& ctx) {
    std::cout << "=== Test 1: KeySwitchExt(ct, true) → KeySwitchDown → compare with ct ===\n";

    const int n = 64;
    std::vector<double> vals(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);

    auto pt = ctx.cc->MakeCKKSPackedPlaintext(vals);
    auto ct = ctx.cc->Encrypt(ctx.pk, pt);

    // Extend to QP, then ModDown back to Q
    auto ctExt = ctx.cc->KeySwitchExt(ct, true);
    auto ctBack = ctx.cc->KeySwitchDown(ctExt);

    auto dec_orig = decrypt(ctx, ct, n);
    auto dec_back = decrypt(ctx, ctBack, n);

    double err = maxAbsErr(dec_orig, dec_back, n);
    std::cout << "  max|err| = " << std::scientific << std::setprecision(6) << err << "\n";
    std::cout << "  " << (err < 1e-3 ? "PASS" : "FAIL") << "\n\n";
}

// ========== Test 2: EvalFastRotationExt + KeySwitchDown == EvalFastRotation ==========
// Compare the two rotation paths for several rotation indices.
static void test2_rotation_paths(const ContextPack& ctx) {
    std::cout << "=== Test 2: EvalFastRotationExt + KeySwitchDown vs EvalFastRotation ===\n";

    const int n = 64;
    std::vector<double> vals(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);

    auto pt = ctx.cc->MakeCKKSPackedPlaintext(vals);
    auto ct = ctx.cc->Encrypt(ctx.pk, pt);

    auto N = ctx.cc->GetRingDimension();
    auto M = 2 * N;
    auto digits = ctx.cc->EvalFastRotationPrecompute(ct);

    std::vector<int> testIndices = {1, 2, 3, 5, 7, 8, 13};
    for (int r : testIndices) {
        // Path A: standard EvalFastRotation (Q basis)
        auto ctA = ctx.cc->EvalFastRotation(ct, r, M, digits);

        // Path B: EvalFastRotationExt (QP basis) → KeySwitchDown
        auto ctExtB = ctx.cc->EvalFastRotationExt(ct, r, digits, true);
        auto ctB = ctx.cc->KeySwitchDown(ctExtB);

        auto decA = decrypt(ctx, ctA, n);
        auto decB = decrypt(ctx, ctB, n);

        double err = maxAbsErr(decA, decB, n);
        std::cout << "  rot=" << std::setw(2) << r
                  << " : max|err| = " << std::scientific << std::setprecision(6) << err
                  << " " << (err < 1e-3 ? "PASS" : "FAIL") << "\n";
    }
    std::cout << "\n";
}

// ========== Test 2.5: QP-basis addition ==========
// (RotExt(ct, 1) + RotExt(ct, 3)) → KeySwitchDown
// vs. Rot(ct, 1) + Rot(ct, 3)
static void test2_5_ext_add(const ContextPack& ctx) {
    std::cout << "=== Test 2.5: QP addition then KeySwitchDown vs Q addition ===\n";

    const int n = 64;
    std::vector<double> vals(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);

    auto pt = ctx.cc->MakeCKKSPackedPlaintext(vals);
    auto ct = ctx.cc->Encrypt(ctx.pk, pt);

    auto N = ctx.cc->GetRingDimension();
    auto M = 2 * N;
    auto digits = ctx.cc->EvalFastRotationPrecompute(ct);

    // Path A: Q-basis rotations + add
    auto rotA1 = ctx.cc->EvalFastRotation(ct, 1, M, digits);
    auto rotA3 = ctx.cc->EvalFastRotation(ct, 3, M, digits);
    ctx.cc->EvalAddInPlace(rotA1, rotA3);

    // Path B: QP-basis rotations + add in QP + KeySwitchDown
    auto rotB1 = ctx.cc->EvalFastRotationExt(ct, 1, digits, true);
    auto rotB3 = ctx.cc->EvalFastRotationExt(ct, 3, digits, true);
    // Add in QP
    auto& cv1 = rotB1->GetElements();
    const auto& cv3 = rotB3->GetElements();
    for (size_t i = 0; i < cv1.size(); i++) cv1[i] += cv3[i];
    auto sumB = ctx.cc->KeySwitchDown(rotB1);

    auto decA = decrypt(ctx, rotA1, n);
    auto decB = decrypt(ctx, sumB, n);

    double err = maxAbsErr(decA, decB, n);
    std::cout << "  Rot(1)+Rot(3) : max|err| = " << std::scientific << std::setprecision(6) << err
              << " " << (err < 1e-3 ? "PASS" : "FAIL") << "\n\n";
}

// ========== Helper: extend Q-basis plaintext to QP via ApproxModUp ==========
[[maybe_unused]] static DCRTPoly ExtendPtToQP(const Plaintext& ptQ,
                              const std::shared_ptr<DCRTPoly::Params>& paramsQP) {
    DCRTPoly ptPoly = ptQ->GetElement<DCRTPoly>();
    size_t sizeQ = ptPoly.GetNumOfElements();
    size_t sizeQP = paramsQP->GetParams().size();
    if (sizeQ >= sizeQP) {
        ptPoly.SetFormat(Format::EVALUATION);
        return ptPoly;
    }
    size_t sizeP = sizeQP - sizeQ;

    std::vector<NativeInteger> moduliQ(sizeQ), rootsQ(sizeQ);
    std::vector<NativeInteger> moduliP(sizeP), rootsP(sizeP);
    for (size_t i = 0; i < sizeQ; i++) {
        moduliQ[i] = paramsQP->GetParams()[i]->GetModulus();
        rootsQ[i]  = paramsQP->GetParams()[i]->GetRootOfUnity();
    }
    for (size_t j = 0; j < sizeP; j++) {
        moduliP[j] = paramsQP->GetParams()[sizeQ + j]->GetModulus();
        rootsP[j]  = paramsQP->GetParams()[sizeQ + j]->GetRootOfUnity();
    }
    uint32_t ringDim = paramsQP->GetRingDimension();
    auto pQ = std::make_shared<DCRTPoly::Params>(2 * ringDim, moduliQ, rootsQ);
    auto pP = std::make_shared<DCRTPoly::Params>(2 * ringDim, moduliP, rootsP);

    std::vector<NativeInteger> QHatInvModq(sizeQ), QHatInvModqPrecon(sizeQ);
    for (size_t i = 0; i < sizeQ; i++) {
        uint64_t qi = moduliQ[i].ConvertToInt<uint64_t>();
        __uint128_t qhat_mod_qi = 1;
        for (size_t k = 0; k < sizeQ; k++) {
            if (k == i) continue;
            qhat_mod_qi = (qhat_mod_qi * (moduliQ[k].ConvertToInt<uint64_t>() % qi)) % qi;
        }
        QHatInvModq[i] = NativeInteger(static_cast<uint64_t>(qhat_mod_qi)).ModInverse(moduliQ[i]);
        QHatInvModqPrecon[i] = QHatInvModq[i].PrepModMulConst(moduliQ[i]);
    }

    std::vector<std::vector<NativeInteger>> QHatModp(sizeQ);
    for (size_t i = 0; i < sizeQ; i++) {
        QHatModp[i].resize(sizeP);
        for (size_t j = 0; j < sizeP; j++) {
            uint64_t pj = moduliP[j].ConvertToInt<uint64_t>();
            __uint128_t prod = 1;
            for (size_t k = 0; k < sizeQ; k++) {
                if (k == i) continue;
                prod = (prod * (moduliQ[k].ConvertToInt<uint64_t>() % pj)) % pj;
            }
            QHatModp[i][j] = NativeInteger(static_cast<uint64_t>(prod));
        }
    }

    const auto BarrettBase128Bit(BigInteger(1).LShiftEq(128));
    std::vector<DoubleNativeInt> modpBarrettMu(sizeP);
    for (size_t j = 0; j < sizeP; j++) {
        modpBarrettMu[j] =
            (BarrettBase128Bit / BigInteger(moduliP[j])).ConvertToInt<DoubleNativeInt>();
    }

    ptPoly.ApproxModUp(pQ, pP, paramsQP, QHatInvModq, QHatInvModqPrecon, QHatModp, modpBarrettMu);
    return ptPoly;
}

// ========== Test 3: QP-basis plaintext multiplication ==========
static void test3_ext_mult_down(const ContextPack& ctx) {
    std::cout << "=== Test 3: QP mult then KeySwitchDown vs Q mult ===\n";

    const int n = 64;
    std::vector<double> vals(n), mask(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);
    for (int i = 0; i < n; i++) mask[i] = (i % 2 == 0) ? 1.0 : 0.0;

    auto ct = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(vals));
    auto ptMask = ctx.cc->MakeCKKSPackedPlaintext(mask);

    auto N = ctx.cc->GetRingDimension();
    auto M = 2 * N;
    auto digits = ctx.cc->EvalFastRotationPrecompute(ct);
    int r = 3;

    auto ctExtB = ctx.cc->EvalFastRotationExt(ct, r, digits, true);
    auto paramsQP = ctExtB->GetElements()[0].GetParams();

    // Create QP plaintext via direct encoding (matching MakeCKKSPackedPlaintext's scaling)
    auto ptQP = MakeQPPlaintext(ctx.cc, paramsQP, mask, n);
    DCRTPoly ptPolyQP = ptQP->GetElement<DCRTPoly>();
    ptPolyQP.SetFormat(Format::EVALUATION);
    std::cout << "  ptPolyQP towers=" << ptPolyQP.GetNumOfElements()
              << " NSD=" << ptQP->GetNoiseScaleDeg() << " SF=" << ptQP->GetScalingFactor() << "\n";

    // Path A: Q-basis rotation + manual Q mult (no rescale)
    auto ctA = ctx.cc->EvalFastRotation(ct, r, M, digits);
    {
        DCRTPoly ptPolyQ = ptMask->GetElement<DCRTPoly>();
        ptPolyQ.SetFormat(Format::EVALUATION);
        auto& cvA = ctA->GetElements();
        for (auto& c : cvA) c *= ptPolyQ;
        ctA->SetNoiseScaleDeg(ctA->GetNoiseScaleDeg() + ptMask->GetNoiseScaleDeg());
        ctA->SetScalingFactor(ctA->GetScalingFactor() * ptMask->GetScalingFactor());
    }

    // Path B: QP rotation * QP plaintext → KeySwitchDown
    auto resultB = ctExtB->Clone();
    {
        auto& cv = resultB->GetElements();
        for (auto& c : cv) c *= ptPolyQP;
        resultB->SetNoiseScaleDeg(resultB->GetNoiseScaleDeg() + ptQP->GetNoiseScaleDeg());
        resultB->SetScalingFactor(resultB->GetScalingFactor() * ptQP->GetScalingFactor());
    }
    auto ctB = ctx.cc->KeySwitchDown(resultB);

    std::cout << "  [PathA] NSD=" << ctA->GetNoiseScaleDeg() << " SF=" << ctA->GetScalingFactor()
              << " towers=" << ctA->GetElements()[0].GetNumOfElements() << "\n";
    std::cout << "  [PathB] NSD=" << ctB->GetNoiseScaleDeg() << " SF=" << ctB->GetScalingFactor()
              << " towers=" << ctB->GetElements()[0].GetNumOfElements() << "\n";

    try {
        auto decA = decrypt(ctx, ctA, n);
        auto decB = decrypt(ctx, ctB, n);
        double err = maxAbsErr(decA, decB, n);
        std::cout << "  rot=" << r << " : max|err| = " << std::scientific << std::setprecision(6) << err
                  << " " << (err < 1e-3 ? "PASS" : "FAIL") << "\n";
    } catch (const std::exception& e) {
        std::cout << "  Decrypt FAILED: " << e.what() << "\n";
    }
    std::cout << "\n";
}

// ========== Test 4: QP-basis accumulation (sum of products) ==========
// sum_i(RotExt(ct, r_i) * pt_i) → KeySwitchDown  vs  sum_i(Rot(ct, r_i) * pt_i)
static void test4_ext_accumulate(const ContextPack& ctx) {
    std::cout << "=== Test 4: QP accumulate then KeySwitchDown vs Q accumulate ===\n";

    const int n = 64;
    std::vector<double> vals(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);

    auto ct = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(vals));

    auto N = ctx.cc->GetRingDimension();
    auto M = 2 * N;
    auto digits = ctx.cc->EvalFastRotationPrecompute(ct);

    std::vector<int> rotations = {1, 2, 3};
    std::vector<std::vector<double>> masks(rotations.size());
    for (size_t idx = 0; idx < rotations.size(); idx++) {
        masks[idx].resize(n);
        for (int i = 0; i < n; i++)
            masks[idx][i] = (i % (idx + 2) == 0) ? 1.0 : 0.0;
    }

    // Get QP params
    auto ctExt0 = ctx.cc->EvalFastRotationExt(ct, rotations[0], digits, true);
    auto paramsQP = ctExt0->GetElements()[0].GetParams();

    // Path A: Q-basis
    Ciphertext<DCRTPoly> sumA;
    for (size_t idx = 0; idx < rotations.size(); idx++) {
        auto rot = ctx.cc->EvalFastRotation(ct, rotations[idx], M, digits);
        auto ptM = ctx.cc->MakeCKKSPackedPlaintext(masks[idx]);
        DCRTPoly ptPolyQ = ptM->GetElement<DCRTPoly>();
        ptPolyQ.SetFormat(Format::EVALUATION);
        auto& cvR = rot->GetElements();
        for (auto& c : cvR) c *= ptPolyQ;
        rot->SetNoiseScaleDeg(rot->GetNoiseScaleDeg() + ptM->GetNoiseScaleDeg());
        rot->SetScalingFactor(rot->GetScalingFactor() * ptM->GetScalingFactor());
        if (idx == 0) sumA = rot;
        else {
            auto& cv1 = sumA->GetElements();
            const auto& cv2 = rot->GetElements();
            for (size_t i = 0; i < cv1.size(); i++) cv1[i] += cv2[i];
        }
    }

    // Path B: QP-basis with direct QP encoding
    Ciphertext<DCRTPoly> sumExtB;
    for (size_t idx = 0; idx < rotations.size(); idx++) {
        auto ctExt = ctx.cc->EvalFastRotationExt(ct, rotations[idx], digits, true);
        auto ptQP = MakeQPPlaintext(ctx.cc, paramsQP, masks[idx], n);
        DCRTPoly ptPolyQP = ptQP->GetElement<DCRTPoly>();
        ptPolyQP.SetFormat(Format::EVALUATION);

        auto termB = ctExt->Clone();
        auto& cv = termB->GetElements();
        for (auto& c : cv) c *= ptPolyQP;
        termB->SetNoiseScaleDeg(termB->GetNoiseScaleDeg() + ptQP->GetNoiseScaleDeg());
        termB->SetScalingFactor(termB->GetScalingFactor() * ptQP->GetScalingFactor());

        if (idx == 0) {
            sumExtB = termB;
        } else {
            auto& cv1 = sumExtB->GetElements();
            const auto& cv2 = termB->GetElements();
            for (size_t i = 0; i < cv1.size(); i++) cv1[i] += cv2[i];
        }
    }
    auto sumB = ctx.cc->KeySwitchDown(sumExtB);

    try {
        auto decA = decrypt(ctx, sumA, n);
        auto decB = decrypt(ctx, sumB, n);
        double err = maxAbsErr(decA, decB, n);
        std::cout << "  accumulated sum : max|err| = "
                  << std::scientific << std::setprecision(6) << err
                  << " " << (err < 1e-3 ? "PASS" : "FAIL") << "\n";
    } catch (const std::exception& e) {
        std::cout << "  FAILED: " << e.what() << "\n";
    }
    std::cout << "\n";
}

// ========== QP-basis plaintext encoding at NSD=1, matching MakeAuxPlaintext ==========
// Encodes with GetScalingFactorReal(level) and NSD=1.
// This is what the bootstrapping code uses for its double hoisting.
static Plaintext MakeQPPlaintextNSD1(const CryptoContext<DCRTPoly>& cc,
                                      const std::shared_ptr<DCRTPoly::Params>& paramsQP,
                                      const std::vector<double>& realValues,
                                      usint slots, uint32_t level) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    double scFact = cryptoParams->GetScalingFactorReal(level);
    size_t noiseScaleDeg = 1;

    std::vector<std::complex<double>> value(realValues.size());
    for (size_t i = 0; i < realValues.size(); i++)
        value[i] = {realValues[i], 0.0};
    value.resize(slots);

    Plaintext p = Plaintext(std::make_shared<CKKSPackedEncoding>(
        paramsQP, cc->GetEncodingParams(), value, noiseScaleDeg, level, scFact, slots));

    DCRTPoly& plainElement = p->GetElement<DCRTPoly>();
    usint N = cc->GetRingDimension();

    std::vector<std::complex<double>> inverse = value;
    inverse.resize(slots);
    DiscreteFourierTransform::FFTSpecialInv(inverse, N * 2);

    double powP = scFact;
    constexpr int32_t MAX_BITS_IN_WORD = 61;
    int32_t logc = 0;
    for (size_t i = 0; i < slots; ++i) {
        inverse[i] *= powP;
        if (inverse[i].real() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].real()))));
            if (logc < logci) logc = logci;
        }
        if (inverse[i].imag() != 0) {
            int32_t logci = static_cast<int32_t>(ceil(log2(std::abs(inverse[i].imag()))));
            if (logc < logci) logc = logci;
        }
    }

    int32_t logValid  = (logc <= MAX_BITS_IN_WORD) ? logc : MAX_BITS_IN_WORD;
    int32_t logApprox = logc - logValid;
    double approxFactor = pow(2, logApprox);

    std::vector<int64_t> temp(2 * slots);
    for (size_t i = 0; i < slots; ++i) {
        double dre = inverse[i].real() / approxFactor;
        double dim = inverse[i].imag() / approxFactor;
        int64_t re = std::llround(dre);
        int64_t im = std::llround(dim);
        temp[i]         = (re < 0) ? Max64BitValue() + re : re;
        temp[i + slots] = (im < 0) ? Max64BitValue() + im : im;
    }

    const auto& nativeParams = plainElement.GetParams()->GetParams();
    int64_t bigBound = Max64BitValue();
    NativeInteger bigValueHf(static_cast<uint64_t>(bigBound >> 1));
    uint32_t dslots = temp.size();
    uint32_t gap = N / dslots;

    for (size_t i = 0; i < nativeParams.size(); i++) {
        NativeVector nativeVec(N, nativeParams[i]->GetModulus());
        NativeInteger modulus = nativeVec.GetModulus();
        NativeInteger diff = NativeInteger(static_cast<uint64_t>(bigBound)) - modulus;
        for (usint k = 0; k < dslots; k++) {
            NativeInteger n(static_cast<uint64_t>(temp[k]));
            if (n > bigValueHf)
                nativeVec[gap * k] = n.ModSub(diff, modulus);
            else
                nativeVec[gap * k] = n.Mod(modulus);
        }
        NativePoly element = plainElement.GetElementAtIndex(i);
        element.SetValues(std::move(nativeVec), Format::COEFFICIENT);
        plainElement.SetElementAtIndex(i, std::move(element));
    }

    if (logApprox > 0) {
        usint numTowers = nativeParams.size();
        std::vector<DCRTPoly::Integer> moduli(numTowers);
        for (usint i = 0; i < numTowers; i++)
            moduli[i] = nativeParams[i]->GetModulus();
        constexpr int32_t MAX_LOG_STEP = 60;
        int32_t remLog = logApprox;
        int32_t logStep = (remLog <= MAX_LOG_STEP) ? remLog : MAX_LOG_STEP;
        DCRTPoly::Integer intStep = uint64_t(1) << logStep;
        std::vector<DCRTPoly::Integer> crtApprox(numTowers, intStep);
        remLog -= logStep;
        while (remLog > 0) {
            logStep = (remLog <= MAX_LOG_STEP) ? remLog : MAX_LOG_STEP;
            intStep = uint64_t(1) << logStep;
            std::vector<DCRTPoly::Integer> crtSF(numTowers, intStep);
            crtApprox = CKKSPackedEncoding::CRTMult(crtApprox, crtSF, moduli);
            remLog -= logStep;
        }
        plainElement = plainElement.Times(crtApprox);
    }

    p->SetFormat(Format::EVALUATION);
    p->SetNoiseScaleDeg(1);
    return p;
}

// ========== Test 5: Full BSGS pattern — multiple giant steps with QP accumulation ==========
// Key insight: pre-rescale ct to NSD=1 before extending to QP, encode plaintexts at NSD=1.
// This matches how OpenFHE's bootstrapping code does double hoisting.
static void test5_bsgs_pattern(const ContextPack& ctx) {
    std::cout << "=== Test 5: Full BSGS pattern (pre-rescale → QP accumulate → KSDown → Rotate → Add) ===\n";

    const int n = 64;
    std::vector<double> vals(n);
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    for (auto& v : vals) v = dis(gen);

    auto ct = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(vals));
    std::cout << "  ct: NSD=" << ct->GetNoiseScaleDeg() << " SF=" << ct->GetScalingFactor()
              << " lvl=" << ct->GetLevel() << " towers=" << ct->GetElements()[0].GetNumOfElements() << "\n";

    auto N = ctx.cc->GetRingDimension();
    auto M = 2 * N;

    // BSGS parameters: 3 baby steps, 2 giant steps
    const int bs = 3;
    const int gs = 2;

    // Generate random masks
    std::vector<std::vector<std::vector<double>>> masks(gs,
        std::vector<std::vector<double>>(bs, std::vector<double>(n)));
    for (int i = 0; i < gs; i++)
        for (int j = 0; j < bs; j++)
            for (int k = 0; k < n; k++)
                masks[i][j][k] = dis(gen);

    // ------ Path A: Q-basis (mimic hoisting version using cc->EvalMult) ------
    auto digitsA = ctx.cc->EvalFastRotationPrecompute(ct);
    std::vector<Ciphertext<DCRTPoly>> babyStepsQ(bs);
    for (int j = 0; j < bs; j++)
        babyStepsQ[j] = ctx.cc->EvalFastRotation(ct, j, M, digitsA);

    std::vector<double> zero(n, 0.0);
    auto sumA = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(zero));

    for (int i = 0; i < gs; i++) {
        auto tmp = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(zero));
        for (int j = 0; j < bs; j++) {
            auto ptM = ctx.cc->MakeCKKSPackedPlaintext(masks[i][j]);
            ctx.cc->EvalAddInPlace(tmp, ctx.cc->EvalMult(ptM, babyStepsQ[j]));
        }
        ctx.cc->EvalAddInPlace(sumA, ctx.cc->EvalRotate(tmp, bs * i));
    }
    std::cout << "  [sumA] NSD=" << sumA->GetNoiseScaleDeg() << " SF=" << sumA->GetScalingFactor()
              << " lvl=" << sumA->GetLevel() << "\n";

    // ------ Path B: QP-basis (double hoisting) ------
    // QP plaintext at NSD=1 → after multiply NSD=3 → ModReduce → NSD=2, level 1
    // This matches EvalMult output without unfair overhead in EvalAdd.
    auto digitsB = ctx.cc->EvalFastRotationPrecompute(ct);
    std::vector<Ciphertext<DCRTPoly>> babyStepsExt(bs);
    babyStepsExt[0] = ctx.cc->KeySwitchExt(ct, true);
    for (int j = 1; j < bs; j++)
        babyStepsExt[j] = ctx.cc->EvalFastRotationExt(ct, j, digitsB, true);
    auto paramsQP = babyStepsExt[0]->GetElements()[0].GetParams();

    uint32_t level = ct->GetLevel();
    auto sumB = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(zero));

    for (int i = 0; i < gs; i++) {
        Ciphertext<DCRTPoly> accExt;
        bool first = true;
        for (int j = 0; j < bs; j++) {
            // Encode plaintext at NSD=1 with GetScalingFactorReal(level)
            auto ptQP = MakeQPPlaintextNSD1(ctx.cc, paramsQP, masks[i][j], n, level);
            DCRTPoly ptPolyQP = ptQP->GetElement<DCRTPoly>();
            ptPolyQP.SetFormat(Format::EVALUATION);

            auto term = babyStepsExt[j]->Clone();
            auto& cv = term->GetElements();
            for (auto& c : cv) c *= ptPolyQP;
            // Update NSD/SF like EvalMultExt does
            term->SetNoiseScaleDeg(term->GetNoiseScaleDeg() + ptQP->GetNoiseScaleDeg());
            term->SetScalingFactor(term->GetScalingFactor() * ptQP->GetScalingFactor());

            if (first) { accExt = term; first = false; }
            else {
                auto& cv1 = accExt->GetElements();
                const auto& cv2 = term->GetElements();
                for (size_t k = 0; k < cv1.size(); k++) cv1[k] += cv2[k];
            }
        }
        auto tmp = ctx.cc->KeySwitchDown(accExt);
        std::cout << "  [B gs=" << i << "] KSDown: NSD=" << tmp->GetNoiseScaleDeg()
                  << " SF=" << tmp->GetScalingFactor() << " lvl=" << tmp->GetLevel()
                  << " towers=" << tmp->GetElements()[0].GetNumOfElements() << "\n";
        // ModReduce: NSD=3→2, level 0→1 (matches EvalMult output, no unfair overhead)
        ctx.cc->ModReduceInPlace(tmp);
        std::cout << "  [B gs=" << i << "] after ModReduce: NSD=" << tmp->GetNoiseScaleDeg()
                  << " SF=" << tmp->GetScalingFactor() << " lvl=" << tmp->GetLevel()
                  << " towers=" << tmp->GetElements()[0].GetNumOfElements() << "\n";
        ctx.cc->EvalAddInPlace(sumB, ctx.cc->EvalRotate(tmp, bs * i));
    }
    std::cout << "  [sumB] NSD=" << sumB->GetNoiseScaleDeg() << " SF=" << sumB->GetScalingFactor()
              << " lvl=" << sumB->GetLevel() << "\n";

    // ------ Compare ------
    try {
        auto decA = decrypt(ctx, sumA, n);
        std::cout << "  PathA decrypt OK. [0:4] = ";
        for (int i = 0; i < 4; i++) std::cout << std::fixed << std::setprecision(6) << decA[i] << " ";
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "  PathA FAILED: " << e.what() << "\n";
    }
    try {
        auto decB = decrypt(ctx, sumB, n);
        std::cout << "  PathB decrypt OK. [0:4] = ";
        for (int i = 0; i < 4; i++) std::cout << std::fixed << std::setprecision(6) << decB[i] << " ";
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "  PathB FAILED: " << e.what() << "\n";
    }
    // Compare if both succeeded
    try {
        auto decA = decrypt(ctx, sumA, n);
        auto decB = decrypt(ctx, sumB, n);
        double err = maxAbsErr(decA, decB, n);
        std::cout << "  max|err| = " << std::scientific << std::setprecision(6) << err
                  << " " << (err < 1e-2 ? "PASS" : "FAIL") << "\n";
    } catch (...) {}
    std::cout << "\n";
}

// ========== Main ==========
int main() {
    std::cout << "Building context (N=2^14, L=7, HYBRID)...\n\n";
    auto ctx = makeContext();

    // Quick debug: check scaling factor
    {
        const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ctx.cc->GetCryptoParameters());
        std::cout << "GetScalingFactorReal(0) = " << cryptoParams->GetScalingFactorReal(0) << "\n";
        std::cout << "GetScalingFactorReal(1) = " << cryptoParams->GetScalingFactorReal(1) << "\n";
        auto ptQ = ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>{1.0});
        std::cout << "MakeCKKSPackedPlaintext: NSD=" << ptQ->GetNoiseScaleDeg()
                  << " SF=" << ptQ->GetScalingFactor() << "\n\n";
    }

    test1_ext_down(ctx);
    test2_rotation_paths(ctx);
    test2_5_ext_add(ctx);
    test3_ext_mult_down(ctx);
    test4_ext_accumulate(ctx);
    test5_bsgs_pattern(ctx);

    std::cout << "All tests complete.\n";
    return 0;
}
