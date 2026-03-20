#include <algorithm>
#include <cmath>
#include "jkls18.h"
#include "math/dftransform.h"
#include "math/nbtheory.h"
#include "encoding/ckkspackedencoding.h"

using namespace lbcrypto;

// ---- Double-hoisting helpers: QP-basis arithmetic ----

// Encode real values into QP-basis plaintext matching MakeCKKSPackedPlaintext's scale/NSD.
static Plaintext MakeQPPlaintext(const CryptoContext<DCRTPoly>& cc,
                                  const std::shared_ptr<DCRTPoly::Params>& paramsQP,
                                  const std::vector<double>& realValues,
                                  uint32_t level = 0) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    double scFact;
    size_t noiseScaleDeg = 1;
    if (level == 0 && cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT) {
        scFact = cryptoParams->GetScalingFactorRealBig(0);
    } else {
        scFact = cryptoParams->GetScalingFactorReal(level);
    }
    usint slots = cc->GetEncodingParams()->GetBatchSize();
    if (slots == 0) slots = cc->GetRingDimension() / 2;

    std::vector<std::complex<double>> value(realValues.size());
    for (size_t i = 0; i < realValues.size(); i++)
        value[i] = {realValues[i], 0.0};
    value.resize(slots);
    for (size_t i = 0; i < value.size(); i++) value[i].imag(0.0);

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
    p->SetNoiseScaleDeg(noiseScaleDeg);
    return p;
}

// Multiply QP-basis ciphertext by QP-basis plaintext element-wise.
// Updates NSD and SF to track the additional plaintext factor (like EvalMultExt).
static Ciphertext<DCRTPoly> MultExtQP(
    ConstCiphertext<DCRTPoly> ctQP,
    const DCRTPoly& ptPolyQP,
    const Plaintext& ptQP)
{
    Ciphertext<DCRTPoly> result = ctQP->Clone();
    auto& cv = result->GetElements();
    for (auto& c : cv) c *= ptPolyQP;
    result->SetNoiseScaleDeg(result->GetNoiseScaleDeg() + ptQP->GetNoiseScaleDeg());
    result->SetScalingFactor(result->GetScalingFactor() * ptQP->GetScalingFactor());
    return result;
}

// Add two QP-basis ciphertexts in place.
static void AddExtInPlace(
    Ciphertext<DCRTPoly>& ct1,
    ConstCiphertext<DCRTPoly> ct2)
{
    auto& cv1 = ct1->GetElements();
    const auto& cv2 = ct2->GetElements();
    for (size_t i = 0; i < cv1.size(); i++) {
        cv1[i] += cv2[i];
    }
}

MATMULT_JKLS18::MATMULT_JKLS18(const CryptoContext<DCRTPoly>& cc,
                           const PublicKey<DCRTPoly>& pk,
                           int dim)
    : m_cc(cc), m_PublicKey(pk), d(dim) {}

MATMULT_JKLS18::MATMULT_JKLS18(int dim)
    : d(dim) {}

std::vector<double> MATMULT_JKLS18::generateSigmaMsk(int k) {
    std::vector<double> msk(d * d, 0);
        if (k >= 0) {
            for (int i = d * k; i < d - k + d * k && i < d * d; ++i) {
                msk[i] = 1.0;
            }
        } else {
            for (int i = 0; i < d * d; i++) {
                if (i < d + d * (d + k) && i >= -k + d * (d + k))
                    msk[i] = 1.0;
                else
                    msk[i] = 0.0;
            }
        }
    return msk;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::sigmaTransform(const Ciphertext<DCRTPoly> &M) {
    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> sigma_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    int bs;
    if (squareRootIntd * squareRootIntd == 0)
        bs = squareRootIntd;
    else
        bs = round(squareRootd);

    Ciphertext<DCRTPoly> babySteps[bs];
    for (int i = 0; i < bs; i++) {
        babySteps[i] = m_cc->EvalRotate(M, i);
    }
    for (int i = 1; i < d - bs * (bs - 1); i++) {
        Plaintext pmsk =
            m_cc->MakeCKKSPackedPlaintext(generateSigmaMsk(-d + i));
        m_cc->EvalAddInPlace(sigma_M, m_cc->EvalMult(m_cc->EvalRotate(M, i - d), pmsk));
    }
    for (int i = -(bs - 1); i < bs; i++) {
        Ciphertext<DCRTPoly> tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);
        for (int j = 0; j < bs; j++) {
            auto msk = generateSigmaMsk(bs * i + j);
            msk = vectorRotate(msk, -bs * i);
            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
            m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(pmsk, babySteps[j]));
        }
        m_cc->EvalAddInPlace(sigma_M, m_cc->EvalRotate(tmp, bs * i));
    }

    return sigma_M;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::sigmaTransformHoisting(const Ciphertext<DCRTPoly> &M) {
    auto n = m_cc->GetRingDimension();
    auto m = 2*n;

    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> sigma_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    int bs;
    if (squareRootIntd * squareRootIntd == 0)
        bs = squareRootIntd;
    else
        bs = round(squareRootd);

    auto pre_M = m_cc->EvalFastRotationPrecompute(M);
    Ciphertext<DCRTPoly> babySteps[bs];
    for (int i = 0; i < bs; i++) {
        babySteps[i] = m_cc->EvalFastRotation(M, i, m, pre_M);
    }
    for (int i = 1; i < d - bs * (bs - 1); i++) {
        Plaintext pmsk =
            m_cc->MakeCKKSPackedPlaintext(generateSigmaMsk(-d + i));
        m_cc->EvalAddInPlace(sigma_M, m_cc->EvalMult(m_cc->EvalRotate(M, i - d), pmsk));
    }
    for (int i = -(bs - 1); i < bs; i++) {
        Ciphertext<DCRTPoly> tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);
        for (int j = 0; j < bs; j++) {
            auto msk = generateSigmaMsk(bs * i + j);
            msk = vectorRotate(msk, -bs * i);
            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
            m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(pmsk, babySteps[j]));
        }
        m_cc->EvalAddInPlace(sigma_M, m_cc->EvalRotate(tmp, bs * i));
    }

    return sigma_M;
}


Ciphertext<DCRTPoly> MATMULT_JKLS18::sigmaTransformLazy(const Ciphertext<DCRTPoly> &M) {
    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> sigma_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    int bs;
    if (squareRootIntd * squareRootIntd == 0)
        bs = squareRootIntd;
    else
        bs = round(squareRootd);

    // Hoisted baby steps: decompose M once, reuse for all baby-step rotations
    auto precomp_M = m_cc->EvalDirectRotatePrecompute(M);
    Ciphertext<DCRTPoly> babySteps[bs];
    for (int i = 0; i < bs; i++) {
        babySteps[i] = m_cc->EvalDirectRotate(M, i, precomp_M);
    }
    for (int i = 1; i < d - bs * (bs - 1); i++) {
        Plaintext pmsk =
            m_cc->MakeCKKSPackedPlaintext(generateSigmaMsk(-d + i));
        m_cc->EvalLazyAddInPlace(sigma_M, m_cc->EvalMult(m_cc->EvalLazyRotate(M, i - d), pmsk));
    }
    for (int i = -(bs - 1); i < bs; i++) {
        Ciphertext<DCRTPoly> tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);
        for (int j = 0; j < bs; j++) {
            auto msk = generateSigmaMsk(bs * i + j);
            msk = vectorRotate(msk, -bs * i);
            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
            m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(pmsk, babySteps[j]));
        }
        m_cc->EvalLazyAddInPlace(sigma_M, m_cc->EvalLazyRotate(tmp, bs * i));
    }
    return m_cc->EvalBatchedKS(sigma_M);
}

std::vector<double> MATMULT_JKLS18::generateTauMsk(int k) {
        std::vector<double> msk(d * d, 0);
        for (int i = k; i < d * d; i += d)
            msk[i] = 1;
        return msk;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::tauTransform(const Ciphertext<DCRTPoly> &M) {
    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> tau_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    if (squareRootIntd * squareRootIntd == d) {
        Ciphertext<DCRTPoly> babySteps[squareRootIntd];
        for (int i = 0; i < squareRootIntd; i++) {
            babySteps[i] = m_cc->EvalRotate(M, d * i);
        }

        for (int i = 0; i < squareRootIntd; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < squareRootIntd; j++) {
                auto msk = generateTauMsk(squareRootIntd * i + j);
                msk = vectorRotate(msk, -squareRootIntd * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalAddInPlace(tau_M, m_cc->EvalRotate(tmp, squareRootIntd * d * i));
        }
    } else {
        int steps = round(squareRootd);

        Ciphertext<DCRTPoly> babySteps[steps];
        for (int i = 0; i < steps; i++) {
            babySteps[i] = m_cc->EvalRotate(M, d * i);
        }

        for (int i = 0; i < d - steps * (steps - 1); i++) {
            Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(
                generateTauMsk(steps * (steps - 1) + i));
            m_cc->EvalAddInPlace(
                tau_M,
                m_cc->EvalMult(m_cc->EvalRotate(M, (steps * (steps - 1) + i) * d),
                                pmsk));
        }
        for (int i = 0; i < steps - 1; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < steps; j++) {
                auto msk = generateTauMsk(steps * i + j);
                msk = vectorRotate(msk, -steps * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalAddInPlace(tau_M, m_cc->EvalRotate(tmp, steps * d * i));
        }
    }

    return tau_M;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::tauTransformHoisting(const Ciphertext<DCRTPoly> &M) {
    auto n = m_cc->GetRingDimension();
    auto m = 2*n;

    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> tau_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    auto pre_M = m_cc->EvalFastRotationPrecompute(M);
    if (squareRootIntd * squareRootIntd == d) {
        Ciphertext<DCRTPoly> babySteps[squareRootIntd];
        for (int i = 0; i < squareRootIntd; i++) {
            babySteps[i] = m_cc->EvalFastRotation(M, d * i, m, pre_M);
        }

        for (int i = 0; i < squareRootIntd; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < squareRootIntd; j++) {
                auto msk = generateTauMsk(squareRootIntd * i + j);
                msk = vectorRotate(msk, -squareRootIntd * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalAddInPlace(tau_M, m_cc->EvalRotate(tmp, squareRootIntd * d * i));
        }
    } else {
        int steps = round(squareRootd);

        Ciphertext<DCRTPoly> babySteps[steps];
        for (int i = 0; i < steps; i++) {
            babySteps[i] = m_cc->EvalFastRotation(M, d * i, m, pre_M);
        }

        for (int i = 0; i < d - steps * (steps - 1); i++) {
            Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(
                generateTauMsk(steps * (steps - 1) + i));
            m_cc->EvalAddInPlace(
                tau_M,
                m_cc->EvalMult(m_cc->EvalRotate(M, (steps * (steps - 1) + i) * d),
                                pmsk));
        }
        for (int i = 0; i < steps - 1; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < steps; j++) {
                auto msk = generateTauMsk(steps * i + j);
                msk = vectorRotate(msk, -steps * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalAddInPlace(tau_M, m_cc->EvalRotate(tmp, steps * d * i));
        }
    }

    return tau_M;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::tauTransformLazy(const Ciphertext<DCRTPoly> &M) {
    std::vector<double> zero(d*d, 0.0);
    auto ptx_zero    = m_cc->MakeCKKSPackedPlaintext(zero);
    Ciphertext<DCRTPoly> tau_M = m_cc->Encrypt(m_PublicKey, ptx_zero);

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    // Hoisted baby steps: decompose M once, reuse for all baby-step rotations
    auto precomp_M = m_cc->EvalDirectRotatePrecompute(M);
    if (squareRootIntd * squareRootIntd == d) {
        Ciphertext<DCRTPoly> babySteps[squareRootIntd];
        for (int i = 0; i < squareRootIntd; i++) {
            babySteps[i] = m_cc->EvalDirectRotate(M, d * i, precomp_M);
        }

        for (int i = 0; i < squareRootIntd; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < squareRootIntd; j++) {
                auto msk = generateTauMsk(squareRootIntd * i + j);
                msk = vectorRotate(msk, -squareRootIntd * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalLazyAddInPlace(tau_M, m_cc->EvalLazyRotate(tmp, squareRootIntd * d * i));
        }
    } else {
        int steps = round(squareRootd);

        Ciphertext<DCRTPoly> babySteps[steps];
        for (int i = 0; i < steps; i++) {
            babySteps[i] = m_cc->EvalDirectRotate(M, d * i, precomp_M);
        }

        for (int i = 0; i < d - steps * (steps - 1); i++) {
            Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(
                generateTauMsk(steps * (steps - 1) + i));
            m_cc->EvalLazyAddInPlace(
                tau_M, m_cc->EvalMult(m_cc->EvalLazyRotate(M, (steps * (steps - 1) + i) * d), pmsk));
        }
        for (int i = 0; i < steps - 1; i++) {
            auto tmp = m_cc->Encrypt(m_PublicKey, ptx_zero);

            for (int j = 0; j < steps; j++) {
                auto msk = generateTauMsk(steps * i + j);
                msk = vectorRotate(msk, -steps * d * i);
                auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk);
                m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babySteps[j], pmsk));
            }
            m_cc->EvalLazyAddInPlace(tau_M, m_cc->EvalLazyRotate(tmp, steps * d * i));
        }
    }
    return m_cc->EvalBatchedKS(tau_M);
}

std::vector<double> MATMULT_JKLS18::generateShiftingMsk(int k) {
    std::vector<double> v(d * d, 0);
    for (int i = k; i < d * d; i += d) {
        for (int j = i; j < i + d - k; ++j) {
            v[j] = 1;
        }
    }
    return v;
}

Ciphertext<lbcrypto::DCRTPoly> MATMULT_JKLS18::columnShifting(const Ciphertext<DCRTPoly> M,
                                                int l) {
    Ciphertext<DCRTPoly> shiftResult;
    if (l == 0)
        return M;
    else {
        std::vector<double> msk = generateShiftingMsk(l);
        Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(msk);

        auto tmp = m_cc->EvalMult(pmsk, M);

        auto M_1 = m_cc->EvalRotate(m_cc->EvalSub(M, tmp), l - d);
        auto M_2 = m_cc->EvalRotate(tmp, l);

        return m_cc->EvalAdd(M_1, M_2);
    }
}
Ciphertext<lbcrypto::DCRTPoly> MATMULT_JKLS18::columnShiftingLazy(const Ciphertext<DCRTPoly> M,int l) {
    if (l == 0)
        return M;
    else {
        std::vector<double> msk = generateShiftingMsk(l);
        Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(msk);

        auto tmp = m_cc->EvalMult(pmsk, M);

        auto M_1 = m_cc->EvalLazyRotate(m_cc->EvalLazySub(M, tmp), l - d);
        auto M_2 = m_cc->EvalLazyRotate(tmp, l);

        return m_cc->EvalBatchedKS(m_cc->EvalLazyAdd(M_1, M_2));
    }
}

std::vector<double> MATMULT_JKLS18::vectorRotate(const std::vector<double> &vec,
                                    int rotateIndex) {
    if (vec.empty())
        return std::vector<double>();

    std::vector<double> result = vec;
    int n = result.size();

    if (rotateIndex > 0) // left rotation
        std::rotate(result.begin(), result.begin() + rotateIndex,
                    result.end());
    else if (rotateIndex < 0) { // right rotation
        rotateIndex += n;
        std::rotate(result.begin(), result.begin() + rotateIndex,
                    result.end());
    }
    return result;
}

// ---------------------- Baseline ----------------------
Ciphertext<DCRTPoly> MATMULT_JKLS18::eval_mult(const Ciphertext<DCRTPoly>& matA,
                                             const Ciphertext<DCRTPoly>& matB) {
    auto sigma_A = sigmaTransform(matA);
    auto tau_B = tauTransform(matB);
    auto matrixC = m_cc->EvalMultAndRelinearize(sigma_A, tau_B);

    for (int i = 1; i < d; i++) {
        auto shifted_A = columnShifting(sigma_A, i);
        tau_B = m_cc->EvalRotate(tau_B, d);
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(shifted_A, tau_B));
    }

    return matrixC;
}

// ---------------------- Hoisting ----------------------
Ciphertext<DCRTPoly> MATMULT_JKLS18::eval_mult_hoist(const Ciphertext<DCRTPoly>& matA,
                                             const Ciphertext<DCRTPoly>& matB) {
    auto n = m_cc->GetRingDimension();
    auto m = 2*n;

    auto sigma_A = sigmaTransformHoisting(matA);
    auto tau_B = tauTransformHoisting(matB);
    auto matrixC = m_cc->EvalMultAndRelinearize(sigma_A, tau_B);

    auto pre_B = m_cc->EvalFastRotationPrecompute(tau_B);
    for (int i = 1; i < d; i++) {
        auto shifted_A = columnShifting(sigma_A, i);
        auto tmp = m_cc->EvalFastRotation(tau_B, d*i, m, pre_B);
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(shifted_A, tmp));
    }
    return matrixC;
}

// ---------------------- Lazy ----------------------
Ciphertext<DCRTPoly> MATMULT_JKLS18::eval_mult_lazy(const Ciphertext<DCRTPoly>& matA,
                                                  const Ciphertext<DCRTPoly>& matB) {
    auto sigma_A = sigmaTransformLazy(matA);
    auto tau_B = tauTransformLazy(matB);
    auto matrixC = m_cc->EvalMultAndRelinearize(sigma_A, tau_B);

    // Hoisted tau_B rotations: decompose once, reuse for d*i
    auto pre_tau_B = m_cc->EvalDirectRotatePrecompute(tau_B);
    for (int i = 1; i < d; i++) {
        auto shifted_A = columnShiftingLazy(sigma_A, i);
        auto tmp = m_cc->EvalDirectRotate(tau_B, d * i, pre_tau_B);
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(shifted_A, tmp));
    }
    return matrixC;
}

// ---------------------- Double Hoisting ----------------------
// Algorithm 6 (Bossuat et al. 2020-1203): Double-hoisting BSGS.
// - Baby steps: rotations in PQ (no ModDown)
// - Inner loop: ct×pt accumulation in PQ
// - Giant step j=0: KeySwitchDownFirstElement for c₀, keep c₁ in PQ
// - Giant step j≠0: KeySwitchDown inner → Q, c₀ via AutomorphismTransform,
//   c₁ via EvalFastRotationExt (PQ, deferred ModDown)
// - Final: single KeySwitchDown on accumulated c₁, add accumulated c₀

Ciphertext<DCRTPoly> MATMULT_JKLS18::sigmaTransformDoubleHoist(const Ciphertext<DCRTPoly> &M) {
    uint32_t N = m_cc->GetRingDimension();
    uint32_t cycOrder = 2 * N;

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);
    int bs;
    if (squareRootIntd * squareRootIntd == 0)
        bs = squareRootIntd;
    else
        bs = round(squareRootd);

    // Pre-rescale M to NSD=1 if needed (mirrors AdjustLevelsAndDepthToOneInPlace in hoisting)
    auto Mr = M->Clone();
    while (Mr->GetNoiseScaleDeg() > 1) {
        m_cc->GetScheme()->ModReduceInternalInPlace(Mr, 1);
    }

    // Phase 1: Decompose once, baby steps in QP basis
    auto digits = m_cc->EvalFastRotationPrecompute(Mr);
    std::vector<Ciphertext<DCRTPoly>> babyStepsExt(bs);
    babyStepsExt[0] = m_cc->KeySwitchExt(Mr, true);
    for (int i = 1; i < bs; i++) {
        babyStepsExt[i] = m_cc->EvalFastRotationExt(Mr, i, digits, true);
    }
    auto paramsQP = babyStepsExt[0]->GetElements()[0].GetParams();

    uint32_t level = Mr->GetLevel();

    // Inner loop helper: accumulate ct×pt products in QP for a given giant step
    auto innerAccum = [&](int giantIdx) -> Ciphertext<DCRTPoly> {
        Ciphertext<DCRTPoly> accExt;
        bool first = true;
        for (int j = 0; j < bs; j++) {
            auto msk = generateSigmaMsk(bs * giantIdx + j);
            msk = vectorRotate(msk, -bs * giantIdx);
            auto ptQP = MakeQPPlaintext(m_cc, paramsQP, msk, level);
            auto term = MultExtQP(babyStepsExt[j], ptQP->GetElement<DCRTPoly>(), ptQP);
            if (first) { accExt = term; first = false; }
            else       AddExtInPlace(accExt, term);
        }
        return accExt;
    };

    // Phase 2: Algorithm 6 — giant step i=0 (no rotation, special case)
    DCRTPoly firstAccum;
    Ciphertext<DCRTPoly> resultExt;
    {
        auto accExt = innerAccum(0);
        firstAccum = m_cc->KeySwitchDownFirstElement(accExt);
        auto elements = accExt->GetElements();
        elements[0].SetValuesToZero();
        accExt->SetElements(std::move(elements));
        resultExt = accExt;
    }

    // Giant steps i ≠ 0
    for (int i = -(bs - 1); i < bs; i++) {
        if (i == 0) continue;
        auto accExt = innerAccum(i);

        // KeySwitchDown inner accumulation → Q basis
        auto inner = m_cc->KeySwitchDown(accExt);

        // c₀: rotate directly via AutomorphismTransform (no key switch needed)
        usint autoIndex = FindAutomorphismIndex2nComplex(bs * i, cycOrder);
        std::vector<usint> map(N);
        PrecomputeAutoMap(N, autoIndex, &map);
        firstAccum += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);

        // c₁: hoisted rotation in QP (deferred ModDown)
        auto innerDigits = m_cc->EvalFastRotationPrecompute(inner);
        AddExtInPlace(resultExt, m_cc->EvalFastRotationExt(inner, bs * i, innerDigits, false));
    }

    // Final: KeySwitchDown accumulated c₁, add accumulated c₀
    auto result = m_cc->KeySwitchDown(resultExt);
    {
        auto elements = result->GetElements();
        elements[0] += firstAccum;
        result->SetElements(std::move(elements));
    }

    // Remainder terms (outside BSGS grid) — standard rotation (use original M)
    for (int i = 1; i < d - bs * (bs - 1); i++) {
        Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(generateSigmaMsk(-d + i));
        m_cc->EvalAddInPlace(result, m_cc->EvalMult(m_cc->EvalRotate(M, i - d), pmsk));
    }

    return result;
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::tauTransformDoubleHoist(const Ciphertext<DCRTPoly> &M) {
    uint32_t N = m_cc->GetRingDimension();
    uint32_t cycOrder = 2 * N;

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);

    // Pre-rescale M to NSD=1 if needed
    auto Mr = M->Clone();
    while (Mr->GetNoiseScaleDeg() > 1) {
        m_cc->GetScheme()->ModReduceInternalInPlace(Mr, 1);
    }

    // Phase 1: Decompose once
    auto digits = m_cc->EvalFastRotationPrecompute(Mr);

    if (squareRootIntd * squareRootIntd == d) {
        int gs = squareRootIntd;

        // Baby steps in QP basis
        std::vector<Ciphertext<DCRTPoly>> babyStepsExt(gs);
        babyStepsExt[0] = m_cc->KeySwitchExt(Mr, true);
        for (int i = 1; i < gs; i++) {
            babyStepsExt[i] = m_cc->EvalFastRotationExt(Mr, d * i, digits, true);
        }
        auto paramsQP = babyStepsExt[0]->GetElements()[0].GetParams();
        uint32_t level = Mr->GetLevel();

        // Inner loop helper
        auto innerAccum = [&](int giantIdx) -> Ciphertext<DCRTPoly> {
            Ciphertext<DCRTPoly> accExt;
            bool first = true;
            for (int j = 0; j < gs; j++) {
                auto msk = generateTauMsk(gs * giantIdx + j);
                msk = vectorRotate(msk, -gs * d * giantIdx);
                auto ptQP = MakeQPPlaintext(m_cc, paramsQP, msk, level);
                auto term = MultExtQP(babyStepsExt[j], ptQP->GetElement<DCRTPoly>(), ptQP);
                if (first) { accExt = term; first = false; }
                else       AddExtInPlace(accExt, term);
            }
            return accExt;
        };

        // Giant step i=0: special case
        DCRTPoly firstAccum;
        Ciphertext<DCRTPoly> resultExt;
        {
            auto accExt = innerAccum(0);
            firstAccum = m_cc->KeySwitchDownFirstElement(accExt);
            auto elements = accExt->GetElements();
            elements[0].SetValuesToZero();
            accExt->SetElements(std::move(elements));
            resultExt = accExt;
        }

        // Giant steps i > 0
        for (int i = 1; i < gs; i++) {
            auto accExt = innerAccum(i);
            auto inner = m_cc->KeySwitchDown(accExt);

            usint autoIndex = FindAutomorphismIndex2nComplex(gs * d * i, cycOrder);
            std::vector<usint> map(N);
            PrecomputeAutoMap(N, autoIndex, &map);
            firstAccum += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);

            auto innerDigits = m_cc->EvalFastRotationPrecompute(inner);
            AddExtInPlace(resultExt, m_cc->EvalFastRotationExt(inner, gs * d * i, innerDigits, false));
        }

        // Final: KeySwitchDown accumulated c₁, add accumulated c₀
        auto result = m_cc->KeySwitchDown(resultExt);
        {
            auto elements = result->GetElements();
            elements[0] += firstAccum;
            result->SetElements(std::move(elements));
        }
        // No ModReduce needed: M was pre-rescaled to NSD=1, MultExtQP gives NSD=2
        return result;

    } else {
        int steps = round(squareRootd);

        // Baby steps in QP basis
        std::vector<Ciphertext<DCRTPoly>> babyStepsExt(steps);
        babyStepsExt[0] = m_cc->KeySwitchExt(Mr, true);
        for (int i = 1; i < steps; i++) {
            babyStepsExt[i] = m_cc->EvalFastRotationExt(Mr, d * i, digits, true);
        }
        auto paramsQP = babyStepsExt[0]->GetElements()[0].GetParams();
        uint32_t level = Mr->GetLevel();

        // Inner loop helper
        auto innerAccum = [&](int giantIdx) -> Ciphertext<DCRTPoly> {
            Ciphertext<DCRTPoly> accExt;
            bool first = true;
            for (int j = 0; j < steps; j++) {
                auto msk = generateTauMsk(steps * giantIdx + j);
                msk = vectorRotate(msk, -steps * d * giantIdx);
                auto ptQP = MakeQPPlaintext(m_cc, paramsQP, msk, level);
                auto term = MultExtQP(babyStepsExt[j], ptQP->GetElement<DCRTPoly>(), ptQP);
                if (first) { accExt = term; first = false; }
                else       AddExtInPlace(accExt, term);
            }
            return accExt;
        };

        // Giant step i=0: special case
        DCRTPoly firstAccum;
        Ciphertext<DCRTPoly> resultExt;
        {
            auto accExt = innerAccum(0);
            firstAccum = m_cc->KeySwitchDownFirstElement(accExt);
            auto elements = accExt->GetElements();
            elements[0].SetValuesToZero();
            accExt->SetElements(std::move(elements));
            resultExt = accExt;
        }

        // Giant steps i > 0
        for (int i = 1; i < steps - 1; i++) {
            auto accExt = innerAccum(i);
            auto inner = m_cc->KeySwitchDown(accExt);

            usint autoIndex = FindAutomorphismIndex2nComplex(steps * d * i, cycOrder);
            std::vector<usint> map(N);
            PrecomputeAutoMap(N, autoIndex, &map);
            firstAccum += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);

            auto innerDigits = m_cc->EvalFastRotationPrecompute(inner);
            AddExtInPlace(resultExt, m_cc->EvalFastRotationExt(inner, steps * d * i, innerDigits, false));
        }

        // Final: KeySwitchDown accumulated c₁, add accumulated c₀
        auto result = m_cc->KeySwitchDown(resultExt);
        {
            auto elements = result->GetElements();
            elements[0] += firstAccum;
            result->SetElements(std::move(elements));
        }
        // No ModReduce needed: M was pre-rescaled to NSD=1, MultExtQP gives NSD=2

        // Remainder terms (outside BSGS grid) — standard rotation
        for (int i = 0; i < d - steps * (steps - 1); i++) {
            Plaintext pmsk = m_cc->MakeCKKSPackedPlaintext(
                generateTauMsk(steps * (steps - 1) + i));
            m_cc->EvalAddInPlace(
                result,
                m_cc->EvalMult(m_cc->EvalRotate(M, (steps * (steps - 1) + i) * d),
                                pmsk));
        }

        return result;
    }
}

Ciphertext<DCRTPoly> MATMULT_JKLS18::eval_mult_double_hoist(const Ciphertext<DCRTPoly>& matA,
                                                             const Ciphertext<DCRTPoly>& matB) {
    auto n = m_cc->GetRingDimension();
    auto m = 2*n;

    auto sigma_A = sigmaTransformDoubleHoist(matA);
    auto tau_B = tauTransformDoubleHoist(matB);
    auto matrixC = m_cc->EvalMultAndRelinearize(sigma_A, tau_B);

    // Main loop: same as hoisting variant (tau_B rotations via hoisting)
    auto pre_B = m_cc->EvalFastRotationPrecompute(tau_B);
    for (int i = 1; i < d; i++) {
        auto shifted_A = columnShifting(sigma_A, i);
        auto tmp = m_cc->EvalFastRotation(tau_B, d*i, m, pre_B);
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(shifted_A, tmp));
    }
    return matrixC;
}


// ---------------------- Plan functions ----------------------
// Note: Caller must call rk.begin() with appropriate slotCount (ringDim/2) before calling these functions
void MATMULT_JKLS18::eval_mult_plan(RotationKeyCollector& rk) const {
    double squareRootd = sqrt(static_cast<double>(d));
    int bs = static_cast<int>(round(squareRootd));

    // sigmaTransform rotations
    for (int i = 0; i < bs; i++) rk.observe(i);
    for (int i = 1; i < d - bs * (bs - 1); i++) rk.observe(i - d);
    for (int i = -(bs - 1); i < bs; i++) rk.observe(bs * i);

    // tauTransform rotations
    int squareRootIntd = static_cast<int>(squareRootd);
    if (squareRootIntd * squareRootIntd == d) {
        for (int i = 0; i < squareRootIntd; i++) rk.observe(d * i);
        for (int i = 0; i < squareRootIntd; i++) rk.observe(squareRootIntd * d * i);
    } else {
        int steps = bs;
        for (int i = 0; i < steps; i++) rk.observe(d * i);
        for (int i = 0; i < d - steps * (steps - 1); i++) rk.observe((steps * (steps - 1) + i) * d);
        for (int i = 0; i < steps - 1; i++) rk.observe(steps * d * i);
    }

    // columnShifting rotations
    for (int l = 1; l < d; l++) {
        rk.observe(l - d);
        rk.observe(l);
    }

    // Main loop: tau_B rotation by d
    rk.observe(d);
}

void MATMULT_JKLS18::eval_mult_hoist_plan(RotationKeyCollector& rk) const {
    // Start with baseline rotation pattern
    eval_mult_plan(rk);

    // Hoisting-specific: main loop uses EvalFastRotation(tau_B, d*i)
    // instead of iterative EvalRotate(tau_B, d), so we need d*i for i in [1, d)
    for (int i = 1; i < d; i++) {
        rk.observe(d * i);
    }
}

void MATMULT_JKLS18::eval_mult_double_hoist_plan(RotationKeyCollector& rk) const {
    // Double hoisting uses the same rotation keys as hoisting
    // (baby steps are in PQ basis so they don't need automorphism keys,
    //  but giant step rotations and remainder terms need keys)
    eval_mult_hoist_plan(rk);
}

void MATMULT_JKLS18::eval_mult_lazy_plan(RotationKeyCollectorLazy& rk) const {
    const int num_slots = d * d;
    rk.begin(num_slots);

    // Helper to create zero ciphertext
    auto makeZeroCT = [&]() {
        std::vector<double> zero(num_slots, 0.0);
        auto pt = m_cc->MakeCKKSPackedPlaintext(zero);
        return m_cc->Encrypt(m_PublicKey, pt);
    };

    double squareRootd = sqrt(static_cast<double>(d));
    int squareRootIntd = static_cast<int>(squareRootd);
    int bs = (squareRootIntd * squareRootIntd == 0) ? squareRootIntd : static_cast<int>(round(squareRootd));

    // Simulate sigmaTransformLazy
    {
        auto sigma_M = makeZeroCT();
        auto M = makeZeroCT();

        // babySteps - each one has individual batched KS
        for (int i = 0; i < bs; i++) {
            auto ct = m_cc->EvalLazyRotate(M, i);
            rk.observeAutoIndices(ct->GetElementKeyIndexVector());  // EvalBatchedKS point
        }

        // Second loop - accumulates lazy rotations into sigma_M
        for (int i = 1; i < d - bs * (bs - 1); i++) {
            // EvalMult with plaintext doesn't affect automorphism
            // We just need to accumulate the lazy rotation
            sigma_M = m_cc->EvalLazyAdd(sigma_M, m_cc->EvalLazyRotate(M, i - d));
        }

        // Third loop - accumulates more lazy rotations into sigma_M
        for (int i = -(bs - 1); i < bs; i++) {
            // tmp stays clean (already batched KS'ed in real code)
            // We just add the rotation to sigma_M
            sigma_M = m_cc->EvalLazyAdd(sigma_M, m_cc->EvalLazyRotate(makeZeroCT(), bs * i));
        }

        // Final EvalBatchedKS on sigma_M - observe all accumulated automorphisms
        rk.observeAutoIndices(sigma_M->GetElementKeyIndexVector());
    }

    // Simulate tauTransformLazy
    {
        auto tau_M = makeZeroCT();
        auto M = makeZeroCT();

        if (squareRootIntd * squareRootIntd == d) {
            // babySteps - each has individual batched KS
            for (int i = 0; i < squareRootIntd; i++) {
                auto ct = m_cc->EvalLazyRotate(M, d * i);
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());  // EvalBatchedKS point
            }

            // Accumulate lazy rotations into tau_M
            for (int i = 0; i < squareRootIntd; i++) {
                tau_M = m_cc->EvalLazyAdd(tau_M, m_cc->EvalLazyRotate(makeZeroCT(), squareRootIntd * d * i));
            }
        } else {
            int steps = bs;

            // babySteps - each has individual batched KS
            for (int i = 0; i < steps; i++) {
                auto ct = m_cc->EvalLazyRotate(M, d * i);
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());  // EvalBatchedKS point
            }

            // Accumulate lazy rotations
            for (int i = 0; i < d - steps * (steps - 1); i++) {
                tau_M = m_cc->EvalLazyAdd(tau_M, m_cc->EvalLazyRotate(M, (steps * (steps - 1) + i) * d));
            }

            for (int i = 0; i < steps - 1; i++) {
                tau_M = m_cc->EvalLazyAdd(tau_M, m_cc->EvalLazyRotate(makeZeroCT(), steps * d * i));
            }
        }

        // Final EvalBatchedKS on tau_M - observe all accumulated automorphisms
        rk.observeAutoIndices(tau_M->GetElementKeyIndexVector());
    }

    // Simulate eval_mult_lazy main loop
    {
        auto sigma_A = makeZeroCT();
        auto tau_B = makeZeroCT();

        for (int i = 1; i < d; i++) {
            // columnShiftingLazy simulation
            // Real code: EvalLazySub, then two rotations, then EvalLazyAdd, then EvalBatchedKS
            {
                auto M = sigma_A;  // Use sigma_A as source
                auto temp1 = makeZeroCT();  // Simulates EvalLazySub result
                auto ct1 = m_cc->EvalLazyRotate(temp1, i - d);
                auto temp2 = makeZeroCT();  // Simulates the masked part
                auto ct2 = m_cc->EvalLazyRotate(temp2, i);
                auto shifted = m_cc->EvalLazyAdd(ct1, ct2);
                rk.observeAutoIndices(shifted->GetElementKeyIndexVector());  // EvalBatchedKS point
            }

            // tau_B = EvalBatchedKS(EvalLazyRotate(tau_B, d))
            // Since tau_B accumulates from previous iterations
            tau_B = m_cc->EvalLazyRotate(tau_B, d);
            rk.observeAutoIndices(tau_B->GetElementKeyIndexVector());
        }
    }
}
