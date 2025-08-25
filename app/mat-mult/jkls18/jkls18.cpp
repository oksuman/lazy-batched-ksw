#include <algorithm>
#include <cmath>
#include "jkls18.h"

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

    Ciphertext<DCRTPoly> babySteps[bs];
    for (int i = 0; i < bs; i++) {
        babySteps[i] = m_cc->EvalBatchedKS(m_cc->EvalLazyRotate(M, i));
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

    if (squareRootIntd * squareRootIntd == d) {
        Ciphertext<DCRTPoly> babySteps[squareRootIntd];
        for (int i = 0; i < squareRootIntd; i++) {
            babySteps[i] = m_cc->EvalBatchedKS(m_cc->EvalLazyRotate(M, d * i));
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
            babySteps[i] = m_cc->EvalBatchedKS(m_cc->EvalLazyRotate(M, d * i));
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
    
    for (int i = 1; i < d; i++) {
        auto shifted_A = columnShiftingLazy(sigma_A, i);
        tau_B = m_cc->EvalBatchedKS(m_cc->EvalLazyRotate(tau_B, d));
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(shifted_A, tau_B));
    }
    return matrixC;
}
