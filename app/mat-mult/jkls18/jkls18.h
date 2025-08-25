#pragma once

#include <openfhe.h>
#include <memory>
#include <vector>

using namespace lbcrypto;

class MATMULT_JKLS18 {
public:
    MATMULT_JKLS18(const CryptoContext<DCRTPoly>& cc,
                 const PublicKey<DCRTPoly>& pk,
                 int dim);
    MATMULT_JKLS18(int dim);

    std::vector<double> generateSigmaMsk(int k);
    Ciphertext<DCRTPoly> sigmaTransform(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> sigmaTransformHoisting(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> sigmaTransformLazy(const Ciphertext<DCRTPoly> &M);

    std::vector<double> generateTauMsk(int k);
    Ciphertext<DCRTPoly> tauTransform(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> tauTransformHoisting(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> tauTransformLazy(const Ciphertext<DCRTPoly> &M);

    std::vector<double> generateShiftingMsk(int k);
    Ciphertext<lbcrypto::DCRTPoly> columnShifting(const Ciphertext<DCRTPoly> M, int l);
    Ciphertext<lbcrypto::DCRTPoly> columnShiftingLazy(const Ciphertext<DCRTPoly> M, int l);

    std::vector<double> vectorRotate(const std::vector<double> &vec, int rotateIndex);

    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& matA,
                                   const Ciphertext<DCRTPoly>& matB);
    Ciphertext<DCRTPoly> eval_mult_hoist(const Ciphertext<DCRTPoly>& matA,
                                   const Ciphertext<DCRTPoly>& matB);

    Ciphertext<DCRTPoly> eval_mult_lazy(const Ciphertext<DCRTPoly>& matA,
                                        const Ciphertext<DCRTPoly>& matB);

private:
    CryptoContext<DCRTPoly> m_cc;
    PublicKey<DCRTPoly>     m_PublicKey;
    int d;          // matrix dimension
};
