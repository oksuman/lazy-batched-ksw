#pragma once

#include <openfhe.h>
#include <memory>
#include <vector>
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

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

    // Double hoisting variants
    Ciphertext<DCRTPoly> sigmaTransformDoubleHoist(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> tauTransformDoubleHoist(const Ciphertext<DCRTPoly> &M);
    Ciphertext<DCRTPoly> eval_mult_double_hoist(const Ciphertext<DCRTPoly>& matA,
                                                const Ciphertext<DCRTPoly>& matB);

    // Plan functions to collect rotation indices
    void eval_mult_plan(RotationKeyCollector& rk) const;
    void eval_mult_hoist_plan(RotationKeyCollector& rk) const;
    void eval_mult_double_hoist_plan(RotationKeyCollector& rk) const;
    void eval_mult_lazy_plan(RotationKeyCollectorLazy& rk) const;

private:
    CryptoContext<DCRTPoly> m_cc;
    PublicKey<DCRTPoly>     m_PublicKey;
    int d;          // matrix dimension
};
