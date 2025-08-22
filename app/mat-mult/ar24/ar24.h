#pragma once

#include <openfhe.h>
#include <memory>
#include <vector>
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

using namespace lbcrypto;

class MATMULT_AR24 {
public:
    MATMULT_AR24(const CryptoContext<DCRTPoly>& cc,
                 const PublicKey<DCRTPoly>& pk,
                 int dim);
    MATMULT_AR24(int dim);

    std::vector<double> generatePhiMsk(int k);
    std::vector<double> generatePsiMsk(int k);

    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& matA,
                                   const Ciphertext<DCRTPoly>& matB);

    Ciphertext<DCRTPoly> eval_mult_lazy(const Ciphertext<DCRTPoly>& matA,
                                        const Ciphertext<DCRTPoly>& matB);


    void eval_mult_plan(RotationKeyCollector& rk) const;
    void eval_mult_lazy_plan(RotationKeyCollectorLazy& rk) const;
    std::vector<int32_t> plan_required_auto_indices() const;

private:
    CryptoContext<DCRTPoly> m_cc;
    PublicKey<DCRTPoly>     m_PublicKey;
    int d;          // matrix dimension
    int max_batch;
    int B;
    int s;
};
