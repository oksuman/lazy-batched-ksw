#pragma once

#include "openfhe.h"
#include <cstdint>
#include <vector>

using namespace lbcrypto;

// Template-free, runtime-dimension implementation
class MATINV_KCL25 {
public:
    MATINV_KCL25(const CryptoContext<DCRTPoly>& cc,
                 const PublicKey<DCRTPoly>& pk,
                 int dim,
                 int iterations,
                 int multDepth);

    Ciphertext<DCRTPoly> eval_inverse(const Ciphertext<DCRTPoly>& M);
    Ciphertext<DCRTPoly> eval_inverse_lazy(const Ciphertext<DCRTPoly>& M);

    // Debug versions that decrypt intermediate values
    Ciphertext<DCRTPoly> eval_inverse_debug(const Ciphertext<DCRTPoly>& M,
                                            const PrivateKey<DCRTPoly>& sk);
    Ciphertext<DCRTPoly> eval_inverse_lazy_debug(const Ciphertext<DCRTPoly>& M,
                                                 const PrivateKey<DCRTPoly>& sk);

private:
    CryptoContext<DCRTPoly> m_cc;
    PublicKey<DCRTPoly>     m_pk;

    int d     = 0;   // dimension
    int r     = 0;   // iteration count
    int depth = 0;   // multDepth (passed from caller)

    // batching parameters
    int max_batch = 0;
    int s         = 0;
    int B         = 0;
    int num_slots = 0;

    // bsgs / precomp params
    int ng = -1;
    int nb = -1;
    int np = -1;

    Ciphertext<DCRTPoly> m_zeroCache;

private:
    void init_params_or_throw();
    Ciphertext<DCRTPoly> createZeroCache() const;
    Ciphertext<DCRTPoly> zeroClone() const;

    std::vector<double> vectorRotate(const std::vector<double>& vec, int rotateIndex) const;

    std::vector<double> generateTransposeMsk(int k) const;
    Ciphertext<DCRTPoly> eval_transpose(Ciphertext<DCRTPoly> M) const;

    Ciphertext<DCRTPoly> eval_trace(Ciphertext<DCRTPoly> M, int batchSize) const;

    std::vector<double> generateMaskVector(int batch_size, int k) const;
    std::vector<double> genDiagVector(int k, int diag_index) const;
    std::vector<double> genBatchDiagVector(int s_local, int k, int diag_index) const;

    Ciphertext<DCRTPoly> vecRotsOpt(const std::vector<Ciphertext<DCRTPoly>>& matrixM, int is) const;

    Ciphertext<DCRTPoly> eval_mult(const Ciphertext<DCRTPoly>& matrixA,
                                   const Ciphertext<DCRTPoly>& matrixB) const;

    Ciphertext<DCRTPoly> eval_mult_lazy(const Ciphertext<DCRTPoly>& matrixA,
                                        const Ciphertext<DCRTPoly>& matrixB) const;

    std::vector<double> initializeIdentityMatrix(int dim) const;
};
