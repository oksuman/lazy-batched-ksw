#include "kcl25.h"
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

MATINV_KCL25::MATINV_KCL25(const CryptoContext<DCRTPoly>& cc,
                           const PublicKey<DCRTPoly>& pk,
                           int dim,
                           int iterations,
                           int multDepth)
    : m_cc(cc), m_pk(pk), d(dim), r(iterations), depth(multDepth) {
    init_params_or_throw();
    m_zeroCache = createZeroCache();
}

void MATINV_KCL25::init_params_or_throw() {
    if (!(d == 4 || d == 8 || d == 16 || d == 32 || d == 64))
        throw std::runtime_error("Unsupported dimension d=" + std::to_string(d));
    if (r <= 0)
        throw std::runtime_error("Invalid iterations r=" + std::to_string(r));
    if (depth <= 0)
        throw std::runtime_error("Invalid multDepth=" + std::to_string(depth));

    max_batch = static_cast<int>(m_cc->GetRingDimension() / 2);
    if (max_batch <= 0)
        throw std::runtime_error("Invalid ring dimension");

    s = std::min(max_batch / d / d, d);
    if (s <= 0) s = 1;

    B         = d / s;
    num_slots = s * d * d;

    // np depends on d (same as original)
    switch (d) {
        case 4:  np = 2; break;
        case 8:  np = 2; break;
        case 16: np = 4; break;
        case 32: np = 4; break;
        case 64: np = 8; break;
        default: np = 2; break;
    }

    // ng/nb depends on max_batch (same table as original)
    if (max_batch == (1 << 13)) {
        switch (d) {
            case 4:  ng = 2; nb = 2;  break;
            case 8:  ng = 2; nb = 4;  break;
            case 16: ng = 4; nb = 4;  np = 4; break;
            case 32: ng = 2; nb = 16; break;
            case 64: ng = 1; nb = 64; break;
            default: ng = -1; nb = -1; break;
        }
    }
    else if (max_batch == (1 << 15)) {
        switch (d) {
            case 4:  ng = 2; nb = 2;  break;
            case 8:  ng = 2; nb = 4;  break;
            case 16: ng = 4; nb = 4;  break;
            case 32: ng = 4; nb = 8;  break;
            case 64: ng = 2; nb = 32; break;
            default: ng = -1; nb = -1; break;
        }
    }
    else if (max_batch == (1 << 16)) {
        switch (d) {
            case 4:  ng = 2; nb = 2;  break;
            case 8:  ng = 2; nb = 4;  break;
            case 16: ng = 4; nb = 4;  break;
            case 32: ng = 4; nb = 8;  break;
            case 64: ng = 4; nb = 16; break;
            default: ng = -1; nb = -1; break;
        }
    }
    else {
        ng = -1;
        nb = -1;
    }

    if (ng <= 0 || nb <= 0)
        throw std::runtime_error("Unsupported (ringDim/2) max_batch=" + std::to_string(max_batch));
}

Ciphertext<DCRTPoly> MATINV_KCL25::createZeroCache() const {
    // Match newCol: d*d sized zero vector
    std::vector<double> zeroVec(static_cast<size_t>(d * d), 0.0);
    auto pt = m_cc->MakeCKKSPackedPlaintext(zeroVec);
    return m_cc->Encrypt(m_pk, pt);
}

Ciphertext<DCRTPoly> MATINV_KCL25::zeroClone() const {
    auto z = m_zeroCache->Clone();
    z->SetSlots(num_slots);
    return z;
}

std::vector<double> MATINV_KCL25::vectorRotate(const std::vector<double>& vec, int rotateIndex) const {
    if (vec.empty()) return std::vector<double>();
    std::vector<double> result = vec;
    int n = static_cast<int>(result.size());

    if (rotateIndex > 0) {
        std::rotate(result.begin(), result.begin() + rotateIndex, result.end());
    } else if (rotateIndex < 0) {
        rotateIndex += n;
        std::rotate(result.begin(), result.begin() + rotateIndex, result.end());
    }
    return result;
}

std::vector<double> MATINV_KCL25::generateTransposeMsk(int k) const {
    std::set<int> indices;
    if (k >= 0) {
        for (int j = 0; j < d - k; j++) indices.insert((d + 1) * j + k);
    } else {
        for (int j = -k; j < d; j++) indices.insert((d + 1) * j + k);
    }
    std::vector<double> msk(static_cast<size_t>(d * d), 0.0);
    for (int idx : indices) msk[static_cast<size_t>(idx)] = 1.0;
    return msk;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_transpose(Ciphertext<DCRTPoly> M) const {
    auto p0 = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(0));
    auto M_transposed = m_cc->EvalMult(M, p0);

    for (int i = 1; i < d; i++) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalRotate(M, (d - 1) * i), pi));
    }
    for (int i = -1; i > -d; i--) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalRotate(M, (d - 1) * i), pi));
    }
    return M_transposed;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_transpose_lazy(Ciphertext<DCRTPoly> M) const {
    auto p0 = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(0));
    auto M_transposed = m_cc->EvalMult(M, p0);

    for (int i = 1; i < d; i++) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalLazyAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalLazyRotate(M, (d - 1) * i), pi));
    }
    for (int i = -1; i > -d; i--) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalLazyAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalLazyRotate(M, (d - 1) * i), pi));
    }
    M_transposed = m_cc->EvalBatchedKS(M_transposed);
    return M_transposed;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_trace(Ciphertext<DCRTPoly> M, int batchSize) const {
    std::vector<double> msk(static_cast<size_t>(batchSize), 0.0);
    for (int i = 0; i < d * d; i += (d + 1)) {
        msk[static_cast<size_t>(i)] = 1.0;
    }

    auto trace = m_cc->EvalMult(M, m_cc->MakeCKKSPackedPlaintext(msk));

    // Sum across batchSize using rotations
    int logB = static_cast<int>(std::log2(static_cast<double>(batchSize)));
    for (int i = 1; i <= logB; i++) {
        m_cc->EvalAddInPlace(trace, m_cc->EvalRotate(trace, batchSize / (1 << i)));
    }
    return trace;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_trace_lazy(Ciphertext<DCRTPoly> M, int batchSize) const {
    std::vector<double> msk(static_cast<size_t>(batchSize), 0.0);
    for (int i = 0; i < d * d; i += (d + 1)) {
        msk[static_cast<size_t>(i)] = 1.0;
    }

    auto trace = m_cc->EvalMult(M, m_cc->MakeCKKSPackedPlaintext(msk));

    // Sum across batchSize using rotations
    int logB = static_cast<int>(std::log2(static_cast<double>(batchSize)));
    for (int i = 1; i <= logB; i++) {
        m_cc->EvalAddInPlace(trace, m_cc->EvalDirectRotate(trace, batchSize / (1 << i)));
    }
    return trace;
}

std::vector<double> MATINV_KCL25::generateMaskVector(int batch_size, int k) const {
    std::vector<double> result(batch_size, 0.0);
    for (int i = k * d * d; i < (k + 1) * d * d; ++i) {
        result[i] = 1.0;
    }
    return result;
}

std::vector<double> MATINV_KCL25::genDiagVector(int k, int diag_index) const  {
    std::vector<double> result(d * d, 0.0);

    if (diag_index < 1 || diag_index > d * d ||
        (diag_index > d && diag_index < d * d - (d - 1))) {
        return result;
    }

    for (int i = 0; i < d; ++i) {
        result[i * d + ((i + k) % d)] = 1.0;
    }

    int rotation = 0;
    bool right_rotation = false;

    if (diag_index <= d) {
        rotation = diag_index - 1;
    } else {
        right_rotation = true;
        rotation = d * d - diag_index + 1;
    }

    if (rotation > 0) {
        for (int i = 0; i < rotation; ++i) {
            for (int j = 0; j < d; ++j) {
                if (right_rotation) {
                    result[j * d + (d - 1 - i)] = 0.0;
                } else {
                    result[j * d + i] = 0.0;
                }
            }
        }
    }

    std::vector<double> rotated(d * d, 0.0);
    for (int i = 0; i < d * d; ++i) {
        int new_pos;
        if (right_rotation) {
            new_pos = (i + rotation) % d + (i / d) * d;
        } else {
            new_pos = (i + d - rotation) % d + (i / d) * d;
        }
        rotated[new_pos] = result[i];
    }

    return rotated;
}


std::vector<double> MATINV_KCL25::genBatchDiagVector(int s, int k, int diag_index) const {
    std::vector<double> result;
    result.reserve(d * d * s);

    for (int i = 0; i < s; ++i) {
        std::vector<double> diag_vector = genDiagVector(k + i, diag_index);
        result.insert(result.end(), diag_vector.begin(), diag_vector.end());
    }

    return result;
}

Ciphertext<DCRTPoly>
MATINV_KCL25::vecRotsOpt(const std::vector<Ciphertext<DCRTPoly>>& matrixM, int is) const {
    auto rotsM = zeroClone();

    // if s < np => loop won't run => returns zero
    for (int j = 0; j < s / np; j++) {
        auto T = zeroClone();

        for (int i = 0; i < np; i++) {
            auto msk = generateMaskVector(num_slots, np * j + i);
            msk = vectorRotate(msk, -is * d * s - j * d * np);

            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk, 1, 0, nullptr, num_slots);
            m_cc->EvalAddInPlace(T, m_cc->EvalMult(matrixM[static_cast<size_t>(i)], pmsk));
        }
        m_cc->EvalAddInPlace(rotsM, m_cc->EvalRotate(T, is * d * s + j * d * np));
    }

    return rotsM;
}

Ciphertext<DCRTPoly>
MATINV_KCL25::vecRotsOptLazy(const std::vector<Ciphertext<DCRTPoly>>& matrixM, int is) const {
    auto rotsM = zeroClone();

    // if s < np => loop won't run => returns zero
    for (int j = 0; j < s / np; j++) {
        auto T = zeroClone();

        for (int i = 0; i < np; i++) {
            auto msk = generateMaskVector(num_slots, np * j + i);
            msk = vectorRotate(msk, -is * d * s - j * d * np);

            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk, 1, 0, nullptr, num_slots);
            m_cc->EvalAddInPlace(T, m_cc->EvalMult(matrixM[static_cast<size_t>(i)], pmsk));
        }
        m_cc->EvalLazyAddInPlace(rotsM, m_cc->EvalLazyRotate(T, is * d * s + j * d * np));
    }
    rotsM = m_cc->EvalBatchedKS(rotsM);
    return rotsM;
}

Ciphertext<DCRTPoly>
MATINV_KCL25::eval_mult(const Ciphertext<DCRTPoly>& matrixA,
                        const Ciphertext<DCRTPoly>& matrixB) const {
    auto matrixC = zeroClone();

    std::vector<Ciphertext<DCRTPoly>> babyStepsOfA(static_cast<size_t>(nb));
    for (int i = 0; i < nb; i++) {
        babyStepsOfA[static_cast<size_t>(i)] = m_cc->EvalRotate(matrixA, i);
    }

    std::vector<Ciphertext<DCRTPoly>> babyStepsOfB;
    if (s >= np) {
        babyStepsOfB.reserve(static_cast<size_t>(np));
        for (int i = 0; i < np; i++) {
            auto t = m_cc->EvalRotate(matrixB, i * d);
            t->SetSlots(num_slots);
            babyStepsOfB.push_back(t);
        }
    }

    for (int i = 0; i < B; i++) {
        auto batched_rotations_B = vecRotsOpt(babyStepsOfB, i);

        auto diagA = zeroClone();
        for (int k = -ng; k < ng; k++) {
            if (k < 0) {
                auto tmp = zeroClone();
                int babyStep = (k == -ng) ? 1 : 0;

                int jStart = d * d + k * nb + 1 + babyStep;
                int jEnd   = d * d + (k + 1) * nb;

                for (int j = jStart; j <= jEnd; j++) {
                    auto rotated_plain_vec = vectorRotate(genBatchDiagVector(s, i * s, j), -k * nb);
                    auto pt = m_cc->MakeCKKSPackedPlaintext(rotated_plain_vec, 1, 0, nullptr, num_slots);
                    m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babyStepsOfA[static_cast<size_t>(babyStep)], pt));
                    babyStep++;
                }
                m_cc->EvalAddInPlace(diagA, m_cc->EvalRotate(tmp, k * nb));
            } else {
                auto tmp = zeroClone();
                int babyStep = 0;

                int jStart = k * nb + 1;
                int jEnd   = (k + 1) * nb;

                for (int j = jStart; j <= jEnd; j++) {
                    auto rotated_plain_vec = vectorRotate(genBatchDiagVector(s, i * s, j), -k * nb);
                    auto pt = m_cc->MakeCKKSPackedPlaintext(rotated_plain_vec, 1, 0, nullptr, num_slots);
                    m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babyStepsOfA[static_cast<size_t>(babyStep)], pt));
                    babyStep++;
                }
                m_cc->EvalAddInPlace(diagA, m_cc->EvalRotate(tmp, k * nb));
            }
        }

        m_cc->EvalAddInPlace(matrixC, m_cc->EvalMult(diagA, batched_rotations_B));
    }

    // Final aggregation: s is power-of-two in these settings
    int logS = static_cast<int>(std::log2(static_cast<double>(s)));
    for (int i = 1; i <= logS; i++) {
        m_cc->EvalAddInPlace(matrixC, m_cc->EvalRotate(matrixC, num_slots / (1 << i)));
    }

    matrixC->SetSlots(d * d);
    return matrixC;
}

Ciphertext<DCRTPoly>
MATINV_KCL25::eval_mult_lazy(const Ciphertext<DCRTPoly>& matrixA,
                             const Ciphertext<DCRTPoly>& matrixB) const {
    auto matrixC = zeroClone();

    std::vector<Ciphertext<DCRTPoly>> babyStepsOfA(static_cast<size_t>(nb));
    for (int i = 0; i < nb; i++) {
        babyStepsOfA[static_cast<size_t>(i)] = m_cc->EvalDirectRotate(matrixA, i);
    }

    std::vector<Ciphertext<DCRTPoly>> babyStepsOfB;
    if (s >= np) {
        babyStepsOfB.reserve(static_cast<size_t>(np));
        for (int i = 0; i < np; i++) {
            auto t = m_cc->EvalDirectRotate(matrixB, i * d);
            t->SetSlots(num_slots);
            babyStepsOfB.push_back(t);
        }
    }

    for (int i = 0; i < B; i++) {
        auto batched_rotations_B = vecRotsOptLazy(babyStepsOfB, i);

        auto diagA = zeroClone();
        for (int k = -ng; k < ng; k++) {
            if (k < 0) {
                auto tmp = zeroClone();
                int babyStep = (k == -ng) ? 1 : 0;

                int jStart = d * d + k * nb + 1 + babyStep;
                int jEnd   = d * d + (k + 1) * nb;

                for (int j = jStart; j <= jEnd; j++) {
                    auto rotated_plain_vec = vectorRotate(genBatchDiagVector(s, i * s, j), -k * nb);
                    auto pt = m_cc->MakeCKKSPackedPlaintext(rotated_plain_vec, 1, 0, nullptr, num_slots);
                    m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babyStepsOfA[static_cast<size_t>(babyStep)], pt));
                    babyStep++;
                }
                m_cc->EvalLazyAddInPlace(diagA, m_cc->EvalLazyRotate(tmp, k * nb));
            } else {
                auto tmp = zeroClone();
                int babyStep = 0;

                int jStart = k * nb + 1;
                int jEnd   = (k + 1) * nb;

                for (int j = jStart; j <= jEnd; j++) {
                    auto rotated_plain_vec = vectorRotate(genBatchDiagVector(s, i * s, j), -k * nb);
                    auto pt = m_cc->MakeCKKSPackedPlaintext(rotated_plain_vec, 1, 0, nullptr, num_slots);
                    m_cc->EvalAddInPlace(tmp, m_cc->EvalMult(babyStepsOfA[static_cast<size_t>(babyStep)], pt));
                    babyStep++;
                }
                m_cc->EvalLazyAddInPlace(diagA, m_cc->EvalLazyRotate(tmp, k * nb));
            }
        }
        diagA = m_cc->EvalBatchedKS(diagA);
        m_cc->EvalAddInPlace(matrixC, m_cc->EvalMult(diagA, batched_rotations_B));
    }

    // Final aggregation: s is power-of-two in these settings
    int logS = static_cast<int>(std::log2(static_cast<double>(s)));
    for (int i = 1; i <= logS; i++) {
        m_cc->EvalAddInPlace(matrixC, m_cc->EvalDirectRotate(matrixC, num_slots / (1 << i)));
    }

    matrixC->SetSlots(d * d);
    return matrixC;
}

std::vector<double> MATINV_KCL25::initializeIdentityMatrix(int dim) const {
    std::vector<double> identity(static_cast<size_t>(dim * dim), 0.0);
    for (int i = 0; i < dim; i++) identity[static_cast<size_t>(i * dim + i)] = 1.0;
    return identity;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_inverse(const Ciphertext<DCRTPoly>& M) {
    auto vI = initializeIdentityMatrix(d);
    Plaintext pI = m_cc->MakeCKKSPackedPlaintext(vI);

    auto M_transposed  = eval_transpose(M);
    auto MM_transposed = eval_mult(M, M_transposed);

    auto trace = eval_trace(MM_transposed, d * d);
    auto trace_reciprocal = m_cc->EvalDivide(trace, (d * d) / 3 - d, (d * d) / 3 + d, 5);

    auto Y     = m_cc->EvalMultAndRelinearize(M_transposed, trace_reciprocal);
    auto A_bar = m_cc->EvalSub(pI, m_cc->EvalMultAndRelinearize(MM_transposed, trace_reciprocal));

    for (int i = 0; i < r - 1; i++) {
        if (d >= 8 && static_cast<int>(Y->GetLevel()) >= depth - 2) {
            A_bar = m_cc->EvalBootstrap(A_bar);
            Y     = m_cc->EvalBootstrap(Y);
        }
        Y     = eval_mult(Y, m_cc->EvalAdd(pI, A_bar));
        A_bar = eval_mult(A_bar, A_bar);
    }
    Y = eval_mult(Y, m_cc->EvalAdd(pI, A_bar));

    return Y;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_inverse_lazy(const Ciphertext<DCRTPoly>& M) {
    auto vI = initializeIdentityMatrix(d);
    Plaintext pI = m_cc->MakeCKKSPackedPlaintext(vI);

    auto M_transposed  = eval_transpose_lazy(M);
    auto MM_transposed = eval_mult_lazy(M, M_transposed);

    auto trace = eval_trace_lazy(MM_transposed, d * d);
    auto trace_reciprocal = m_cc->EvalDivide(trace, (d * d) / 3 - d, (d * d) / 3 + d, 5);

    auto Y     = m_cc->EvalMultAndRelinearize(M_transposed, trace_reciprocal);
    auto A_bar = m_cc->EvalSub(pI, m_cc->EvalMultAndRelinearize(MM_transposed, trace_reciprocal));

    for (int i = 0; i < r - 1; i++) {
        if (d >= 8 && static_cast<int>(Y->GetLevel()) >= depth - 2) {
            A_bar = m_cc->EvalBootstrap(A_bar);
            Y     = m_cc->EvalBootstrap(Y);
        }
        Y     = eval_mult_lazy(Y, m_cc->EvalAdd(pI, A_bar));
        A_bar = eval_mult_lazy(A_bar, A_bar);
    }
    Y = eval_mult_lazy(Y, m_cc->EvalAdd(pI, A_bar));

    return Y;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_inverse_debug(const Ciphertext<DCRTPoly>& M,
                                                       const PrivateKey<DCRTPoly>& sk) {
    auto vI = initializeIdentityMatrix(d);
    Plaintext pI = m_cc->MakeCKKSPackedPlaintext(vI);

    std::cout << "\n=== eval_inverse_debug: Starting ===" << std::endl;

    // Decrypt input M
    Plaintext ptx_M;
    m_cc->Decrypt(sk, M, &ptx_M);
    ptx_M->SetLength(d * d);
    auto vec_M = ptx_M->GetRealPackedValue();
    std::cout << "Input M (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_M[i] << " ";
    }
    std::cout << "\nM level: " << M->GetLevel() << std::endl;

    // Transpose
    auto M_transposed = eval_transpose(M);
    Plaintext ptx_Mt;
    m_cc->Decrypt(sk, M_transposed, &ptx_Mt);
    ptx_Mt->SetLength(d * d);
    auto vec_Mt = ptx_Mt->GetRealPackedValue();
    std::cout << "\nM_transposed (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Mt[i] << " ";
    }
    std::cout << "\nM_transposed level: " << M_transposed->GetLevel() << std::endl;

    // MM^T
    auto MM_transposed = eval_mult(M, M_transposed);
    Plaintext ptx_MMt;
    m_cc->Decrypt(sk, MM_transposed, &ptx_MMt);
    ptx_MMt->SetLength(d * d);
    auto vec_MMt = ptx_MMt->GetRealPackedValue();
    std::cout << "\nMM_transposed (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_MMt[i] << " ";
    }
    std::cout << "\nMM_transposed level: " << MM_transposed->GetLevel() << std::endl;

    // Trace
    auto trace = eval_trace(MM_transposed, d * d);
    Plaintext ptx_trace;
    m_cc->Decrypt(sk, trace, &ptx_trace);
    ptx_trace->SetLength(d * d);
    auto vec_trace = ptx_trace->GetRealPackedValue();
    std::cout << "\ntrace (first element): " << vec_trace[0];
    std::cout << "\ntrace level: " << trace->GetLevel() << std::endl;

    // Trace reciprocal
    auto trace_reciprocal = m_cc->EvalDivide(trace, (d * d) / 3 - d, (d * d) / 3 + d, 5);
    Plaintext ptx_tr_recip;
    m_cc->Decrypt(sk, trace_reciprocal, &ptx_tr_recip);
    ptx_tr_recip->SetLength(d * d);
    auto vec_tr_recip = ptx_tr_recip->GetRealPackedValue();
    std::cout << "\ntrace_reciprocal (first element): " << vec_tr_recip[0];
    std::cout << "\ntrace_reciprocal level: " << trace_reciprocal->GetLevel() << std::endl;

    // Initial Y and A_bar
    auto Y = m_cc->EvalMultAndRelinearize(M_transposed, trace_reciprocal);
    Plaintext ptx_Y_init;
    m_cc->Decrypt(sk, Y, &ptx_Y_init);
    ptx_Y_init->SetLength(d * d);
    auto vec_Y_init = ptx_Y_init->GetRealPackedValue();
    std::cout << "\nInitial Y (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Y_init[i] << " ";
    }
    std::cout << "\nInitial Y level: " << Y->GetLevel() << std::endl;

    auto A_bar = m_cc->EvalSub(pI, m_cc->EvalMultAndRelinearize(MM_transposed, trace_reciprocal));
    Plaintext ptx_Abar_init;
    m_cc->Decrypt(sk, A_bar, &ptx_Abar_init);
    ptx_Abar_init->SetLength(d * d);
    auto vec_Abar_init = ptx_Abar_init->GetRealPackedValue();
    std::cout << "\nInitial A_bar (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Abar_init[i] << " ";
    }
    std::cout << "\nInitial A_bar level: " << A_bar->GetLevel() << std::endl;

    // Iterations
    for (int i = 0; i < r - 1; i++) {
        std::cout << "\n--- Iteration " << i << " ---" << std::endl;

        if (d >= 8 && static_cast<int>(Y->GetLevel()) >= depth - 2) {
            std::cout << "Bootstrapping (Y level=" << Y->GetLevel()
                      << ", A_bar level=" << A_bar->GetLevel() << ")" << std::endl;
            A_bar = m_cc->EvalBootstrap(A_bar);
            Y = m_cc->EvalBootstrap(Y);
            std::cout << "After bootstrap: Y level=" << Y->GetLevel()
                      << ", A_bar level=" << A_bar->GetLevel() << std::endl;
        }

        auto I_plus_Abar = m_cc->EvalAdd(pI, A_bar);
        Y = eval_mult(Y, I_plus_Abar);

        Plaintext ptx_Y_iter;
        m_cc->Decrypt(sk, Y, &ptx_Y_iter);
        ptx_Y_iter->SetLength(d * d);
        auto vec_Y_iter = ptx_Y_iter->GetRealPackedValue();
        std::cout << "Y after mult (first " << std::min(d*d, 8) << " elements): ";
        for (int j = 0; j < std::min(d*d, 8); j++) {
            std::cout << vec_Y_iter[j] << " ";
        }
        std::cout << "\nY level: " << Y->GetLevel() << std::endl;

        A_bar = eval_mult(A_bar, A_bar);

        Plaintext ptx_Abar_iter;
        m_cc->Decrypt(sk, A_bar, &ptx_Abar_iter);
        ptx_Abar_iter->SetLength(d * d);
        auto vec_Abar_iter = ptx_Abar_iter->GetRealPackedValue();
        std::cout << "A_bar after mult (first " << std::min(d*d, 8) << " elements): ";
        for (int j = 0; j < std::min(d*d, 8); j++) {
            std::cout << vec_Abar_iter[j] << " ";
        }
        std::cout << "\nA_bar level: " << A_bar->GetLevel() << std::endl;
    }

    // Final iteration
    std::cout << "\n--- Final iteration ---" << std::endl;
    Y = eval_mult(Y, m_cc->EvalAdd(pI, A_bar));

    Plaintext ptx_Y_final;
    m_cc->Decrypt(sk, Y, &ptx_Y_final);
    ptx_Y_final->SetLength(d * d);
    auto vec_Y_final = ptx_Y_final->GetRealPackedValue();
    std::cout << "Final Y (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Y_final[i] << " ";
    }
    std::cout << "\nFinal Y level: " << Y->GetLevel() << std::endl;
    std::cout << "=== eval_inverse_debug: Complete ===" << std::endl;

    return Y;
}

Ciphertext<DCRTPoly> MATINV_KCL25::eval_inverse_lazy_debug(const Ciphertext<DCRTPoly>& M,
                                                            const PrivateKey<DCRTPoly>& sk) {
    auto vI = initializeIdentityMatrix(d);
    Plaintext pI = m_cc->MakeCKKSPackedPlaintext(vI);

    std::cout << "\n=== eval_inverse_lazy_debug: Starting ===" << std::endl;

    // Decrypt input M
    Plaintext ptx_M;
    m_cc->Decrypt(sk, M, &ptx_M);
    ptx_M->SetLength(d * d);
    auto vec_M = ptx_M->GetRealPackedValue();
    std::cout << "Input M (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_M[i] << " ";
    }
    std::cout << "\nM level: " << M->GetLevel() << std::endl;

    // Transpose
    auto M_transposed = eval_transpose_lazy(M);
    Plaintext ptx_Mt;
    m_cc->Decrypt(sk, M_transposed, &ptx_Mt);
    ptx_Mt->SetLength(d * d);
    auto vec_Mt = ptx_Mt->GetRealPackedValue();
    std::cout << "\nM_transposed (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Mt[i] << " ";
    }
    std::cout << "\nM_transposed level: " << M_transposed->GetLevel() << std::endl;

    // MM^T (lazy)
    auto MM_transposed = eval_mult_lazy(M, M_transposed);
    Plaintext ptx_MMt;
    m_cc->Decrypt(sk, MM_transposed, &ptx_MMt);
    ptx_MMt->SetLength(d * d);
    auto vec_MMt = ptx_MMt->GetRealPackedValue();
    std::cout << "\nMM_transposed (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_MMt[i] << " ";
    }
    std::cout << "\nMM_transposed level: " << MM_transposed->GetLevel() << std::endl;

    // Trace
    auto trace = eval_trace_lazy(MM_transposed, d * d);
    Plaintext ptx_trace;
    m_cc->Decrypt(sk, trace, &ptx_trace);
    ptx_trace->SetLength(d * d);
    auto vec_trace = ptx_trace->GetRealPackedValue();
    std::cout << "\ntrace (first element): " << vec_trace[0];
    std::cout << "\ntrace level: " << trace->GetLevel() << std::endl;

    // Trace reciprocal
    auto trace_reciprocal = m_cc->EvalDivide(trace, (d * d) / 3 - d, (d * d) / 3 + d, 5);
    Plaintext ptx_tr_recip;
    m_cc->Decrypt(sk, trace_reciprocal, &ptx_tr_recip);
    ptx_tr_recip->SetLength(d * d);
    auto vec_tr_recip = ptx_tr_recip->GetRealPackedValue();
    std::cout << "\ntrace_reciprocal (first element): " << vec_tr_recip[0];
    std::cout << "\ntrace_reciprocal level: " << trace_reciprocal->GetLevel() << std::endl;

    // Initial Y and A_bar
    auto Y = m_cc->EvalMultAndRelinearize(M_transposed, trace_reciprocal);
    Plaintext ptx_Y_init;
    m_cc->Decrypt(sk, Y, &ptx_Y_init);
    ptx_Y_init->SetLength(d * d);
    auto vec_Y_init = ptx_Y_init->GetRealPackedValue();
    std::cout << "\nInitial Y (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Y_init[i] << " ";
    }
    std::cout << "\nInitial Y level: " << Y->GetLevel() << std::endl;

    auto A_bar = m_cc->EvalSub(pI, m_cc->EvalMultAndRelinearize(MM_transposed, trace_reciprocal));
    Plaintext ptx_Abar_init;
    m_cc->Decrypt(sk, A_bar, &ptx_Abar_init);
    ptx_Abar_init->SetLength(d * d);
    auto vec_Abar_init = ptx_Abar_init->GetRealPackedValue();
    std::cout << "\nInitial A_bar (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Abar_init[i] << " ";
    }
    std::cout << "\nInitial A_bar level: " << A_bar->GetLevel() << std::endl;

    // Iterations
    for (int i = 0; i < r - 1; i++) {
        std::cout << "\n--- Iteration " << i << " ---" << std::endl;

        if (d >= 8 && static_cast<int>(Y->GetLevel()) >= depth - 2) {
            std::cout << "Bootstrapping (Y level=" << Y->GetLevel()
                      << ", A_bar level=" << A_bar->GetLevel() << ")" << std::endl;
            A_bar = m_cc->EvalBootstrap(A_bar);
            Y = m_cc->EvalBootstrap(Y);
            std::cout << "After bootstrap: Y level=" << Y->GetLevel()
                      << ", A_bar level=" << A_bar->GetLevel() << std::endl;
        }

        auto I_plus_Abar = m_cc->EvalAdd(pI, A_bar);
        Y = eval_mult_lazy(Y, I_plus_Abar);

        Plaintext ptx_Y_iter;
        m_cc->Decrypt(sk, Y, &ptx_Y_iter);
        ptx_Y_iter->SetLength(d * d);
        auto vec_Y_iter = ptx_Y_iter->GetRealPackedValue();
        std::cout << "Y after mult (first " << std::min(d*d, 8) << " elements): ";
        for (int j = 0; j < std::min(d*d, 8); j++) {
            std::cout << vec_Y_iter[j] << " ";
        }
        std::cout << "\nY level: " << Y->GetLevel() << std::endl;

        A_bar = eval_mult_lazy(A_bar, A_bar);

        Plaintext ptx_Abar_iter;
        m_cc->Decrypt(sk, A_bar, &ptx_Abar_iter);
        ptx_Abar_iter->SetLength(d * d);
        auto vec_Abar_iter = ptx_Abar_iter->GetRealPackedValue();
        std::cout << "A_bar after mult (first " << std::min(d*d, 8) << " elements): ";
        for (int j = 0; j < std::min(d*d, 8); j++) {
            std::cout << vec_Abar_iter[j] << " ";
        }
        std::cout << "\nA_bar level: " << A_bar->GetLevel() << std::endl;
    }

    // Final iteration
    std::cout << "\n--- Final iteration ---" << std::endl;
    Y = eval_mult_lazy(Y, m_cc->EvalAdd(pI, A_bar));

    Plaintext ptx_Y_final;
    m_cc->Decrypt(sk, Y, &ptx_Y_final);
    ptx_Y_final->SetLength(d * d);
    auto vec_Y_final = ptx_Y_final->GetRealPackedValue();
    std::cout << "Final Y (first " << std::min(d*d, 16) << " elements): ";
    for (int i = 0; i < std::min(d*d, 16); i++) {
        std::cout << vec_Y_final[i] << " ";
    }
    std::cout << "\nFinal Y level: " << Y->GetLevel() << std::endl;
    std::cout << "=== eval_inverse_lazy_debug: Complete ===" << std::endl;

    return Y;
}

// ============================================================================
// Plan functions: collect required rotation indices
// ============================================================================

// Helper to create zero ciphertext for plan functions
Ciphertext<DCRTPoly> MATINV_KCL25::createZeroCT() const {
    // Create zero vector with d*d size (matching crypto context batch size)
    std::vector<double> zero(d * d, 0.0);
    auto pt = m_cc->MakeCKKSPackedPlaintext(zero);
    auto ct = m_cc->Encrypt(m_pk, pt);
    // Set slots to num_slots for batching (replicates the d*d data s times)
    ct->SetSlots(num_slots);
    return ct;
}

// Baseline plan functions - just collect rotation indices
void MATINV_KCL25::eval_transpose_plan(RotationKeyCollector& rk) const {
    for (int i = 1; i < d; i++) {
        rk.observe((d - 1) * i);
    }
    for (int i = -1; i > -d; i--) {
        rk.observe((d - 1) * i);
    }
}

void MATINV_KCL25::eval_trace_plan(RotationKeyCollector& rk, int batchSize) const {
    int logB = static_cast<int>(std::log2(static_cast<double>(batchSize)));
    for (int i = 1; i <= logB; i++) {
        rk.observe(batchSize / (1 << i));
    }
}

void MATINV_KCL25::vecRotsOpt_plan(RotationKeyCollector& rk, int is) const {
    for (int j = 0; j < s / np; j++) {
        rk.observe(is * d * s + j * d * np);
    }
}

void MATINV_KCL25::eval_mult_plan(RotationKeyCollector& rk) const {
    // Baby steps for A
    for (int i = 0; i < nb; i++) {
        rk.observe(i);
    }
    // Baby steps for B
    if (s >= np) {
        for (int i = 0; i < np; i++) {
            rk.observe(i * d);
        }
    }
    // vecRotsOpt
    for (int i = 0; i < B; i++) {
        vecRotsOpt_plan(rk, i);
    }
    // Giant steps
    for (int k = -ng; k < ng; k++) {
        if (k != 0) {
            rk.observe(k * nb);
        }
    }
    // Final aggregation
    int logS = static_cast<int>(std::log2(static_cast<double>(s)));
    for (int i = 1; i <= logS; i++) {
        rk.observe(num_slots / (1 << i));
    }
}

void MATINV_KCL25::eval_inverse_plan(RotationKeyCollector& rk) {
    eval_transpose_plan(rk);
    eval_mult_plan(rk);  // M * M^T
    eval_trace_plan(rk, d * d);
    // Loop: 2r-1 multiplications
    for (int i = 0; i < 2 * r - 1; i++) {
        eval_mult_plan(rk);
    }
}

// Lazy plan functions - execute on dummy ciphertexts to collect automorphism indices
void MATINV_KCL25::eval_transpose_lazy_plan(RotationKeyCollectorLazy& rk) const {
    auto M = createZeroCT();
    auto p0 = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(0));
    auto M_transposed = m_cc->EvalMult(M, p0);

    for (int i = 1; i < d; i++) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalLazyAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalLazyRotate(M, (d - 1) * i), pi));
    }
    for (int i = -1; i > -d; i--) {
        auto pi = m_cc->MakeCKKSPackedPlaintext(generateTransposeMsk(i));
        m_cc->EvalLazyAddInPlace(M_transposed, m_cc->EvalMult(m_cc->EvalLazyRotate(M, (d - 1) * i), pi));
    }
    // Collect accumulated automorphism indices before EvalBatchedKS
    rk.observeAutoIndices(M_transposed->GetElementKeyIndexVector());
}

void MATINV_KCL25::eval_trace_lazy_plan(RotationKeyCollectorLazy& rk, int batchSize) const {
    std::vector<int32_t> directRotations;
    // eval_trace_lazy uses EvalDirectRotate - these need immediate keys
    int logB = static_cast<int>(std::log2(static_cast<double>(batchSize)));
    for (int i = 1; i <= logB; i++) {
        directRotations.push_back(batchSize / (1 << i));
    }
    // Use observeRotationIndices to properly convert rotation indices to automorphism indices
    rk.observeRotationIndices(m_cc, directRotations);
}

void MATINV_KCL25::vecRotsOptLazy_plan(RotationKeyCollectorLazy& rk, const std::vector<Ciphertext<DCRTPoly>>& matrixM, int is) const {
    auto rotsM = zeroClone();
    for (int j = 0; j < s / np; j++) {
        auto T = zeroClone();
        for (int i = 0; i < np; i++) {
            auto msk = generateMaskVector(num_slots, np * j + i);
            msk = vectorRotate(msk, -is * d * s - j * d * np);
            auto pmsk = m_cc->MakeCKKSPackedPlaintext(msk, 1, 0, nullptr, num_slots);
            m_cc->EvalAddInPlace(T, m_cc->EvalMult(matrixM[i], pmsk));
        }
        m_cc->EvalLazyAddInPlace(rotsM, m_cc->EvalLazyRotate(T, is * d * s + j * d * np));
    }
    // Collect before EvalBatchedKS
    rk.observeAutoIndices(rotsM->GetElementKeyIndexVector());
}

void MATINV_KCL25::eval_mult_lazy_plan(RotationKeyCollectorLazy& rk) const {
    std::vector<int32_t> directRotations;

    auto matrixA = createZeroCT();
    auto matrixB = createZeroCT();
    auto matrixC = zeroClone();

    // Baby steps use EvalDirectRotate - collect these indices
    for (int i = 0; i < nb; i++) {
        directRotations.push_back(i);
    }
    if (s >= np) {
        for (int i = 0; i < np; i++) {
            directRotations.push_back(i * d);
        }
    }

    // Create dummy baby steps for vecRotsOptLazy
    std::vector<Ciphertext<DCRTPoly>> babyStepsOfB;
    if (s >= np) {
        for (int i = 0; i < np; i++) {
            babyStepsOfB.push_back(createZeroCT());
        }
    }

    // vecRotsOptLazy and giant steps use EvalLazyRotate
    for (int i = 0; i < B; i++) {
        if (s >= np) {
            vecRotsOptLazy_plan(rk, babyStepsOfB, i);
        }

        auto diagA = zeroClone();
        for (int k = -ng; k < ng; k++) {
            auto tmp = zeroClone();
            // ... (masking operations don't need rotation keys)
            m_cc->EvalLazyAddInPlace(diagA, m_cc->EvalLazyRotate(tmp, k * nb));
        }
        // Collect before EvalBatchedKS in eval_mult_lazy
        rk.observeAutoIndices(diagA->GetElementKeyIndexVector());
    }

    // Final aggregation uses EvalDirectRotate
    int logS = static_cast<int>(std::log2(static_cast<double>(s)));
    for (int i = 1; i <= logS; i++) {
        directRotations.push_back(num_slots / (1 << i));
    }

    // Use observeRotationIndices to properly convert rotation indices to automorphism indices
    rk.observeRotationIndices(m_cc, directRotations);
}

void MATINV_KCL25::eval_inverse_lazy_plan(RotationKeyCollectorLazy& rk) {
    eval_transpose_lazy_plan(rk);
    eval_mult_lazy_plan(rk);  // M * M^T
    eval_trace_lazy_plan(rk, d * d);
    // Loop: 2r-1 multiplications
    for (int i = 0; i < 2 * r - 1; i++) {
        eval_mult_lazy_plan(rk);
    }
}
