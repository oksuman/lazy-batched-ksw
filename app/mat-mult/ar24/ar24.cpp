#include "ar24.h"
#include <algorithm>
#include <cmath>

MATMULT_AR24::MATMULT_AR24(const CryptoContext<DCRTPoly>& cc,
                           const PublicKey<DCRTPoly>& pk,
                           int dim)
    : m_cc(cc), m_PublicKey(pk), d(dim) {
    this->max_batch = static_cast<int>(m_cc->GetRingDimension()) / 2;
    this->s = std::min(max_batch / (d * d), d);
    if (s <= 0) s = 1;
    this->B = d / s;
}

MATMULT_AR24::MATMULT_AR24(int dim)
    : d(dim), max_batch(0), B(0), s(0) {}

// k: column index 0..d-1
std::vector<double> MATMULT_AR24::generatePhiMsk(int k) {
    std::vector<double> msk(d * d * s, 0.0);
    for (int i = k; i < d * d * s; i += d) {
        msk[i] = 1.0;
    }
    return msk;
}

// k: row index 0..d-1
std::vector<double> MATMULT_AR24::generatePsiMsk(int k) {
    std::vector<double> msk(d * d * s, 0.0);
    for (int i = 0; i < s; i++) {
        for (int j = i * d * d + k * d; j < i * d * d + k * d + d; j++) {
            msk[j] = 1.0;
        }
    }
    return msk;
}
// ---------------------- Baseline (real) ----------------------
Ciphertext<DCRTPoly> MATMULT_AR24::eval_mult(const Ciphertext<DCRTPoly>& matA,
                                             const Ciphertext<DCRTPoly>& matB) {
    const int num_slots = d * d * s;

    std::vector<double> zero(num_slots, 0.0);
    auto ptxC    = m_cc->MakeCKKSPackedPlaintext(zero);
    auto matrixC = m_cc->Encrypt(m_PublicKey, ptxC);

    auto matrixA = matA->Clone();
    auto matrixB = matB->Clone();

    std::vector<Ciphertext<DCRTPoly>> Tilde_A(B), Tilde_B(B);


    for (int i = 0; i < static_cast<int>(std::log2(s)); i++) {
        auto tmp = m_cc->EvalRotate(matrixA, (1 << i) - d * d * (1 << i));
        m_cc->EvalAddInPlace(matrixA, tmp);
    }

    for (int i = 0; i < static_cast<int>(std::log2(s)); i++) {
        auto tmp = m_cc->EvalRotate(matrixB, d * (1 << i) - d * d * (1 << i));
        m_cc->EvalAddInPlace(matrixB, tmp);
    }

    for (int i = 0; i < B; i++) {
        auto phi_si = m_cc->MakeCKKSPackedPlaintext(
            generatePhiMsk(s * i), 1, 0, nullptr, num_slots);

        auto tmp = m_cc->EvalMult(matrixA, phi_si);

    
        tmp = m_cc->EvalRotate(tmp, s * i);
        for (int j = 0; j < static_cast<int>(std::log2(d)); j++) {
            m_cc->EvalAddInPlace(tmp, m_cc->EvalRotate(tmp, -(1 << j)));
        }
        Tilde_A[i] = tmp;
    }

    for (int i = 0; i < B; i++) {
        auto psi_si = m_cc->MakeCKKSPackedPlaintext(
            generatePsiMsk(s * i), 1, 0, nullptr, num_slots);

        auto tmp = m_cc->EvalMult(matrixB, psi_si);
      
        tmp = m_cc->EvalRotate(tmp, s * i * d);
        for (int j = 0; j < static_cast<int>(std::log2(d)); j++) {
            m_cc->EvalAddInPlace(tmp, m_cc->EvalRotate(tmp, -(1 << j) * d));
        }
        Tilde_B[i] = tmp;
    }

    for (int i = 0; i < B; i++) {
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(Tilde_A[i], Tilde_B[i]));
    }

    // Step 7: final fold on matrixC
    for (int i = 0; i < static_cast<int>(std::log2(s)); i++) {
        m_cc->EvalAddInPlace(matrixC, m_cc->EvalRotate(matrixC, (d * d) * (1 << i)));
    }

    matrixC->SetSlots(d * d);
    return matrixC;
}

// ---------------------- Baseline (plan-only) ----------------------
void MATMULT_AR24::eval_mult_plan(RotationKeyCollector& rk) const {
    const int num_slots = d * d * s;                 // use batch size
    rk.begin(num_slots, /*lazy=*/false);             // FIX: was ringSlots

    const int Ls = static_cast<int>(std::log2(s));
    const int Ld = static_cast<int>(std::log2(d));

    // A preprocessing
    for (int i = 0; i < Ls; ++i)
        rk.observe((1 << i) - d * d * (1 << i));

    // B preprocessing
    for (int i = 0; i < Ls; ++i)
        rk.observe(d * (1 << i) - d * d * (1 << i));

    // Tilde_A[i]
    for (int i = 0; i < B; ++i) {
        rk.observe(s * i);
        for (int j = 0; j < Ld; ++j)
            rk.observe(-(1 << j));
    }

    // Tilde_B[i]
    for (int i = 0; i < B; ++i) {
        rk.observe(s * i * d);
        for (int j = 0; j < Ld; ++j)
            rk.observe(-(1 << j) * d);
    }

    // Final folding into C
    for (int i = 0; i < Ls; ++i)
        rk.observe((d * d) * (1 << i));
}



// ---------------------- Lazy (real) ----------------------
Ciphertext<DCRTPoly> MATMULT_AR24::eval_mult_lazy(const Ciphertext<DCRTPoly>& matA,
                                                  const Ciphertext<DCRTPoly>& matB) {
    const int num_slots = d * d * s;

    std::vector<double> zero(num_slots, 0.0);
    auto ptxC    = m_cc->MakeCKKSPackedPlaintext(zero);
    auto matrixC = m_cc->Encrypt(m_PublicKey, ptxC);

    auto matrixA = matA->Clone();
    auto matrixB = matB->Clone();

    std::vector<Ciphertext<DCRTPoly>> Tilde_A(B), Tilde_B(B);

    {
        int rot_count = 0;
        const int rounds = static_cast<int>(std::log2(s));
        for (int i = 0; i < rounds; i++) {
            auto tmp = m_cc->EvalLazyRotate(matrixA, (1 << i) - d * d * (1 << i));
            m_cc->EvalLazyAddInPlace(matrixA, tmp);
            rot_count++;

            const bool last = (i == rounds - 1);
            if (rot_count == 3 || last) {
                matrixA = m_cc->EvalBatchedKS(matrixA);
                rot_count = 0;
            }
        }
    }

    {
        int rot_count = 0;
        const int rounds = static_cast<int>(std::log2(s));
        for (int i = 0; i < rounds; i++) {
            auto tmp = m_cc->EvalLazyRotate(matrixB, d * (1 << i) - d * d * (1 << i));
            m_cc->EvalLazyAddInPlace(matrixB, tmp);
            rot_count++;

            const bool last = (i == rounds - 1);
            if (rot_count == 3 || last) {
                matrixB = m_cc->EvalBatchedKS(matrixB);
                rot_count = 0;
            }
        }
    }

    for (int i = 0; i < B; i++) {
        auto phi_si = m_cc->MakeCKKSPackedPlaintext(
            generatePhiMsk(s * i), 1, 0, nullptr, num_slots);

        auto tmp = m_cc->EvalMult(matrixA, phi_si);

        {
            int rot_count = 0;

            // initial rotate by s*i counts as one rotation
            tmp = m_cc->EvalLazyRotate(tmp, s * i);
            rot_count++;

            const int rounds = static_cast<int>(std::log2(d));
            for (int j = 0; j < rounds; j++) {
                // bundle (rotate + add) counts as one rotation
                auto r = m_cc->EvalLazyRotate(tmp, -(1 << j));
                m_cc->EvalLazyAddInPlace(tmp, r);
                rot_count++;

                const bool last = (j == rounds - 1);
                if (rot_count == 3 || last) {
                    tmp = m_cc->EvalBatchedKS(tmp);
                    rot_count = 0;
                }
            }
            // ensure canonical form
            tmp = m_cc->EvalBatchedKS(tmp);
        }
        Tilde_A[i] = tmp;
    }

    for (int i = 0; i < B; i++) {
        auto psi_si = m_cc->MakeCKKSPackedPlaintext(
            generatePsiMsk(s * i), 1, 0, nullptr, num_slots);

        auto tmp = m_cc->EvalMult(matrixB, psi_si);
        {
            int rot_count = 0;

            // initial rotate by s*i*d counts as one rotation
            tmp = m_cc->EvalLazyRotate(tmp, s * i * d);
            rot_count++;

            const int rounds = static_cast<int>(std::log2(d));
            for (int j = 0; j < rounds; j++) {
                auto r = m_cc->EvalLazyRotate(tmp, -(1 << j) * d);
                m_cc->EvalLazyAddInPlace(tmp, r);
                rot_count++;

                const bool last = (j == rounds - 1);
                if (rot_count == 3 || last) {
                    tmp = m_cc->EvalBatchedKS(tmp);
                    rot_count = 0;
                }
            }
            // ensure canonical form
            tmp = m_cc->EvalBatchedKS(tmp);
        }
        Tilde_B[i] = tmp;
    }

    for (int i = 0; i < B; i++) {
        m_cc->EvalAddInPlace(
            matrixC, m_cc->EvalMultAndRelinearize(Tilde_A[i], Tilde_B[i]));
    }

    {
        int rot_count = 0;
        const int rounds = static_cast<int>(std::log2(s));
        for (int i = 0; i < rounds; i++) {
            matrixC = m_cc->EvalLazyAdd(matrixC, m_cc->EvalLazyRotate(matrixC, (d * d) * (1 << i)));
            rot_count++;

            const bool last = (i == rounds - 1);
            if (rot_count == 3 || last) {
                matrixC = m_cc->EvalBatchedKS(matrixC);
                rot_count = 0;
            }
        }
        // ensure canonical form anyway
        matrixC = m_cc->EvalBatchedKS(matrixC);
    }

    matrixC->SetSlots(d * d);
    return matrixC;
}

// ---------------------- Lazy (plan-only) ----------------------
void MATMULT_AR24::eval_mult_lazy_plan(RotationKeyCollectorLazy& rk) const {
    const int num_slots = d * d * s;
    rk.begin(num_slots);

    const int Ls = static_cast<int>(std::log2(s));
    const int Ld = static_cast<int>(std::log2(d));
    const int K  = 3; // observe every 3 lazy rotations (or at the end of the loop)

    auto makeZeroCT = [&]() {
        std::vector<double> zero(num_slots, 0.0);
        auto pt = m_cc->MakeCKKSPackedPlaintext(zero);
        return m_cc->Encrypt(m_PublicKey, pt);
    };

    // Step 1: same algorithm; observe where EvalBatchedKS would be called.
    {
        auto ct = makeZeroCT();
        int rot_count = 0;
        for (int i = 0; i < Ls; ++i) {
            ct = m_cc->EvalLazyAdd(ct, m_cc->EvalLazyRotate(ct, (1 << i) - d * d * (1 << i)));
            ++rot_count;
            const bool last = (i == Ls - 1);
            if (rot_count == K || last) {
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());
                rot_count = 0;
            }
        }
    }

    // Step 2
    {
        auto ct = makeZeroCT();
        int rot_count = 0;
        for (int i = 0; i < Ls; ++i) {
            ct = m_cc->EvalLazyAdd(ct, m_cc->EvalLazyRotate(ct, d * (1 << i) - d * d * (1 << i)));
            ++rot_count;
            const bool last = (i == Ls - 1);
            if (rot_count == K || last) {
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());
                rot_count = 0;
            }
        }
    }

    // Step 4 path: initial rotate (s*i) counts as one; each (rotate + add) bundle counts as one.
    for (int i = 0; i < B; ++i) {
        auto ct = makeZeroCT();
        int rot_count = 0;

        // initial rotate: s*i
        ct = m_cc->EvalLazyRotate(ct, s * i);
        ++rot_count;
        if (rot_count == K) {
            rk.observeAutoIndices(ct->GetElementKeyIndexVector());
            rot_count = 0;
        }

        for (int j = 0; j < Ld; ++j) {
            ct = m_cc->EvalLazyAdd(ct, m_cc->EvalLazyRotate(ct, -(1 << j)));
            ++rot_count;
            const bool last = (j == Ld - 1);
            if (rot_count == K || last) {
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());
                rot_count = 0;
            }
        }
    }

    // Step 6 path
    for (int i = 0; i < B; ++i) {
        auto ct = makeZeroCT();
        int rot_count = 0;

        // initial rotate: s*i*d
        ct = m_cc->EvalLazyRotate(ct, s * i * d);
        ++rot_count;
        if (rot_count == K) {
            rk.observeAutoIndices(ct->GetElementKeyIndexVector());
            rot_count = 0;
        }

        for (int j = 0; j < Ld; ++j) {
            ct = m_cc->EvalLazyAdd(ct, m_cc->EvalLazyRotate(ct, -(1 << j) * d));
            ++rot_count;
            const bool last = (j == Ld - 1);
            if (rot_count == K || last) {
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());
                rot_count = 0;
            }
        }
    }

    // Step 7
    {
        auto ct = makeZeroCT();
        int rot_count = 0;
        for (int i = 0; i < Ls; ++i) {
            ct = m_cc->EvalLazyAdd(ct, m_cc->EvalLazyRotate(ct, (d * d) * (1 << i)));
            ++rot_count;
            const bool last = (i == Ls - 1);
            if (rot_count == K || last) {
                rk.observeAutoIndices(ct->GetElementKeyIndexVector());
                rot_count = 0;
            }
        }
    }
}
