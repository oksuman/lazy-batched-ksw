#include "cbc_conv.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

// ---------- Constructors ----------

CBC_Conv::CBC_Conv(const CryptoContext<DCRTPoly>& cc,
                   const PublicKey<DCRTPoly>& pk,
                   int H, int W, int c_in, int c_out, int gap)
    : m_cc(cc), m_pk(pk), m_H(H), m_W(W), m_cin(c_in), m_cout(c_out), m_gap(gap)
{
    // Compute 9 rotation offsets: -(j-1)*W*gap - (k-1)*gap  for j,k in {0,1,2}
    m_offsets.reserve(9);
    for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
            m_offsets.push_back(-(j - 1) * m_W * m_gap - (k - 1) * m_gap);
}

CBC_Conv::CBC_Conv(int H, int W, int c_in, int c_out, int gap)
    : m_H(H), m_W(W), m_cin(c_in), m_cout(c_out), m_gap(gap)
{
    m_offsets.reserve(9);
    for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
            m_offsets.push_back(-(j - 1) * m_W * m_gap - (k - 1) * m_gap);
}

// ---------- Kernel encoding ----------

void CBC_Conv::encodeKernels(
    const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
    const std::vector<double>& bias)
{
    int numSlots = m_H * m_W;

    m_kernels.resize(m_cout);
    for (int i = 0; i < m_cout; i++) {
        m_kernels[i].resize(m_cin);
        for (int ch = 0; ch < m_cin; ch++) {
            m_kernels[i][ch].resize(9);
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    double w = weights[i][ch][j][k];
                    std::vector<double> v(numSlots, w);
                    m_kernels[i][ch][j * 3 + k] = m_cc->MakeCKKSPackedPlaintext(v);
                }
            }
        }
    }

    m_bias.resize(m_cout);
    for (int i = 0; i < m_cout; i++) {
        std::vector<double> v(numSlots, bias[i]);
        m_bias[i] = m_cc->MakeCKKSPackedPlaintext(v);
    }
}

// ---------- CBC ConvBN ----------

std::vector<Ciphertext<DCRTPoly>> CBC_Conv::eval_baseline(
    const std::vector<Ciphertext<DCRTPoly>>& input)
{
    auto t_rot_start = std::chrono::steady_clock::now();

    std::vector<std::vector<Ciphertext<DCRTPoly>>> rotated(
        m_cin, std::vector<Ciphertext<DCRTPoly>>(9));
    for (int ch = 0; ch < m_cin; ch++) {
        for (int f = 0; f < 9; f++) {
            if (m_offsets[f] == 0)
                rotated[ch][f] = input[ch];
            else
                rotated[ch][f] = m_cc->EvalRotate(input[ch], m_offsets[f]);
        }
    }

    auto t_rot_end = std::chrono::steady_clock::now();
    double ms_rot = std::chrono::duration<double, std::milli>(t_rot_end - t_rot_start).count();

    double ms_pmult = 0.0, ms_add = 0.0;
    std::vector<Ciphertext<DCRTPoly>> results(m_cout);
    for (int out = 0; out < m_cout; out++) {
        Ciphertext<DCRTPoly> result;
        bool first = true;
        for (int ch = 0; ch < m_cin; ch++) {
            for (int f = 0; f < 9; f++) {
                auto tp0 = std::chrono::steady_clock::now();
                auto term = m_cc->EvalMult(rotated[ch][f],
                                           m_kernels[out][ch][f]);
                auto tp1 = std::chrono::steady_clock::now();
                ms_pmult += std::chrono::duration<double, std::milli>(tp1 - tp0).count();

                if (first) { result = term; first = false; }
                else {
                    auto ta0 = std::chrono::steady_clock::now();
                    m_cc->EvalAddInPlace(result, term);
                    auto ta1 = std::chrono::steady_clock::now();
                    ms_add += std::chrono::duration<double, std::milli>(ta1 - ta0).count();
                }
            }
        }
        results[out] = result;
    }

    std::cout << "      [baseline breakdown] rot=" << std::fixed << std::setprecision(1)
              << ms_rot << " ms, pmult=" << ms_pmult << " ms, add=" << ms_add << " ms"
              << ", total=" << (ms_rot + ms_pmult + ms_add) << " ms\n";

    return results;
}

std::vector<Ciphertext<DCRTPoly>> CBC_Conv::eval_hoisted(
    const std::vector<Ciphertext<DCRTPoly>>& input)
{
    auto n = m_cc->GetRingDimension();
    auto M = 2 * n;

    // Phase 1: compute all rotations once (hoisted EvalFastRotation)
    std::vector<std::vector<Ciphertext<DCRTPoly>>> rotated(
        m_cin, std::vector<Ciphertext<DCRTPoly>>(9));
    for (int ch = 0; ch < m_cin; ch++) {
        auto digits = m_cc->EvalFastRotationPrecompute(input[ch]);
        for (int f = 0; f < 9; f++) {
            if (m_offsets[f] == 0)
                rotated[ch][f] = input[ch];
            else
                rotated[ch][f] = m_cc->EvalFastRotation(
                    input[ch], m_offsets[f], M, digits);
        }
    }

    std::vector<Ciphertext<DCRTPoly>> results(m_cout);
    for (int out = 0; out < m_cout; out++) {
        Ciphertext<DCRTPoly> result;
        bool first = true;
        for (int ch = 0; ch < m_cin; ch++) {
            for (int f = 0; f < 9; f++) {
                auto term = m_cc->EvalMult(rotated[ch][f],
                                           m_kernels[out][ch][f]);
                if (first) { result = term; first = false; }
                else       m_cc->EvalAddInPlace(result, term);
            }
        }
        results[out] = result;
    }
    return results;
}

std::vector<Ciphertext<DCRTPoly>> CBC_Conv::eval_lazy(
    const std::vector<Ciphertext<DCRTPoly>>& input)
{
    // Phase 1: lazy rotate all inputs once (no key-switching)
    std::vector<std::vector<Ciphertext<DCRTPoly>>> rotated(
        m_cin, std::vector<Ciphertext<DCRTPoly>>(9));
    for (int ch = 0; ch < m_cin; ch++) {
        for (int f = 0; f < 9; f++) {
            if (m_offsets[f] == 0)
                rotated[ch][f] = input[ch];
            else
                rotated[ch][f] = m_cc->EvalLazyRotate(input[ch], m_offsets[f]);
        }
    }

    // Phase 2: per output channel, accumulate lazily then batched KS
    std::vector<Ciphertext<DCRTPoly>> results(m_cout);
    for (int out = 0; out < m_cout; out++) {
        Ciphertext<DCRTPoly> acc;
        bool first = true;
        for (int ch = 0; ch < m_cin; ch++) {
            for (int f = 0; f < 9; f++) {
                auto term = m_cc->EvalMult(rotated[ch][f],
                                           m_kernels[out][ch][f]);
                if (first) { acc = term; first = false; }
                else       m_cc->EvalLazyAddInPlace(acc, term);
            }
        }
        results[out] = m_cc->EvalBatchedKS(acc);
    }
    return results;
}

std::vector<Ciphertext<DCRTPoly>> CBC_Conv::eval_twostage(
    const std::vector<Ciphertext<DCRTPoly>>& input)
{
    // Two-stage optimization:
    // Stage 1: Row shifts ({-W*gap, 0, +W*gap}) via hoisted EvalDirectRotate (shared across output channels)
    // Stage 2: Col shifts ({-gap, 0, +gap}) via lazy EvalLazyRotate + EvalBatchedKS(m=2)

    int row_offset = m_W * m_gap;  // row shift amount

    // Stage 1: Compute row-shifted ciphertexts for each input channel
    // row_cts[ch][0] = rotate by +W*gap (j=0), row_cts[ch][1] = identity (j=1), row_cts[ch][2] = rotate by -W*gap (j=2)
    std::vector<std::vector<Ciphertext<DCRTPoly>>> row_cts(
        m_cin, std::vector<Ciphertext<DCRTPoly>>(3));
    for (int ch = 0; ch < m_cin; ch++) {
        auto precomp = m_cc->EvalDirectRotatePrecompute(input[ch]);
        row_cts[ch][0] = m_cc->EvalDirectRotate(input[ch], row_offset, precomp);   // j=0: +W*gap
        row_cts[ch][1] = input[ch];                                                 // j=1: identity
        row_cts[ch][2] = m_cc->EvalDirectRotate(input[ch], -row_offset, precomp);   // j=2: -W*gap
    }

    // Pre-compute col shifts for all row-shifted ciphertexts (shared across output channels)
    // col_shifted[ch][j*3 + k] corresponds to row j, col k
    // Col offsets: k=0 → +gap, k=1 → 0 (identity), k=2 → -gap
    std::vector<std::vector<Ciphertext<DCRTPoly>>> col_shifted(
        m_cin, std::vector<Ciphertext<DCRTPoly>>(9));
    for (int ch = 0; ch < m_cin; ch++) {
        for (int j = 0; j < 3; j++) {
            col_shifted[ch][j * 3 + 0] = m_cc->EvalLazyRotate(row_cts[ch][j], m_gap);    // k=0: +gap
            col_shifted[ch][j * 3 + 1] = row_cts[ch][j];                                  // k=1: identity
            col_shifted[ch][j * 3 + 2] = m_cc->EvalLazyRotate(row_cts[ch][j], -m_gap);   // k=2: -gap
        }
    }

    // Stage 2: Per output channel, ptMult + lazy accumulate + BatchKS(m=2)
    std::vector<Ciphertext<DCRTPoly>> results(m_cout);
    for (int out = 0; out < m_cout; out++) {
        Ciphertext<DCRTPoly> acc;
        bool first = true;
        for (int ch = 0; ch < m_cin; ch++) {
            for (int f = 0; f < 9; f++) {
                auto term = m_cc->EvalMult(col_shifted[ch][f], m_kernels[out][ch][f]);
                if (first) { acc = term; first = false; }
                else       m_cc->EvalLazyAddInPlace(acc, term);
            }
        }
        results[out] = m_cc->EvalBatchedKS(acc);  // BatchKS with m=2 (only {-gap, +gap} keys)
    }
    return results;
}

// ---------- Plaintext convolution ----------

std::vector<std::vector<double>> CBC_Conv::eval_plain(
    const std::vector<std::vector<double>>& input_channels,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
    const std::vector<double>& bias) const
{
    int numSlots = m_H * m_W;
    std::vector<std::vector<double>> output(m_cout, std::vector<double>(numSlots, 0.0));

    for (int out = 0; out < m_cout; out++) {
        for (int ch = 0; ch < m_cin; ch++) {
            for (int i = 0; i < numSlots; i++) {
                for (int f = 0; f < 9; f++) {
                    int src = ((i + m_offsets[f]) % numSlots + numSlots) % numSlots;
                    output[out][i] += input_channels[ch][src]
                                      * weights[out][ch][f / 3][f % 3];
                }
            }
        }
        for (int i = 0; i < numSlots; i++)
            output[out][i] += bias[out];
    }
    return output;
}

// ---------- Proper 2D convolution (ground truth) ----------

std::vector<std::vector<double>> CBC_Conv::eval_plain_2d(
    const std::vector<std::vector<double>>& input_channels,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
    const std::vector<double>& bias) const
{
    std::vector<std::vector<double>> output(m_cout,
        std::vector<double>(m_H * m_W, 0.0));

    for (int out = 0; out < m_cout; out++) {
        for (int ch = 0; ch < m_cin; ch++) {
            for (int r = 0; r < m_H; r++) {
                for (int c = 0; c < m_W; c++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int dk = -1; dk <= 1; dk++) {
                            int nr = r + dj, nc = c + dk;
                            if (nr >= 0 && nr < m_H && nc >= 0 && nc < m_W) {
                                // kernel index: j = 1-dj, k = 1-dk
                                // (matches offset = -(j-1)*W - (k-1) = dj*W + dk)
                                int j = 1 - dj, k = 1 - dk;
                                output[out][r * m_W + c] +=
                                    input_channels[ch][nr * m_W + nc]
                                    * weights[out][ch][j][k];
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < m_H * m_W; i++)
            output[out][i] += bias[out];
    }
    return output;
}

// ---------- Plan functions ----------

void CBC_Conv::eval_baseline_plan(RotationKeyCollector& rk) const {
    for (int f = 0; f < 9; f++)
        if (m_offsets[f] != 0)
            rk.observe(m_offsets[f]);
}

void CBC_Conv::eval_hoisted_plan(RotationKeyCollector& rk) const {
    for (int f = 0; f < 9; f++)
        if (m_offsets[f] != 0)
            rk.observe(m_offsets[f]);
}

void CBC_Conv::eval_lazy_plan(RotationKeyCollectorLazy& rk) const {
    int numSlots = m_H * m_W;
    rk.begin(numSlots);

    std::vector<double> zero(numSlots, 0.0);
    auto pt = m_cc->MakeCKKSPackedPlaintext(zero);
    auto dummy = m_cc->Encrypt(m_pk, pt);

    // Simulate lazy rotations for 1 input channel (all channels share same offsets)
    Ciphertext<DCRTPoly> acc = dummy;
    for (int f = 0; f < 9; f++) {
        if (m_offsets[f] != 0) {
            auto rotated = m_cc->EvalLazyRotate(dummy, m_offsets[f]);
            m_cc->EvalLazyAddInPlace(acc, rotated);
        }
    }
    rk.observeAutoIndices(acc->GetElementKeyIndexVector());
}

void CBC_Conv::eval_twostage_plan(RotationKeyCollector& rk) const {
    int row_offset = m_W * m_gap;

    // Stage 1: Row shifts use EvalDirectRotate (needs lazy keys for ±row_offset)
    rk.observe(row_offset);
    rk.observe(-row_offset);

    // Stage 2: Col shifts use EvalLazyRotate (needs lazy keys for ±gap)
    rk.observe(m_gap);
    rk.observe(-m_gap);
}
