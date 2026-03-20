#pragma once

#include <openfhe.h>
#include <memory>
#include <vector>
#include "rotation_collector_base.h"
#include "rotation_collector_lazy.h"

using namespace lbcrypto;

class CBC_Conv {
public:
    CBC_Conv(const CryptoContext<DCRTPoly>& cc,
             const PublicKey<DCRTPoly>& pk,
             int H, int W, int c_in, int c_out,
             int gap = 1);
    CBC_Conv(int H, int W, int c_in, int c_out, int gap = 1);

    void encodeKernels(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
        const std::vector<double>& bias);

    // CBC ConvBN (Algorithm 1) — 4 variants
    std::vector<Ciphertext<DCRTPoly>> eval_baseline(
        const std::vector<Ciphertext<DCRTPoly>>& input);
    std::vector<Ciphertext<DCRTPoly>> eval_hoisted(
        const std::vector<Ciphertext<DCRTPoly>>& input);
    std::vector<Ciphertext<DCRTPoly>> eval_lazy(
        const std::vector<Ciphertext<DCRTPoly>>& input);
    std::vector<Ciphertext<DCRTPoly>> eval_twostage(
        const std::vector<Ciphertext<DCRTPoly>>& input);

    // Plaintext convolution (circular, matches HE behavior)
    std::vector<std::vector<double>> eval_plain(
        const std::vector<std::vector<double>>& input_channels,
        const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
        const std::vector<double>& bias) const;

    // Proper 2D convolution with zero-padding (ground truth)
    std::vector<std::vector<double>> eval_plain_2d(
        const std::vector<std::vector<double>>& input_channels,
        const std::vector<std::vector<std::vector<std::vector<double>>>>& weights,
        const std::vector<double>& bias) const;

    // Plan functions to collect rotation indices
    void eval_baseline_plan(RotationKeyCollector& rk) const;
    void eval_hoisted_plan(RotationKeyCollector& rk) const;
    void eval_lazy_plan(RotationKeyCollectorLazy& rk) const;
    void eval_twostage_plan(RotationKeyCollector& rk) const;

    const std::vector<int>& offsets() const { return m_offsets; }

private:
    CryptoContext<DCRTPoly> m_cc;
    PublicKey<DCRTPoly>     m_pk;
    int m_H, m_W, m_cin, m_cout, m_gap;
    std::vector<int> m_offsets;  // 9 rotation offsets
    std::vector<std::vector<std::vector<Plaintext>>> m_kernels;  // [c_out][c_in][9]
    std::vector<Plaintext> m_bias;
};
