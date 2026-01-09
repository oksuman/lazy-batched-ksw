#pragma once
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <mutex>
#include <cstdint>

class RotationKeyCollectorLazy {
public:
    void begin(int slotCount) {
        std::lock_guard<std::mutex> g(mu_);
        slotCount_ = slotCount;
        autoIdx_.clear();
    }

    void observeAutoIndices(const std::vector<int32_t>& v) {
        std::lock_guard<std::mutex> g(mu_);
        for (int32_t a : v) autoIdx_.insert(a);
    }

    template <class CryptoContextT, class PrivateKeyT>
    void generate(const CryptoContextT& cc, const PrivateKeyT& sk) {
        std::vector<int32_t> exps;
        {
            std::lock_guard<std::mutex> g(mu_);
            const int ringDim = static_cast<int>(cc->GetRingDimension());
            const int M = 2 * ringDim;
            if (ringDim <= 0) return;
            std::vector<int32_t> pow5(ringDim), log5(M, -1);
            pow5[0] = 1 % M;
            log5[pow5[0]] = 0;
            for (int t = 1; t < ringDim; ++t) {
                pow5[t] = static_cast<int32_t>((static_cast<int64_t>(pow5[t - 1]) * 5) % M);
                if (log5[pow5[t]] < 0) log5[pow5[t]] = t;
            }
            exps.reserve(autoIdx_.size() + 1);
            exps.push_back(0);
            for (int32_t a : autoIdx_) {
                if (a <= 0 || a >= M) continue;
                int32_t t = log5[a];
                if (t >= 0) exps.push_back(t);
            }
        }
        if (exps.empty()) return;
        std::sort(exps.begin(), exps.end());
        exps.erase(std::unique(exps.begin(), exps.end()), exps.end());
        cc->EvalLazyRotateKeyGen(sk, exps);
    }

    void clearCollected() {
        std::lock_guard<std::mutex> g(mu_);
        autoIdx_.clear();
    }

    std::vector<int32_t> getCollectedAutoIndices() const {
        std::lock_guard<std::mutex> g(mu_);
        std::vector<int32_t> out(autoIdx_.begin(), autoIdx_.end());
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

private:
    mutable std::mutex mu_;
    int slotCount_{0};
    std::unordered_set<int32_t> autoIdx_;
};
