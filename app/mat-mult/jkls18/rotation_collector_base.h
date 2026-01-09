#pragma once
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <mutex>
#include <cstdint>

class RotationKeyCollector {
public:
    void begin(int slotCount, bool lazy) {
        std::lock_guard<std::mutex> g(mu_);
        slotCount_ = slotCount;
        lazy_ = lazy;
        indices_.clear();
        deps_.clear();
    }

    void observe(int k) {
        std::lock_guard<std::mutex> g(mu_);
        indices_.insert(normSigned(k));
    }

    int newStream() {
        std::lock_guard<std::mutex> g(mu_);
        deps_.emplace_back();
        deps_.back().insert(0);
        return static_cast<int>(deps_.size()) - 1;
    }

    template <class CryptoContextT, class PrivateKeyT>
    void generate(const CryptoContextT& cc, const PrivateKeyT& sk) {
        std::vector<int32_t> rots;
        {
            std::lock_guard<std::mutex> g(mu_);
            if (indices_.empty()) return;
            rots.reserve(indices_.size());
            for (int32_t k : indices_) rots.push_back(normSigned(k));
        }
        if (rots.empty()) return;
        std::sort(rots.begin(), rots.end());
        rots.erase(std::unique(rots.begin(), rots.end()), rots.end());
        if (lazy_) cc->EvalLazyRotateKeyGen(sk, rots);
        else       cc->EvalRotateKeyGen(sk, rots);
    }

    void clearCollected() {
        std::lock_guard<std::mutex> g(mu_);
        indices_.clear();
    }

private:
    int32_t normSigned(int v) const {
        if (slotCount_ <= 0) return 0;
        int r = v % slotCount_;
        if (r < 0) r += slotCount_;
        if (r > slotCount_ / 2) r -= slotCount_;
        return r;
    }

    mutable std::mutex mu_;
    int slotCount_{0};
    bool lazy_{false};
    std::unordered_set<int32_t> indices_;
    std::vector<std::unordered_set<int32_t>> deps_;
};
