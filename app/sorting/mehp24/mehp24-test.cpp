#include "sorting.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include "openfhe.h"
#include "rotation_collector_lazy.h"
#include "bench_stages.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

using namespace lbcrypto;
using benchstages::StageTotals;
using benchstages::stages_reset;
using benchstages::stages_snapshot;

// Enable both paths
#define RUN_FG       1
#define RUN_FG_LAZY  1

// ================= User-tunable =================
static const int      kTrials          = 10;
static const unsigned kSeed            = 1337u;   // keep for reproducibility if RNG used elsewhere
static const int      kNs[]            = {8, 16, 32, 64};
// static const int      kNs[]            = {16, 32, 64, 128, 256};
static const char*    kCSV             = "sort_bench_results.csv";
static const usint    kRingDim         = (1u << 17);
static const int      kTrialCooldownMs = 150;
static const int      kBetweenNMs      = 1500;
// ================================================

#if defined(__GLIBC__)
  #include <malloc.h>
  static inline void TrimMalloc() { malloc_trim(0); }
#else
  static inline void TrimMalloc() {}
#endif

static inline void CooldownMs(int ms) {
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

static inline void HardReset() {
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    TrimMalloc();
}

struct Stats { double mean_s{0.0}; double std_s{0.0}; };
struct Acc   { double max_abs_err{0.0}; double mse{0.0}; };
struct Row   {
    size_t n;
    Stats st;
    Acc ac;
    std::array<double,7> stages_mean_s{};
};

static Stats meanStd(const std::vector<double>& v){
    if (v.empty()) return {};
    double m = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
    double var = 0.0;
    for (double x : v) var += (x - m) * (x - m);
    var /= double(v.size());
    return {m, std::sqrt(var)};
}

static Acc accuracy(const std::vector<double>& a, const std::vector<double>& b){
    size_t n = std::min(a.size(), b.size());
    double mx = 0.0;
    long double sse = 0.0L;
    for (size_t i = 0; i < n; ++i) {
        double e = std::abs(a[i] - b[i]);
        if (e > mx) mx = e;
        sse += (long double)e * e;
    }
    return {mx, (double)(sse / (long double)std::max<size_t>(1, n))};
}

static std::vector<double> loadUniqueCSV(size_t n, const std::string& path = "data/points1d.csv"){
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("open: " + path);
    std::string text((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
    fin.close();
    for (char& c : text) if (c == ',') c = ' ';
    std::istringstream iss(text);
    std::unordered_set<std::string> seen; seen.reserve(n * 2);
    std::vector<double> out; out.reserve(n);
    std::string tok;
    while (iss >> tok) {
        if (!seen.insert(tok).second) continue;
        try {
            double v = std::stod(tok);
            out.push_back(v);
            if (out.size() == n) break;
        } catch (...) {}
    }
    if (out.size() < n) throw std::runtime_error("not enough unique values");
    return out;
}

struct Ctx {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly>     pk;
    PrivateKey<DCRTPoly>    sk;
    size_t slots{0};
};

static void paramsFG(size_t n, size_t& dg_c, size_t& df_c, size_t& dg_i, size_t& df_i, usint& depth, usint& slots){
    dg_c = 3; df_c = 2; dg_i = (size_t)((std::log2((double)n) + 1) / 2); df_i = 2;
    depth = 4 * (dg_c + df_c + dg_i + df_i) + 4;
    slots = (usint)(n * n);
}

static Ctx makeCtxFG(size_t n){
    size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
    paramsFG(n, dg_c, df_c, dg_i, df_i, depth, slots);

    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(depth);
    P.SetBatchSize(slots);
    P.SetScalingModSize(40);
    P.SetRingDim(kRingDim);
    P.SetKeySwitchTechnique(HYBRID);
    P.SetNumLargeDigits(5);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    // Baseline: pre-generate rotation keys
    auto indices = getRotationIndices(n);
    cc->EvalRotateKeyGen(kp.secretKey, indices);

    return {cc, kp.publicKey, kp.secretKey, slots};
}

static Ctx makeCtxFGLazy(size_t n){
    size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
    paramsFG(n, dg_c, df_c, dg_i, df_i, depth, slots);

    CCParams<CryptoContextCKKSRNS> P;
    P.SetMultiplicativeDepth(depth);
    P.SetBatchSize(slots);
    P.SetScalingModSize(40);
    P.SetRingDim(kRingDim);
    P.SetKeySwitchTechnique(BATCHED);
    P.SetNumLargeDigits(5);

    auto cc = GenCryptoContext(P);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    RotationKeyCollectorLazy rc; rc.begin((int)slots);
    Plan_sortFGLazy(cc, kp.publicKey, n, dg_c, df_c, dg_i, df_i, rc);
    rc.generate(cc, kp.secretKey);

    return {cc, kp.publicKey, kp.secretKey, slots};
}

static Ciphertext<DCRTPoly> run_sortFG(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk,
                                       const std::vector<double>& v, size_t n){
    size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
    paramsFG(n, dg_c, df_c, dg_i, df_i, depth, slots);
    auto ct = cc->Encrypt(pk, cc->MakeCKKSPackedPlaintext(v));
    return sortFG(ct, n, dg_c, df_c, dg_i, df_i);
}

static Ciphertext<DCRTPoly> run_sortFGLazy(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk,
                                           const std::vector<double>& v, size_t n){
    size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
    paramsFG(n, dg_c, df_c, dg_i, df_i, depth, slots);
    auto ct = cc->Encrypt(pk, cc->MakeCKKSPackedPlaintext(v));
    return sortFGLazy(ct, n, dg_c, df_c, dg_i, df_i);
}

static void printParams(const std::string& tag,
                        const CryptoContext<DCRTPoly>& cc,
                        size_t n,
                        size_t dg_c, size_t df_c, size_t dg_i, size_t df_i, usint depth, usint slots,
                        const char* ksTech) {
    std::cout << "[" << tag << "] n=" << n
              << " | slots=" << slots
              << " | ringDim=" << cc->GetRingDimension()
              << " | depth=" << depth
              << " | (dg_c,df_c,dg_i,df_i)=(" << dg_c << "," << df_c << "," << dg_i << "," << df_i << ")"
              << " | KeySwitch=" << ksTech
              << "\n";
}

template <typename FnRun>
static std::tuple<Stats, Acc, std::array<double,7>>
bench_suite(const std::string& label, size_t n, int trials, const Ctx& ctx, const FnRun& runFn,
            const std::vector<double>& v, const std::vector<double>& v_sorted) {
    // warm-up
    {
        auto ct = runFn(ctx.cc, ctx.pk, v, n);
        Plaintext pt; ctx.cc->Decrypt(ctx.sk, ct, &pt);
        pt->SetLength(n * n);
        (void)pt->GetRealPackedValue();
    }

    std::vector<double> times_s; times_s.reserve(trials);
    double mx = 0.0; long double mseSum = 0.0L;
    std::array<long double,7> stages_sum_s{}; stages_sum_s.fill(0.0L);

    for (int t = 0; t < trials; ++t) {
        stages_reset();

        const auto t0 = std::chrono::steady_clock::now();
        std::chrono::steady_clock::time_point t1;

        {
            auto ct = runFn(ctx.cc, ctx.pk, v, n);
            t1 = std::chrono::steady_clock::now();

            Plaintext pt;
            ctx.cc->Decrypt(ctx.sk, ct, &pt);
            pt->SetLength(n * n);

            std::vector<double> M = pt->GetRealPackedValue();
            std::vector<double> out(n);
            for (size_t i = 0; i < n; ++i) out[i] = M[i * n];

            auto acc = accuracy(out, v_sorted);
            if (acc.max_abs_err > mx) mx = acc.max_abs_err;
            mseSum += acc.mse;

            StageTotals snap = stages_snapshot(); // ms
            for (int k = 0; k < 7; ++k) stages_sum_s[k] += (long double)(snap.ms[k] / 1000.0);
        }

        const double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
        times_s.push_back(elapsed_s);

        // cool between trials to help the allocator release free pages
        CooldownMs(kTrialCooldownMs);
        TrimMalloc();
    }

    Stats st = meanStd(times_s);
    Acc ac{mx, (double)(mseSum / trials)};
    std::array<double,7> stages_mean{};
    for (int k = 0; k < 7; ++k) stages_mean[k] = (double)(stages_sum_s[k] / (long double)trials);

    std::cout << std::fixed << std::setprecision(6)
              << label << " n=" << n << " : mean " << st.mean_s << " s, std " << st.std_s
              << " s | max|err| " << std::setprecision(6) << ac.max_abs_err
              << ", mse " << ac.mse << "\n";

    return {st, ac, stages_mean};
}

#ifdef RUN_FG
static std::vector<Row> run_fg(const std::vector<int>& ns){
    std::vector<Row> rows; rows.reserve(ns.size());
    std::cout << "==== sortFG (Baseline / HYBRID) ====\n";
    for (int n : ns) {
        auto v = loadUniqueCSV((size_t)n);
        auto v_sorted = v; std::sort(v_sorted.begin(), v_sorted.end());

        Ctx ctx = makeCtxFG((size_t)n);
        size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
        paramsFG((size_t)n, dg_c, df_c, dg_i, df_i, depth, slots);
        printParams("FG", ctx.cc, (size_t)n, dg_c, df_c, dg_i, df_i, depth, slots, "HYBRID");

        auto [st, ac, stage_means] = bench_suite("sortFG", (size_t)n, kTrials, ctx, run_sortFG, v, v_sorted);
        rows.push_back({(size_t)n, st, ac, stage_means});

        // strong cleanup between different n
        ctx.cc->ClearEvalAutomorphismKeys();
        ctx.cc->ClearEvalMultKeys();
        ctx.pk.reset();
        ctx.sk.reset();
        ctx.cc.reset();
        HardReset();
        CooldownMs(kBetweenNMs);
    }
    return rows;
}
#endif

#ifdef RUN_FG_LAZY
static std::vector<Row> run_fg_lazy(const std::vector<int>& ns){
    std::vector<Row> rows; rows.reserve(ns.size());
    std::cout << "==== sortFGLazy (Lazy-aware / BATCHED) ====\n";
    for (int n : ns) {
        auto v = loadUniqueCSV((size_t)n);
        auto v_sorted = v; std::sort(v_sorted.begin(), v_sorted.end());

        Ctx ctx = makeCtxFGLazy((size_t)n);
        size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
        paramsFG((size_t)n, dg_c, df_c, dg_i, df_i, depth, slots);
        printParams("FGLazy", ctx.cc, (size_t)n, dg_c, df_c, dg_i, df_i, depth, slots, "BATCHED");

        auto [st, ac, stage_means] = bench_suite("sortFGLazy", (size_t)n, kTrials, ctx, run_sortFGLazy, v, v_sorted);
        rows.push_back({(size_t)n, st, ac, stage_means});

        // strong cleanup between different n
        ctx.cc->ClearEvalAutomorphismKeys();
        ctx.cc->ClearEvalMultKeys();
        ctx.pk.reset();
        ctx.sk.reset();
        ctx.cc.reset();
        HardReset();
        CooldownMs(kBetweenNMs);
    }
    return rows;
}
#endif

int main(){
    std::vector<int> ns(std::begin(kNs), std::end(kNs));

#ifdef RUN_FG
    auto fg = run_fg(ns);
#endif
#ifdef RUN_FG_LAZY
    auto fgl = run_fg_lazy(ns);
#endif

    std::ofstream csv(kCSV);
    csv << std::fixed << std::setprecision(6);
    csv << "n";
#ifdef RUN_FG
    csv << ",fg_mean_s,fg_std_s,fg_max_abs_err,fg_mse";
#endif
#ifdef RUN_FG_LAZY
    csv << ",fgl_mean_s,fgl_std_s,fgl_max_abs_err,fgl_mse";
#endif
#if defined(RUN_FG) && defined(RUN_FG_LAZY)
    csv << ",speedup_fgl_x";
#endif
    csv << ",fg_stage1_mean_s,fg_stage2_mean_s,fg_stage3_mean_s,fg_stage4_mean_s,fg_stage5_mean_s,fg_stage6_mean_s,fg_stage7_mean_s";
    csv << ",fgl_stage1_mean_s,fgl_stage2_mean_s,fgl_stage3_mean_s,fgl_stage4_mean_s,fgl_stage5_mean_s,fgl_stage6_mean_s,fgl_stage7_mean_s";
    csv << ",dg_c,df_c,dg_i,df_i,depth,ringDim,slots\n";

    std::cout << "==== Summary ====\n";
    for (size_t i = 0; i < ns.size(); ++i) {
        int n = ns[i];
        csv << n;

#ifdef RUN_FG
        const auto& a = fg[i];
        csv << "," << a.st.mean_s << "," << a.st.std_s
            << "," << a.ac.max_abs_err << "," << a.ac.mse;
#endif
#ifdef RUN_FG_LAZY
        const auto& b = fgl[i];
        csv << "," << b.st.mean_s << "," << b.st.std_s
            << "," << b.ac.max_abs_err << "," << b.ac.mse;
#endif
#if defined(RUN_FG) && defined(RUN_FG_LAZY)
        {
            double speedup = fg[i].st.mean_s / std::max(1e-12, fgl[i].st.mean_s);
            csv << "," << speedup;
            std::cout << "n=" << n
                      << " | fg "  << fg[i].st.mean_s
                      << " s | fgl " << fgl[i].st.mean_s
                      << " s | speedup " << speedup << "x\n";
        }
#endif

#ifdef RUN_FG
        for (int k = 0; k < 7; ++k) csv << "," << fg[i].stages_mean_s[k];
#else
        for (int k = 0; k < 7; ++k) csv << ",NA";
#endif
#ifdef RUN_FG_LAZY
        for (int k = 0; k < 7; ++k) csv << "," << fgl[i].stages_mean_s[k];
#else
        for (int k = 0; k < 7; ++k) csv << ",NA";
#endif

        size_t dg_c, df_c, dg_i, df_i; usint depth, slots;
        paramsFG((size_t)n, dg_c, df_c, dg_i, df_i, depth, slots);
        csv << "," << dg_c << "," << df_c << "," << dg_i << "," << df_i
            << "," << depth << "," << kRingDim << "," << slots << "\n";
    }
    csv.close();

    std::cout << "Results written to " << kCSV << "\n";

    // final hard reset before exit
    HardReset();
    return 0;
}
