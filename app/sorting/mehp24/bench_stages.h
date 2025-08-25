#pragma once
#include <array>
#include <chrono>

namespace benchstages {

enum : int {
    ST_REPLICATE_ROW = 1,
    ST_TRANSPOSE_ROW = 2,
    ST_REPLICATE_COL = 3,
    ST_COMPARE_ADV   = 4,
    ST_SUM_ROWS      = 5,
    ST_INDICATOR     = 6,
    ST_SUM_COLS      = 7
};

struct StageTotals { std::array<double,7> ms{}; };

inline thread_local std::array<double,7> g_ms{};
inline thread_local std::array<bool,7>   g_open{};
inline thread_local std::array<std::chrono::steady_clock::time_point,7> g_t0{};

inline void stages_reset() { g_ms.fill(0.0); g_open.fill(false); }

inline void stages_begin(int id) {
    const int i = id - 1; if (i < 0 || i >= 7) return;
    g_open[i] = true; g_t0[i] = std::chrono::steady_clock::now();
}
inline void stages_end(int id) {
    const int i = id - 1; if (i < 0 || i >= 7 || !g_open[i]) return;
    auto t1 = std::chrono::steady_clock::now(); g_open[i] = false;
    g_ms[i] += std::chrono::duration<double,std::milli>(t1 - g_t0[i]).count();
}

inline StageTotals stages_snapshot() { StageTotals s; s.ms = g_ms; return s; }

} // namespace benchstages
