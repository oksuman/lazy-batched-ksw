#!/bin/bash

# ============================================================
# Matrix Multiplication Benchmark Script for param1 (N=2^14)
# ============================================================
#
# This script runs matrix multiplication benchmarks for both
# ar24 and jkls18 algorithms with the following settings:
#
# PARAMETER SET (param1):
#   - Ring dimension: N = 2^14 (16384)
#   - Multiplicative depth: L = 7
#   - Scaling mod size: 34 bits
#   - First mod size: 46 bits
#   - dnum (NumLargeDigits): 3
#   - Security level: HEStd_128_classic
#
# EXPERIMENTS:
#   - AR24: baseline (HYBRID) + lazy (BATCHED)
#   - JKLS18: baseline + lazy + hoisting
#   - Each with OMP_NUM_THREADS = 1 and 32
#   - Matrix dimensions: 8, 16, 32, 64
#
# OUTPUT:
#   Results are saved in param1_results/ directory:
#     - ar24_param1_threads1.csv
#     - ar24_param1_threads32.csv
#     - jkls18_param1_threads1.csv
#     - jkls18_param1_threads32.csv
#
#   Each CSV contains:
#     - Timing statistics (mean_ms, std_ms)
#     - Accuracy metrics (max_abs_err, mse)
#     - Speedup comparisons
#
# USAGE:
#   chmod +x run_param1_experiments.sh
#   ./run_param1_experiments.sh
#
#
# PERFORMANCE ISOLATION:
#   - 15-second cooldown between experiments to minimize interference
#   - Each test includes internal warm-up runs
#   - Consider running in single-user mode for most consistent results
#
# ============================================================

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPENFHE_LIB="${SCRIPT_DIR}/../../openfhe-development/build/lib"

echo "========================================"
echo "Matrix Multiplication Benchmark Suite"
echo "Parameter Set: param1 (N=2^14)"
echo "========================================"
echo ""
echo "Performance isolation settings:"
echo "  - 15s cooldown between experiments"
echo "  - Each test includes warm-up runs"
echo "  - Total expected runtime: ~15-35 minutes"
echo ""
echo "Waiting 5 seconds for system stabilization..."
sleep 5
echo ""

# Create results directory
RESULTS_DIR="${SCRIPT_DIR}/param1_results"
mkdir -p "${RESULTS_DIR}"

# ============================================================
# AR24 Experiments
# ============================================================
echo "[1/6] Building ar24..."
cd "${SCRIPT_DIR}/ar24"
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake .. > /dev/null 2>&1
make -j > /dev/null 2>&1
echo "  ar24 build complete"
echo ""

echo "[2/6] Running ar24 with OMP_NUM_THREADS=1..."
cd "${SCRIPT_DIR}/ar24/build"
OMP_NUM_THREADS=1 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./ar24-test
mv ar24_bench_results.csv "${RESULTS_DIR}/ar24_param1_threads1.csv"
echo "  Results saved to: ${RESULTS_DIR}/ar24_param1_threads1.csv"
echo "  Cooling down for 15 seconds..."
sleep 15
echo ""

echo "[3/6] Running ar24 with OMP_NUM_THREADS=32..."
cd "${SCRIPT_DIR}/ar24/build"
OMP_NUM_THREADS=32 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./ar24-test
mv ar24_bench_results.csv "${RESULTS_DIR}/ar24_param1_threads32.csv"
echo "  Results saved to: ${RESULTS_DIR}/ar24_param1_threads32.csv"
echo "  Cooling down for 15 seconds..."
sleep 15
echo ""

# ============================================================
# JKLS18 Experiments
# ============================================================
echo "[4/6] Building jkls18..."
cd "${SCRIPT_DIR}/jkls18"
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake .. > /dev/null 2>&1
make -j > /dev/null 2>&1
echo "  jkls18 build complete"
echo ""

echo "[5/6] Running jkls18 with OMP_NUM_THREADS=1..."
cd "${SCRIPT_DIR}/jkls18/build"
OMP_NUM_THREADS=1 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./jkls18-test
mv jkls18_bench_results.csv "${RESULTS_DIR}/jkls18_param1_threads1.csv"
echo "  Results saved to: ${RESULTS_DIR}/jkls18_param1_threads1.csv"
echo "  Cooling down for 15 seconds..."
sleep 15
echo ""

echo "[6/6] Running jkls18 with OMP_NUM_THREADS=32..."
cd "${SCRIPT_DIR}/jkls18/build"
OMP_NUM_THREADS=32 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./jkls18-test
mv jkls18_bench_results.csv "${RESULTS_DIR}/jkls18_param1_threads32.csv"
echo "  Results saved to: ${RESULTS_DIR}/jkls18_param1_threads32.csv"
echo ""

# ============================================================
# Summary
# ============================================================
echo "========================================"
echo "All experiments completed!"
echo "========================================"
echo ""
echo "Results directory: ${RESULTS_DIR}"
echo ""
echo "Generated files:"
echo "  - ar24_param1_threads1.csv    (AR24: baseline + lazy, 1 thread)"
echo "  - ar24_param1_threads32.csv   (AR24: baseline + lazy, 32 threads)"
echo "  - jkls18_param1_threads1.csv  (JKLS18: baseline + lazy + hoisting, 1 thread)"
echo "  - jkls18_param1_threads32.csv (JKLS18: baseline + lazy + hoisting, 32 threads)"
echo ""
echo "Note: Each CSV contains results for dimensions: 8, 16, 32, 64"
echo ""
