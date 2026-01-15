# ============================================================
# Batched Key-Switching Benchmark Script
# ============================================================
#
# This script runs batched key-switching benchmarks comparing
# HYBRID (baseline) vs BATCHED (lazy) key-switching techniques.
#
# PARAMETER SETS:
#   - N=2^14: ringDim=16384, multDepth=7, scalingBits=34, firstModBits=46
#   - N=2^15: ringDim=32768, multDepth=13, scalingBits=40, firstModBits=51
#   - N=2^16: ringDim=65536, multDepth=24, scalingBits=45, firstModBits=56
#
# EXPERIMENTS:
#   - Rotation counts k: {2, 4, 8, 16, 32, 64}
#   - 100 trials per (parameter set, k) combination
#   - OMP_NUM_THREADS = 1 (single-threaded) and 8 (can be changed depending on system) (multi-threaded)
#
# OUTPUT:
#   Results saved in batched_ks_results/ directory:
#     - timings_threads1.csv  (single-threaded results)
#     - timings_threads8.csv (multi-threaded results)
#
#   Each CSV contains for every trial:
#     - Parameter info: preset, N, depth_L, scalingBits, firstModBits
#     - Modulus info: P_k, sum_logP_bits, P_bits (individual), Q_l, sum_logQ_bits, Q_bits
#     - Timing: ms_baseline, ms_batched, speedup
#     - KS phases: modup/inner/moddown times and call counts
#     - Accuracy: err_baseline, err_batched, err_between
#     - Memory: ciphertext and rotation key sizes
#
# USAGE:
#   chmod +x run_batched_ks_experiments.sh
#   ./run_batched_ks_experiments.sh
#
#
# ============================================================

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPENFHE_LIB="${SCRIPT_DIR}/../openfhe-development/build/lib"
BUILD_DIR="${SCRIPT_DIR}/build"
EXECUTABLE="${BUILD_DIR}/measure-batched-ks"

TRIALS=100

echo "========================================"
echo "Batched Key-Switching Benchmark Suite"
echo "========================================"
echo ""
echo "Configuration:"
echo "  - Trials per (preset, k): ${TRIALS}"
echo "  - Parameter sets: N=2^14, N=2^15, N=2^16"
echo "  - Rotation counts k: {2, 4, 8, 16, 32, 64}"
echo "  - Thread configurations: 1, 8"
echo ""


# Create results directory
RESULTS_DIR="${SCRIPT_DIR}/batched_ks_results"
mkdir -p "${RESULTS_DIR}"

# ============================================================
# Build measure-batched-ks
# ============================================================
echo "[1/3] Building measure-batched-ks..."
cd "${SCRIPT_DIR}"
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake .. > /dev/null 2>&1
make -j > /dev/null 2>&1
echo "  Build complete"
echo ""

if [ ! -f "${EXECUTABLE}" ]; then
    echo "Error: ${EXECUTABLE} not found after build"
    exit 1
fi

echo "Waiting 5 seconds for system stabilization..."
sleep 5
echo ""

# ============================================================
# Single-threaded experiment (OMP_NUM_THREADS=1)
# ============================================================
echo "[2/3] Running single-threaded experiment (OMP_NUM_THREADS=1)..."
echo "  This will run 3 parameter sets × 6 k-values × ${TRIALS} trials"
echo ""

cd "${BUILD_DIR}"
OMP_NUM_THREADS=1 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./measure-batched-ks ${TRIALS}

if [ -f "timings.csv" ]; then
    mv timings.csv "${RESULTS_DIR}/timings_threads1.csv"
    echo ""
    echo "  Single-threaded results saved to: ${RESULTS_DIR}/timings_threads1.csv"
else
    echo "  Warning: timings.csv not found"
fi

echo "  Cooling down for 20 seconds..."
sleep 20
echo ""

# ============================================================
# Multi-threaded experiment (OMP_NUM_THREADS=8)
# ============================================================
echo "[3/3] Running multi-threaded experiment (OMP_NUM_THREADS=8)..."
echo "  This will run 3 parameter sets × 6 k-values × ${TRIALS} trials"
echo "  Estimated time: 15-45 minutes"
echo ""

cd "${BUILD_DIR}"
OMP_NUM_THREADS=8 LD_LIBRARY_PATH="${OPENFHE_LIB}:$LD_LIBRARY_PATH" ./measure-batched-ks ${TRIALS}

if [ -f "timings.csv" ]; then
    mv timings.csv "${RESULTS_DIR}/timings_threads8.csv"
    echo ""
    echo "  Multi-threaded results saved to: ${RESULTS_DIR}/timings_threads8.csv"
else
    echo "  Warning: timings.csv not found"
fi

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
echo "  - timings_threads1.csv   (single-threaded: OMP_NUM_THREADS=1)"
echo "  - timings_threads8.csv   (multi-threaded: OMP_NUM_THREADS=8)"
echo ""
echo "Each CSV contains ${TRIALS} trials for:"
echo "  - 3 parameter sets (N=2^14, 2^15, 2^16)"
echo "  - 6 rotation counts k (2, 4, 8, 16, 32, 64)"
echo ""
echo "CSV columns include:"
echo "  - Timing: ms_baseline, ms_batched, speedup"
echo "  - Moduli: P_k, P_bits, Q_l, Q_bits (showing q_i and p_i values used)"
echo "  - KS phases: modup/inner/moddown times and call counts"
echo "  - Accuracy: err_baseline, err_batched, err_between"
echo "  - Memory: ciphertext and rotation key sizes"
echo ""
