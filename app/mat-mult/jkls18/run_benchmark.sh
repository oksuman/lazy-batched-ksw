#!/bin/bash

# Automated benchmark script for JKLS18 matrix multiplication
# Tests sequential (1 thread) vs parallel (32 threads) execution
# Runs across 3 parameter sets and 4 dimensions

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
LIB_PATH="$SCRIPT_DIR/../../../openfhe-development/build/lib"
RESULTS_DIR="$SCRIPT_DIR/results"

# Create results directory
mkdir -p "$RESULTS_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=====================================${NC}"
echo -e "${BLUE}JKLS18 Matrix Multiplication Benchmark${NC}"
echo -e "${BLUE}=====================================${NC}"
echo ""

# Build the project
echo -e "${YELLOW}Building jkls18-test...${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake .. > cmake.log 2>&1
if [ $? -ne 0 ]; then
    echo -e "${RED}CMake failed! Check cmake.log${NC}"
    exit 1
fi

make -j > make.log 2>&1
if [ $? -ne 0 ]; then
    echo -e "${RED}Make failed! Check make.log${NC}"
    exit 1
fi
echo -e "${GREEN}Build successful!${NC}"
echo ""

# Check if executable exists
if [ ! -f "$BUILD_DIR/jkls18-test" ]; then
    echo -e "${RED}Error: jkls18-test executable not found!${NC}"
    exit 1
fi

# Timestamp for this run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Function to run benchmark with specific thread count
run_benchmark() {
    local THREADS=$1
    local RUN_NAME=$2
    local OUTPUT_FILE="$RESULTS_DIR/jkls18_results_${RUN_NAME}_threads${THREADS}_${TIMESTAMP}.csv"

    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Running with OMP_NUM_THREADS=${THREADS}${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""

    export OMP_NUM_THREADS=$THREADS
    export LD_LIBRARY_PATH="$LIB_PATH:$LD_LIBRARY_PATH"

    # Run the benchmark
    cd "$BUILD_DIR"
    ./jkls18-test 2>&1 | tee "$RESULTS_DIR/jkls18_output_${RUN_NAME}_threads${THREADS}_${TIMESTAMP}.log"

    # Move the generated CSV to results directory with timestamp
    if [ -f "jkls18_bench_results.csv" ]; then
        mv jkls18_bench_results.csv "$OUTPUT_FILE"
        echo -e "${GREEN}Results saved to: $OUTPUT_FILE${NC}"
    else
        echo -e "${RED}Warning: jkls18_bench_results.csv not found!${NC}"
    fi

    echo ""
}

# Run sequential (1 thread)
run_benchmark 1 "sequential"

echo -e "${YELLOW}Waiting 5 seconds before next run...${NC}"
sleep 5

# Run parallel (32 threads)
run_benchmark 32 "parallel"

echo ""
echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}All benchmarks completed!${NC}"
echo -e "${GREEN}=====================================${NC}"
echo ""
echo -e "Results saved in: ${BLUE}$RESULTS_DIR${NC}"
echo ""
echo "Files generated:"
ls -lh "$RESULTS_DIR"/*${TIMESTAMP}* 2>/dev/null || echo "No files found"

echo ""
echo -e "${YELLOW}Summary:${NC}"
echo "  - Sequential (1 thread):  results_sequential_threads1_${TIMESTAMP}.csv"
echo "  - Parallel (32 threads):  results_parallel_threads32_${TIMESTAMP}.csv"
echo ""
echo "You can compare the results to see the speedup from parallelization."
