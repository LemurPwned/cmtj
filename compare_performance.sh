#!/bin/bash
# Compare performance between two git commits
# Usage: ./compare_performance.sh [baseline_commit] [optimized_commit]
#   If no commits specified, compares HEAD~1 with HEAD

set -e

BASELINE=${1:-HEAD~1}
OPTIMIZED=${2:-HEAD}
WORKSPACE="/workspace"
RESULTS_DIR="$WORKSPACE/benchmark_results"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}CMTJ Performance Comparison${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "Baseline:  ${YELLOW}$BASELINE${NC}"
echo -e "Optimized: ${YELLOW}$OPTIMIZED${NC}"
echo ""

# Create results directory
mkdir -p "$RESULTS_DIR"

# Save current branch/commit
CURRENT=$(git -C "$WORKSPACE" rev-parse HEAD)
echo -e "${GREEN}Saving current state...${NC}"

# Function to cleanup on exit
cleanup() {
    echo -e "\n${GREEN}Restoring original state...${NC}"
    git -C "$WORKSPACE" checkout "$CURRENT" 2>/dev/null || true
    git -C "$WORKSPACE" submodule update --init --recursive 2>/dev/null || true
}
trap cleanup EXIT

# Function to build and benchmark a commit
benchmark_commit() {
    local COMMIT=$1
    local NAME=$2
    
    echo -e "\n${BLUE}=== Building $NAME ($COMMIT) ===${NC}"
    
    # Checkout commit
    git -C "$WORKSPACE" checkout "$COMMIT"
    git -C "$WORKSPACE" submodule update --init --recursive
    
    # Uninstall previous version
    python3 -m pip uninstall -y cmtj 2>/dev/null || true
    
    # Build and install
    echo "Building..."
    cd "$WORKSPACE"
    python3 -m pip install -e . -q
    
    if [ $? -ne 0 ]; then
        echo -e "${YELLOW}Warning: Build had warnings, but continuing...${NC}"
    fi
    
    # Run benchmark
    echo -e "${GREEN}Running benchmarks...${NC}"
    python3 "$WORKSPACE/quick_benchmark.py" | tee "$RESULTS_DIR/${NAME}_results.txt"
    
    # Extract total time for comparison
    grep "Total benchmark time:" "$RESULTS_DIR/${NAME}_results.txt" | \
        awk '{print $4}' | sed 's/s$//' > "$RESULTS_DIR/${NAME}_total.txt"
}

# Benchmark baseline
benchmark_commit "$BASELINE" "baseline"

# Benchmark optimized
benchmark_commit "$OPTIMIZED" "optimized"

# Compare results
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}COMPARISON${NC}"
echo -e "${BLUE}========================================${NC}"

BASELINE_TIME=$(cat "$RESULTS_DIR/baseline_total.txt")
OPTIMIZED_TIME=$(cat "$RESULTS_DIR/optimized_total.txt")

# Calculate speedup using awk for floating point arithmetic
SPEEDUP=$(awk "BEGIN {printf \"%.2f\", $BASELINE_TIME / $OPTIMIZED_TIME}")

echo -e "Baseline total:  ${YELLOW}${BASELINE_TIME}s${NC}"
echo -e "Optimized total: ${GREEN}${OPTIMIZED_TIME}s${NC}"
echo -e "Speedup:         ${GREEN}${SPEEDUP}x${NC}"

# Calculate percentage improvement
IMPROVEMENT=$(awk "BEGIN {printf \"%.1f\", (($BASELINE_TIME - $OPTIMIZED_TIME) / $BASELINE_TIME) * 100}")
echo -e "Improvement:     ${GREEN}${IMPROVEMENT}%${NC}"

echo -e "\n${BLUE}Detailed results saved in: ${NC}$RESULTS_DIR"
echo -e "${BLUE}========================================${NC}"
