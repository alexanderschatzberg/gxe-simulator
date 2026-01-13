#!/bin/bash
# ============================================================================
# test_uv_setup.sh
#
# Quick test to verify UV cache setup without submitting a SLURM job
# Run this on Perlmutter: ./test_uv_setup.sh
# ============================================================================

set -e

echo "=========================================="
echo "Testing UV setup on Perlmutter"
echo "=========================================="
echo ""

# Simulate SLURM environment
export SLURM_SUBMIT_DIR=$(pwd)
export SLURM_JOB_ID=test_$$
export SLURM_ARRAY_TASK_ID=0

echo "1. Current directory: $(pwd)"
echo ""

# Set up UV cache directory (same as sbatch script)
echo "2. Setting up UV cache directory..."
unset UV_CACHE_DIR
unset UV_TOOL_DIR
unset UV_PYTHON_INSTALL_DIR
export UV_CACHE_DIR=$SLURM_SUBMIT_DIR/.uv-cache
mkdir -p $UV_CACHE_DIR
echo "   UV_CACHE_DIR: $UV_CACHE_DIR"
echo "   Directory exists: $([ -d "$UV_CACHE_DIR" ] && echo 'YES' || echo 'NO')"
echo "   Directory writable: $([ -w "$UV_CACHE_DIR" ] && echo 'YES' || echo 'NO')"
echo ""

# Check for system uv config that might override
echo "3. Checking for system UV config files..."
if [ -f ~/.config/uv/uv.toml ]; then
    echo "   WARNING: Found ~/.config/uv/uv.toml"
    echo "   Contents:"
    cat ~/.config/uv/uv.toml | sed 's/^/     /'
    echo ""
    echo "   This might override UV_CACHE_DIR!"
else
    echo "   No ~/.config/uv/uv.toml found (good)"
fi
echo ""

# Test uv with a simple command
echo "4. Testing UV with a simple Python command..."
echo "   Running: uv run python --version"
if uv run python --version 2>&1; then
    echo "   ✓ SUCCESS: UV is working!"
else
    echo "   ✗ FAILED: UV encountered an error"
    exit 1
fi
echo ""

# Create a temporary test script to verify imports work
echo "5. Testing Python package imports (msprime, numpy, pandas)..."
cat > /tmp/test_imports_$$.py <<'EOF'
import sys
print(f"Python: {sys.version}")
try:
    import msprime
    print(f"msprime: {msprime.__version__}")
except ImportError as e:
    print(f"ERROR importing msprime: {e}")
    sys.exit(1)

try:
    import numpy as np
    print(f"numpy: {np.__version__}")
except ImportError as e:
    print(f"ERROR importing numpy: {e}")
    sys.exit(1)

try:
    import pandas as pd
    print(f"pandas: {pd.__version__}")
except ImportError as e:
    print(f"ERROR importing pandas: {e}")
    sys.exit(1)

print("All imports successful!")
EOF

if uv run python /tmp/test_imports_$$.py 2>&1; then
    echo "   ✓ SUCCESS: All required packages are importable!"
else
    echo "   ✗ FAILED: Package import error"
    rm -f /tmp/test_imports_$$.py
    exit 1
fi
rm -f /tmp/test_imports_$$.py
echo ""

# Test actual msprime simulation (tiny one)
echo "6. Testing minimal msprime simulation..."
cat > /tmp/test_msprime_$$.py <<'EOF'
import msprime
print("Running minimal msprime simulation (10 individuals, 10kb)...")
ts = msprime.sim_ancestry(
    samples=10,
    sequence_length=10000,
    recombination_rate=1e-8,
    population_size=10000,
    random_seed=42
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=43)
print(f"  Generated {ts.num_sites} variants")
print("  Msprime simulation successful!")
EOF

if uv run python /tmp/test_msprime_$$.py 2>&1; then
    echo "   ✓ SUCCESS: Msprime can run simulations!"
else
    echo "   ✗ FAILED: Msprime simulation error"
    rm -f /tmp/test_msprime_$$.py
    exit 1
fi
rm -f /tmp/test_msprime_$$.py
echo ""

echo "=========================================="
echo "✓ All tests passed!"
echo "=========================================="
echo ""
echo "Your environment is correctly configured."
echo "You can now submit jobs with:"
echo "  ./submit_dataset_generation.sh"
echo ""
