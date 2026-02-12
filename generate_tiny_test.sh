#!/bin/bash
set -e  # Exit on error

# ============================================================================
# generate_tiny_test.sh
#
# Generates a TINY test dataset that will run quickly on a laptop
# - N=100 individuals
# - seq_length=100,000 (~100 SNPs after MAF filtering)
# - L=1 environmental factor
#
# Uses direct PLINK output (no VCF intermediate)
# ============================================================================

# Check for required tools
command -v uv >/dev/null 2>&1 || { echo "Error: uv is required but not installed"; exit 1; }

# Configuration
BASE_DIR="sim_data"
LABEL="tiny"
N=100
SEQ_LENGTH=100000
L=1
RECOMB_RATE="1e-7"
MUTATION_RATE="1e-8"

# ============================================================================
# Helper Functions
# ============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ============================================================================
# Main Generation
# ============================================================================

dataset_dir="${BASE_DIR}/${LABEL}_N${N}_L${L}"
log "Generating TINY test dataset: $dataset_dir"
log "This should complete in under a minute on a laptop"
echo ""

mkdir -p "$dataset_dir"

# Set seed for reproducibility
seed=42

# Step 1: Generate PLINK files directly (bypasses VCF)
log "[1/4] Generating genotypes directly to PLINK (N=$N, seq_length=$SEQ_LENGTH)..."
bed_file="${dataset_dir}/genotypes.bed"

uv run python scripts/simulate_genotype.py \
    --n_individuals "$N" \
    --seq_length "$SEQ_LENGTH" \
    --recomb_rate "$RECOMB_RATE" \
    --plink_prefix "${dataset_dir}/genotypes" \
    --maf_filter 0.05
echo ""

# Steps 2-3: Generate environment and covariates IN PARALLEL
log "[2/4] Generating environment and covariates (in parallel)..."
env_file="${dataset_dir}/environment.csv"
cov_file="${dataset_dir}/covariates.csv"

cat > /tmp/gen_env_tiny.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($((seed + 1000)))
n = $N
L = $L

env_data = {}
for i in range(L):
    env_data[f"E{i+1}"] = np.random.normal(0, 1, size=n)

env_df = pd.DataFrame(env_data)
env_df.to_csv("$env_file", index=False)
print(f"Generated {L} environmental factors for {n} individuals")
EOF

cat > /tmp/gen_cov_tiny.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($((seed + 2000)))
n = $N
num_covars = 3

covariates = []
for i in range(num_covars):
    if i % 3 == 0:
        cov = np.random.normal(loc=0, scale=10, size=n) / 20
    elif i % 3 == 1:
        cov = np.random.binomial(n=1, p=0.48, size=n) / 2
    else:
        cov = np.random.uniform(low=5, high=10, size=n) / 10
    covariates.append(cov)

cov_df = pd.DataFrame({
    f"C{i+1}": covariates[i] for i in range(num_covars)
})

cov_df.to_csv("$cov_file", index=False)
print(f"Generated {num_covars} covariates for {n} individuals")
EOF

uv run python /tmp/gen_env_tiny.py &
pid_env=$!
uv run python /tmp/gen_cov_tiny.py &
pid_cov=$!
wait $pid_env
wait $pid_cov
rm /tmp/gen_env_tiny.py /tmp/gen_cov_tiny.py
echo ""

# Step 4: Generate phenotype
log "[3/4] Generating phenotype..."
pheno_file="${dataset_dir}/phenotype.csv"
prefix="${dataset_dir}/genotypes"

cat > /tmp/gen_pheno_tiny.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($((seed + 3000)))

# Read PLINK FAM to get number of individuals
with open("${prefix}.fam", "r") as f:
    n = sum(1 for line in f)

# Read environment
env_df = pd.read_csv("$env_file")
L = env_df.shape[1]

# Read genotypes (simplified: read from BIM to get number of SNPs)
with open("${prefix}.bim", "r") as f:
    m = sum(1 for line in f)

print(f"N={n} individuals, M={m} SNPs, L={L} environmental factors")

# Simple simulation: y = genetic effect + GxE effect + noise
genetic_effect = np.random.normal(0, 1, size=n)

# Simulate GxE effects
gxe_effect = np.zeros(n)
for l in range(L):
    gxe_effect += env_df.iloc[:, l].values * np.random.normal(0, 0.5, size=n)

# Add noise
noise = np.random.normal(0, 1, size=n)

# Total heritability controlled
h2_g = 0.3
h2_gxe = 0.1
h2_total = h2_g + h2_gxe

genetic_effect = genetic_effect / np.std(genetic_effect) * np.sqrt(h2_g)
gxe_effect = gxe_effect / np.std(gxe_effect) * np.sqrt(h2_gxe)
noise = noise / np.std(noise) * np.sqrt(1 - h2_total)

y = genetic_effect + gxe_effect + noise

# Save phenotype
pheno_df = pd.DataFrame({"y": y})
pheno_df.to_csv("$pheno_file", index=False)
print(f"Phenotype written to $pheno_file")
EOF

uv run python /tmp/gen_pheno_tiny.py
rm /tmp/gen_pheno_tiny.py
echo ""

# ============================================================================
# Summary
# ============================================================================

log "[4/4] COMPLETE! Tiny test dataset generated at: $dataset_dir"
echo ""
log "Dataset summary:"
log "  Individuals: $(wc -l < ${dataset_dir}/genotypes.fam | tr -d ' ')"
log "  SNPs (after MAF filter): $(wc -l < ${dataset_dir}/genotypes.bim | tr -d ' ')"
log "  Environmental factors: $L"
log "  Files generated:"
ls -lh "$dataset_dir" | tail -n +2 | awk '{print "    " $9 " (" $5 ")"}'
echo ""
log "To inspect the data:"
log "  head $dataset_dir/phenotype.csv"
log "  head $dataset_dir/environment.csv"
log "  head $dataset_dir/covariates.csv"
