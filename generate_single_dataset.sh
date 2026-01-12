#!/bin/bash
set -e  # Exit on error
set -u  # Exit on undefined variable

# ============================================================================
# generate_single_dataset.sh
#
# Generates a single GxE profiling dataset
# Designed to be called from SLURM array job
#
# Usage: ./generate_single_dataset.sh <label> <N> <seq_length> <L> <rep>
# ============================================================================

if [ $# -ne 5 ]; then
    echo "Usage: $0 <label> <N> <seq_length> <L> <rep>"
    echo "Example: $0 small 5000 250000000 1 0"
    exit 1
fi

LABEL=$1
N=$2
SEQ_LENGTH=$3
L=$4
REP=$5

# Configuration
BASE_DIR="sim_data"
RECOMB_RATE="1e-7"
MUTATION_RATE="1e-8"

# ============================================================================
# Helper Functions
# ============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Generate VCF using msprime
generate_vcf() {
    local n_individuals=$1
    local seq_length=$2
    local seed=$3
    local output_vcf=$4

    log "Generating VCF: N=$n_individuals, L=$seq_length, seed=$seed"

    cat > /tmp/sim_genotype_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py <<EOF
import msprime
from pathlib import Path
import sys

n_individuals = $n_individuals
sequence_length = $seq_length
mutation_rate = $MUTATION_RATE
recombination_rate = $RECOMB_RATE
random_seed = $seed
vcf_path = Path("$output_vcf")

ts = msprime.sim_ancestry(
    samples=n_individuals,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    population_size=10_000,
    ploidy=2,
    random_seed=random_seed
)

ts = msprime.sim_mutations(
    ts,
    rate=mutation_rate,
    random_seed=random_seed + 1
)

print(f"Simulated {ts.num_individuals} individuals", file=sys.stderr)
print(f"Simulated {ts.num_sites} variants", file=sys.stderr)

individual_names = [f"indiv{i+1}" for i in range(ts.num_individuals)]
vcf_path.parent.mkdir(parents=True, exist_ok=True)

with open(vcf_path, "w") as f:
    ts.write_vcf(f, individual_names=individual_names)

print(f"VCF written to {vcf_path}", file=sys.stderr)
EOF

    uv run python /tmp/sim_genotype_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
    rm /tmp/sim_genotype_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
}

# Convert VCF to PLINK binary format
vcf_to_plink() {
    local vcf_path=$1
    local output_prefix=$2

    log "Converting VCF to PLINK format"

    # Check if plink exists, download if not
    if [ ! -f "./plink" ]; then
        log "Downloading PLINK..."
        curl -o /tmp/plink_${SLURM_JOB_ID}.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
        unzip -p /tmp/plink_${SLURM_JOB_ID}.zip plink > ./plink
        chmod a+x ./plink
        rm /tmp/plink_${SLURM_JOB_ID}.zip
    fi

    ./plink --vcf "$vcf_path" --maf 0.05 --make-bed --out "$output_prefix" --allow-extra-chr
}

# Generate environment file with L factors
generate_environment() {
    local n_individuals=$1
    local n_env=$2
    local seed=$3
    local output_file=$4

    log "Generating $n_env environmental factors"

    cat > /tmp/gen_env_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($seed)
n = $n_individuals
L = $n_env

# Generate L environmental factors, each N(0,1)
env_data = {}
for i in range(L):
    env_data[f"E{i+1}"] = np.random.normal(0, 1, size=n)

env_df = pd.DataFrame(env_data)
env_df.to_csv("$output_file", index=False)
print(f"Generated {L} environmental factors for {n} individuals")
EOF

    uv run python /tmp/gen_env_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
    rm /tmp/gen_env_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
}

# Generate covariates
generate_covariates() {
    local n_individuals=$1
    local seed=$2
    local output_file=$3

    log "Generating covariates"

    cat > /tmp/gen_cov_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($seed)
n = $n_individuals
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

cov_df.to_csv("$output_file", index=False)
print(f"Generated {num_covars} covariates for {n} individuals")
EOF

    uv run python /tmp/gen_cov_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
    rm /tmp/gen_cov_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
}

# Generate phenotype (simple additive model for now)
generate_phenotype() {
    local bed_file=$1
    local env_file=$2
    local seed=$3
    local output_file=$4

    log "Generating phenotype"

    # Extract prefix from bed file
    local prefix="${bed_file%.bed}"

    cat > /tmp/gen_pheno_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py <<EOF
import numpy as np
import pandas as pd

np.random.seed($seed)

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
# For profiling purposes, we just need realistic-looking data

# Simulate some genetic effects (would normally read from .bed)
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
pheno_df.to_csv("$output_file", index=False)
print(f"Phenotype written to $output_file")
EOF

    uv run python /tmp/gen_pheno_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
    rm /tmp/gen_pheno_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.py
}

# ============================================================================
# Main Generation
# ============================================================================

# Create dataset directory
dataset_dir="${BASE_DIR}/${LABEL}_N${N}_L${L}_rep${REP}"
log "Generating: $dataset_dir"

mkdir -p "$dataset_dir"

# Set seed based on configuration (ensures reproducibility)
seed=$((1000 * N + 100 * L + 10 * REP + 42))

# Step 1: Generate VCF
vcf_file="${dataset_dir}/genotype.vcf"
if [ ! -f "$vcf_file" ]; then
    generate_vcf "$N" "$SEQ_LENGTH" "$seed" "$vcf_file"
else
    log "VCF already exists, skipping"
fi

# Step 2: Convert to PLINK format
bed_file="${dataset_dir}/genotypes.bed"
if [ ! -f "$bed_file" ]; then
    vcf_to_plink "$vcf_file" "${dataset_dir}/genotypes"
else
    log "PLINK files already exist, skipping"
fi

# Step 3: Generate environment
env_file="${dataset_dir}/environment.csv"
if [ ! -f "$env_file" ]; then
    generate_environment "$N" "$L" "$((seed + 1000))" "$env_file"
else
    log "Environment file already exists, skipping"
fi

# Step 4: Generate covariates
cov_file="${dataset_dir}/covariates.csv"
if [ ! -f "$cov_file" ]; then
    generate_covariates "$N" "$((seed + 2000))" "$cov_file"
else
    log "Covariates file already exists, skipping"
fi

# Step 5: Generate phenotype
pheno_file="${dataset_dir}/phenotype.csv"
if [ ! -f "$pheno_file" ]; then
    generate_phenotype "$bed_file" "$env_file" "$((seed + 3000))" "$pheno_file"
else
    log "Phenotype file already exists, skipping"
fi

# Optional: Clean up intermediate VCF to save space
rm -f "$vcf_file"

log "Completed: $dataset_dir"

# Print summary
log "Dataset summary:"
log "  Genotypes: $(wc -l < ${dataset_dir}/genotypes.bim) SNPs"
log "  Individuals: $(wc -l < ${dataset_dir}/genotypes.fam)"
log "  Env factors: $L"
log "  Phenotype observations: $(tail -n +2 ${dataset_dir}/phenotype.csv | wc -l)"
