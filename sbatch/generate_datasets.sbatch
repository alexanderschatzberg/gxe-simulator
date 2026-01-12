#!/bin/bash
#SBATCH --account=m4646
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=gen_gxe_data
#SBATCH --output=logs/gen_data_%A_%a.out
#SBATCH --error=logs/gen_data_%A_%a.err
#SBATCH --array=0-8
#SBATCH --export=ALL,UV_CACHE_DIR=/pscratch/sd/q/qys/uv-cache

# ============================================================================
# generate_datasets.sbatch
#
# SLURM array job to generate all 9 profiling datasets in parallel
# Each array task generates one dataset configuration
# ============================================================================

mkdir -p logs

# Dataset configurations
# Format: "label N seq_length"
declare -A SIZE_CONFIGS=(
    [0]="small 5000 250000000"       # ~100k SNPs
    [1]="medium 50000 1250000000"    # ~500k SNPs
    [2]="large 200000 2500000000"    # ~1M SNPs
)

ENV_COUNTS=(1 4 10)

# Map array task ID to configuration
# 9 total: 3 sizes × 3 L values
SIZE_IDX=$((SLURM_ARRAY_TASK_ID / 3))
L_IDX=$((SLURM_ARRAY_TASK_ID % 3))

read -r LABEL N SEQ_LENGTH <<< "${SIZE_CONFIGS[$SIZE_IDX]}"
L=${ENV_COUNTS[$L_IDX]}

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Configuration: ${LABEL}_N${N}_L${L}"
echo "Parameters: N=$N, SEQ_LENGTH=$SEQ_LENGTH, L=$L"
echo "Job started at $(date), timestamp $(date +%s)"
echo "=========================================="

# Navigate to submission directory
cd $SLURM_SUBMIT_DIR

# Ensure script is executable
chmod +x generate_single_dataset.sh

# Set number of threads for numpy/etc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run generation
./generate_single_dataset.sh "$LABEL" "$N" "$SEQ_LENGTH" "$L"

echo "=========================================="
echo "Job completed at $(date), timestamp $(date +%s)"
echo "=========================================="
