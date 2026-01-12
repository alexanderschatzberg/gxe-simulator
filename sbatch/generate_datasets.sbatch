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
#SBATCH --array=0-6
#SBATCH --export=ALL,UV_CACHE_DIR=/pscratch/sd/q/qys/uv-cache

# ============================================================================
# generate_datasets.sbatch
#
# SLURM array job to generate all 7 profiling datasets in parallel
# Each array task generates one dataset configuration
# Note: Skipping large_N200000_L1 and large_N200000_L4 (too slow)
# ============================================================================

mkdir -p logs

# Explicit configuration mapping
# Format: "label N seq_length L"
declare -A CONFIGS=(
    [0]="small 5000 250000000 1"
    [1]="small 5000 250000000 4"
    [2]="small 5000 250000000 10"
    [3]="medium 50000 1250000000 1"
    [4]="medium 50000 1250000000 4"
    [5]="medium 50000 1250000000 10"
    [6]="large 200000 2500000000 10"
)

# Parse configuration for this task
read -r LABEL N SEQ_LENGTH L <<< "${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

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
