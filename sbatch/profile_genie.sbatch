#!/bin/bash
#SBATCH --account=m4646
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=256
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --job-name=genie_sim
#SBATCH --output=logs/genie_%A_%a.out
#SBATCH --error=logs/genie_%A_%a.err
#SBATCH --array=0-3
#SBATCH --export=ALL,UV_CACHE_DIR=/pscratch/sd/q/qys/uv-cache

mkdir -p logs

# Format: "n_individuals seq_length"
PARAMS=(
    "100000 100000000"
    "200000 100000000"
    "100000 200000000"
    "200000 200000000"
)

IFS=' ' read -r N_IND SEQ_LEN <<< "${PARAMS[$SLURM_ARRAY_TASK_ID]}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Parameters: N_INDIVIDUALS=$N_IND, SEQ_LENGTH=$SEQ_LEN"
echo "Job started at $(date), timestamp $(date +%s)"
echo "=========================================="

cd $SLURM_SUBMIT_DIR
make N_INDIVIDUALS=$N_IND SEQ_LENGTH=$SEQ_LEN all

echo "Job completed at $(date), timestamp $(date +%s)"
