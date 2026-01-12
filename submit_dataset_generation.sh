#!/bin/bash
# ============================================================================
# submit_dataset_generation.sh
#
# Helper script to submit dataset generation jobs to SLURM
# ============================================================================

set -e

# Create logs directory if it doesn't exist
mkdir -p logs

# Make scripts executable
chmod +x generate_single_dataset.sh
chmod +x sbatch/generate_datasets.sbatch

echo "=========================================="
echo "Submitting dataset generation array job"
echo "=========================================="
echo ""
echo "This will generate 27 datasets in parallel:"
echo "  - 3 size levels (small, medium, large)"
echo "  - 3 L values (1, 4, 10 environmental factors)"
echo "  - 3 replicates per configuration"
echo ""
echo "Configurations:"
echo "  small:  N=5,000,    ~100k SNPs,  L=1/4/10, rep=0/1/2"
echo "  medium: N=50,000,   ~500k SNPs,  L=1/4/10, rep=0/1/2"
echo "  large:  N=200,000,  ~1M SNPs,    L=1/4/10, rep=0/1/2"
echo ""
echo "Output directory: sim_data/"
echo "Logs directory: logs/"
echo ""

# Submit the array job
JOB_ID=$(sbatch sbatch/generate_datasets.sbatch | awk '{print $4}')

echo "Submitted array job: $JOB_ID"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/gen_data_${JOB_ID}_*.out"
echo ""
echo "Check status summary:"
echo "  sacct -j $JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS"
echo ""
echo "After completion, verify datasets:"
echo "  find sim_data -name 'genotypes.bed' | wc -l  # should be 27"
echo ""
