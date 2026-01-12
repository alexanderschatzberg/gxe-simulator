# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project simulates and profiles Gene-Environment (GxE) interaction analyses using the GENIE software. The pipeline generates synthetic genotype data, creates annotations and environmental factors, simulates phenotypes with GxE effects, and profiles GENIE's performance with various thread counts.

## Commands

### Python Environment

The project uses `uv` for Python package management with Python 3.14:

```bash
# Run Python scripts (use uv run)
uv run python scripts/simulate_genotype.py --n_individuals 1000 --seq_length 1000000 --recomb_rate 1e-7 --vcf_path /path/to/output.vcf
```

### Generate Profiling Datasets

Generate standardized datasets for profiling with varying N (individuals), M (SNPs), and L (environmental factors):

**On Perlmutter (recommended for large datasets):**
```bash
# Submit SLURM array job to generate all 7 datasets in parallel
./submit_dataset_generation.sh

# Monitor progress
squeue -u $USER
tail -f logs/gen_data_*.out

# Check completion status
sacct -j <JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

**On local machine:**
```bash
# Generate all profiling datasets (7 total, runs serially)
./generate_profiling_data.sh

# Or generate a single dataset
./generate_single_dataset.sh small 5000 250000000 1
```

This creates datasets in `sim_data/` with structure: `sim_data/{label}_N{N}_L{L}/`

Dataset sizes:
- **small**: N=5,000, ~100k SNPs (seq_length=250M), L=1/4/10
- **medium**: N=50,000, ~500k SNPs (seq_length=1.25G), L=1/4/10
- **large**: N=200,000, ~1M SNPs (seq_length=2.5G), **L=10 only** (L=1/4 skipped due to runtime)

Each dataset contains:
- `genotypes.bed/bim/fam` - PLINK binary format (MAF > 0.05)
- `phenotype.csv` - Simulated phenotype (header: y)
- `environment.csv` - Environmental factors (header: E1,E2,...,EL)
- `covariates.csv` - Covariates (header: C1,C2,C3)

**SLURM array job details:**
- Runs 7 tasks in parallel (array=0-6)
- Each task: 1 node, 16 CPUs, 64GB RAM, 24h time limit
- Uses UV_CACHE_DIR for shared uv package cache
- VCF files deleted after conversion to save space
- Task mapping: 0-2 (small L=1/4/10), 3-5 (medium L=1/4/10), 6 (large L=10)

### Make Targets

All scripts should be run from the project root directory. The Makefile orchestrates the entire simulation and profiling pipeline:

```bash
# Run full pipeline with custom parameters
make N_INDIVIDUALS=100000 SEQ_LENGTH=10000000 all

# Download PLINK tools
make plink
make plink2

# Clean generated data
make clean
```

The `all` target generates flame graphs for GENIE runs with 1, 2, 4, 8, 16, 32, and 64 threads.

### SLURM Batch Jobs

```bash
# Submit profiling job array to SLURM
sbatch sbatch/profile_genie.sbatch
```

The sbatch script runs 4 parameter combinations: varying N_INDIVIDUALS (100k, 200k) and SEQ_LENGTH (100M, 200M).

## Architecture

### Pipeline Stages

The Makefile defines a multi-stage pipeline with explicit dependencies:

1. **Genotype Simulation** (`simulate_genotype.py`)
   - Uses `msprime` to simulate coalescent genealogies with recombination
   - Outputs VCF format genotypes
   - Configurable: population size (10k), mutation rate (1e-8), recombination rate

2. **PLINK Format Conversion**
   - Converts VCF to PLINK binary format (.bed/.bim/.fam)
   - Filters variants by MAF > 0.05

3. **Allele Frequency & LD Calculation**
   - `plink2 --freq`: generates `.afreq` file
   - `plink2 --r2-phased`: computes pairwise LD in `.vcor` format
   - LD window: 1000kb, unlimited variant pairs, all R² values (no threshold)

4. **MAF-LD Score Aggregation** (`vcor_to_maf_ld.py`)
   - Streams large vcor files in chunks (1M rows)
   - Computes per-variant LD scores (sum of R² with all other variants + self)
   - Joins with MAF data
   - Output: `maf_ld.txt` with MAF and LD columns

5. **Annotations, Environment & Parameters** (`simulate_annot_env_param.py`)
   - Binary annotations: 20% of variants marked as functional
   - Environmental factor: single N(0,1) variable per individual
   - GxE simulation parameters: p_causal=0.1, total_h2=0.01, runs 3 sets of 5 simulations each

6. **Covariate Simulation** (`simulate_covars.py`)
   - Generates 5 covariates with mixed distributions (normal, binary, uniform)
   - Output format: tab-separated with FID/IID columns

7. **Phenotype Simulation** (external `Simulator_gxe`)
   - C++ tool (not in repo, under `external/Simulator/build/`)
   - Uses genotypes, annotations, environment, MAF-LD scores
   - Parameters: k=10 (genetic components), jn=50 (jackknife replicates)
   - Outputs phenotypes to `pheno_gxe/` directory

8. **GENIE Profiling**
   - C++ tool (not in repo, under `external/GENIE/build/`)
   - Runs with model: `G+GxE+NxE` (genetic + gene-environment + noise-environment)
   - Parameters: 10 eigenvectors, 50 jackknife replicates
   - Profiled using `perf record` with DWARF call graphs at 100 Hz
   - Generates flame graphs via FlameGraph toolkit (expected under `external/FlameGraph/`)

### Key Data Flow

- All intermediate and output files go to: `/pscratch/sd/q/qys/genie/data_{N_INDIVIDUALS}_{SEQ_LENGTH}/`
- The Makefile uses `.SECONDARY` to preserve intermediate profiling files
- Scripts expect to be invoked with absolute paths or from project root

### Dependencies

**Python packages** (in pyproject.toml):
- msprime: coalescent simulation
- numpy: numerical operations
- pandas: data manipulation

**External tools** (not included):
- PLINK 1.9: VCF to binary format conversion, MAF filtering
- PLINK 2.0: allele frequency and LD calculations
- Simulator_gxe: phenotype simulation with GxE effects
- GENIE: gene-environment interaction analysis
- FlameGraph: performance visualization (stackcollapse-perf.pl, flamegraph.pl)
- perf: Linux profiling tool

**System requirements**:
- Linux (references paths like `/global/homes/`, `/usr/bin/time`, perf tools)
- Expected to run on NERSC compute infrastructure (references SLURM, pscratch filesystem)

### Script Details

**vcor_to_maf_ld.py**: Handles multiple vcor formats (ID_A/ID_B vs SNP_A/SNP_B, R2 vs R variants). Critical for large-scale LD calculations as it streams data to avoid memory issues.

**simulate_annot_env_param.py**: The param.gxe.txt file has header `p_casual ld_ex maf_ex min_maf max_maf total_h2 num_simul` and writes 3 identical parameter rows.
