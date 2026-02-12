import argparse
import msprime
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser(description="Simulate genotype data.")
parser.add_argument("--n_individuals", required=True, type=int, help="Number of individuals to simulate.")
parser.add_argument("--seq_length", required=True, type=int, help="Length of genomic sequence.")
parser.add_argument("--recomb_rate", required=True, type=float, default=1e-7, help="Recombination rate.")
parser.add_argument("--vcf_path", type=str, default=None, help="Output path to the VCF file (legacy mode).")
parser.add_argument("--plink_prefix", type=str, default=None, help="Output PLINK prefix (direct BED/BIM/FAM, bypasses VCF).")
parser.add_argument("--maf_filter", type=float, default=0.05, help="Minor allele frequency filter (default 0.05).")
args = parser.parse_args()

if args.vcf_path is None and args.plink_prefix is None:
    raise ValueError("Must specify either --vcf_path or --plink_prefix")

n_individuals = args.n_individuals
sequence_length = args.seq_length
mutation_rate = 1e-8
recombination_rate = args.recomb_rate
random_seed = 42

print(f"Simulating ancestry: N={n_individuals}, seq_length={sequence_length}")
ts = msprime.sim_ancestry(
    samples=n_individuals,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    population_size=10_000,
    ploidy=2,
    random_seed=random_seed
)

print(f"Adding mutations (rate={mutation_rate})")
ts = msprime.sim_mutations(
    ts,
    rate=mutation_rate,
    random_seed=random_seed + 1
)

print(f"Simulated {ts.num_individuals} individuals")
print(f"Simulated {ts.num_sites} variants (segregating sites)")

if args.vcf_path is not None:
    # Legacy mode: write VCF (slow for large N)
    vcf_path = Path(args.vcf_path)
    individual_names = [f"indiv{i+1}" for i in range(ts.num_individuals)]
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    with open(vcf_path, "w") as f:
        ts.write_vcf(f, individual_names=individual_names)
    print(f"VCF written to {vcf_path}")

if args.plink_prefix is not None:
    # Direct mode: write PLINK BED/BIM/FAM from tree sequence (no VCF intermediate)
    prefix = Path(args.plink_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    n_samples = ts.num_individuals  # diploid individuals
    n_haplotypes = ts.num_samples   # 2 * n_individuals

    # Pre-compute allele frequencies to filter by MAF
    print(f"Computing allele frequencies and filtering MAF > {args.maf_filter}...")
    variant_ids = []
    variant_positions = []
    variant_alleles = []
    genotype_matrix_rows = []

    for var in ts.variants():
        # Compute allele frequency of alt allele
        alt_freq = np.mean(var.genotypes)
        maf = min(alt_freq, 1.0 - alt_freq)
        if maf < args.maf_filter:
            continue

        variant_ids.append(var.site.id)
        variant_positions.append(int(var.site.position))
        variant_alleles.append(var.alleles)

        # Reshape haplotypes to diploid genotypes (sum of two haplotypes per individual)
        haps = var.genotypes.astype(np.int8)
        dosage = haps[0::2] + haps[1::2]  # shape: (n_individuals,)
        genotype_matrix_rows.append(dosage)

    n_variants = len(variant_ids)
    print(f"After MAF filter: {n_variants} variants (removed {ts.num_sites - n_variants})")

    # Write FAM file (one line per individual)
    fam_path = str(prefix) + ".fam"
    with open(fam_path, "w") as f:
        for i in range(n_samples):
            fid = f"indiv{i+1}"
            iid = f"indiv{i+1}"
            f.write(f"{fid} {iid} 0 0 0 -9\n")
    print(f"Wrote {fam_path} ({n_samples} individuals)")

    # Write BIM file (one line per variant)
    bim_path = str(prefix) + ".bim"
    with open(bim_path, "w") as f:
        for idx in range(n_variants):
            chrom = "1"
            pos = variant_positions[idx]
            alleles = variant_alleles[idx]
            ref = alleles[0] if len(alleles) > 0 else "A"
            alt = alleles[1] if len(alleles) > 1 else "G"
            snp_id = f"1:{pos}"
            f.write(f"{chrom}\t{snp_id}\t0\t{pos}\t{alt}\t{ref}\n")
    print(f"Wrote {bim_path} ({n_variants} variants)")

    # Write BED file (PLINK binary genotype format)
    # BED format: magic bytes + packed genotypes (2 bits per genotype, SNP-major)
    # PLINK encoding: 00=homozygous A1, 01=missing, 10=heterozygous, 11=homozygous A2
    # For dosage 0 (hom ref) -> 00, dosage 1 (het) -> 10, dosage 2 (hom alt) -> 11
    bed_path = str(prefix) + ".bed"
    bytes_per_snp = (n_samples + 3) // 4  # ceil(n_samples / 4)

    # PLINK encoding: dosage 0 -> 0b00, 1 -> 0b10, 2 -> 0b11, missing -> 0b01
    plink_encode = np.array([0b00, 0b10, 0b11, 0b01], dtype=np.uint8)

    # Bit shift table for positions within each byte (4 genotypes per byte)
    bit_shifts = np.array([0, 2, 4, 6], dtype=np.uint8)

    with open(bed_path, "wb") as f:
        # Magic bytes: 0x6C, 0x1B, 0x01 (SNP-major mode)
        f.write(bytes([0x6C, 0x1B, 0x01]))

        for row_idx in range(n_variants):
            dosage = genotype_matrix_rows[row_idx]
            dosage_clipped = np.clip(dosage, 0, 2)
            encoded = plink_encode[dosage_clipped]

            # Pad to multiple of 4 for vectorized packing
            padded_len = bytes_per_snp * 4
            padded = np.zeros(padded_len, dtype=np.uint8)
            padded[:n_samples] = encoded
            reshaped = padded.reshape(bytes_per_snp, 4)
            byte_arr = np.bitwise_or.reduce(reshaped << bit_shifts, axis=1)

            f.write(byte_arr.tobytes())

    print(f"Wrote {bed_path} ({n_variants} variants x {n_samples} individuals)")
    print(f"Direct PLINK output complete (VCF bypass)")
