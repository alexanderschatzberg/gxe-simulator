import argparse
import msprime
from pathlib import Path

# command-line args for n_individuals, sequence_length, and recombination_rate

parser = argparse.ArgumentParser(description="Simulate genotype data.")
parser.add_argument("--n_individuals", required=True, type=int, help="Number of individuals to simulate.")
parser.add_argument("--seq_length", required=True, type=int, help="Length of genomic sequence.")
parser.add_argument("--recomb_rate", required=True, type=float, default=1e-7, help="Recombination rate.")
parser.add_argument("--vcf_path", type=str, required=True, help="Output path to the VCF file.")
args = parser.parse_args()

n_individuals = args.n_individuals
sequence_length = args.seq_length
mutation_rate = 1e-8
recombination_rate = args.recomb_rate
vcf_path = Path(args.vcf_path)
random_seed = 42

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

print(f"Simulated {ts.num_individuals} individuals")
print(f"Simulated {ts.num_sites} variants (segregating sites)")

individual_names = [f"indiv{i+1}" for i in range(ts.num_individuals)]

vcf_path.parent.mkdir(parents=True, exist_ok=True)

with open(vcf_path, "w") as f:
    ts.write_vcf(f, individual_names=individual_names)

print(f"VCF written to {vcf_path}")
