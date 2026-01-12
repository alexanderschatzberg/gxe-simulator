import numpy as np

# command-line arg to get output dir and bim file

import argparse
parser = argparse.ArgumentParser(description="Simulate annotation and environmental data.")
parser.add_argument("--data_dir", type=str, required=True, help="Directory where PLINK files are located.")
parser.add_argument("--plink_prefix", type=str, required=True, help="Prefix for PLINK files.")
parser.add_argument("--num_simul", type=int, required=True, help="Number of simulations to run.")
args = parser.parse_args()

data_dir = args.data_dir
plink_prefix = args.plink_prefix
num_simul = args.num_simul

bim_file = f"{data_dir}/{plink_prefix}.bim"

# get line count of bim file
with open(bim_file, "r") as f:
    num_sites = sum(1 for line in f)

annotations = np.random.binomial(1, 0.2, size=num_sites)

np.savetxt(f"{data_dir}/annotations.txt", annotations, fmt="%d")
print(f"Wrote binary annotations to {data_dir}annotations.txt")

fam_file = f"{data_dir}{plink_prefix}.fam"
with open(fam_file, "r") as f:
    n_individuals = sum(1 for line in f)

E = np.random.normal(0, 1, size=(n_individuals, 1))
np.savetxt(f"{data_dir}env.txt", E, fmt="%.6f")
print(f"Wrote environmental data for {n_individuals} individuals to {data_dir}env.txt")

p_causal = 0.1
ld_ex = 0.0
maf_ex = 0.0
min_maf = 0.0
max_maf = 0.5  
total_h2 = 0.01

with open(f"{data_dir}/param.gxe.txt", "w") as f:
    f.write(f"p_casual ld_ex maf_ex min_maf max_maf total_h2 num_simul\n")
    for _ in range(3): 
        f.write(f"{p_causal} {ld_ex} {maf_ex} {min_maf} {max_maf} {total_h2} {num_simul}\n")
print(f"Wrote simulation parameters to {data_dir}/param.gxe.txt")

# pretty border string
border = "=" * 50
print(border)

print("Simulation data generation complete.")
