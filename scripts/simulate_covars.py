import numpy as np
import pandas as pd

import argparse
parser = argparse.ArgumentParser(description="Simulate environmental covariates.")
parser.add_argument("--fam_file", type=str, required=True, help="Path to the FAM file.")
parser.add_argument("--num_covars", type=int, required=True, help="Number of environmental covariates to simulate.")
args = parser.parse_args()

fam_file = args.fam_file    
num_covars = args.num_covars

with open(fam_file, "r") as f:
    n_individuals = sum(1 for line in f)
n = n_individuals

np.random.seed(42)
covariates = []
# for i in range(num_covars):
#     if i % 2 == 0:
#         cov = np.random.normal(loc=50, scale=10, size=n) / 20
#     else:
#         cov = np.random.binomial(n=1, p=0.48, size=n) / 2
#     covariates.append(cov)

# covariates, some negative, some binary, some uniform
for i in range(num_covars):
    if i % 3 == 0:
        cov = np.random.normal(loc=0, scale=10, size=n) / 20
    elif i % 3 == 1:
        cov = np.random.binomial(n=1, p=0.48, size=n) / 2
    else:
        cov = np.random.uniform(low=5, high=10, size=n) / 10
    covariates.append(cov)

env_covariates = pd.DataFrame({
    f"cov{i+1}": covariates[i] for i in range(num_covars)
})

# cov1 = np.random.normal(loc=50, scale=10, size=N)
# cov2 = np.random.binomial(n=1, p=0.48, size=N)
# cov3 = np.random.binomial(n=1, p=0.3, size=N)
# cov4 = np.random.uniform(low=5, high=30, size=N)
# cov5 = np.random.binomial(n=1, p=0.7, size=N)

# env_covariates = pd.DataFrame({
#     "cov1": cov1,
#     "cov2": cov2,
#     "cov3": cov3,
#     "cov4": cov4,
#     "cov5": cov5
# })

output_file = "simul.cov"

output_dir = "/".join(fam_file.split("/")[:-1])
output_path = f"{output_dir}/{output_file}"

# write column headers: FID, IID, cov names
env_covariates.insert(0, 'FID', range(1, n + 1))
env_covariates.insert(1, 'IID', range(1, n + 1))
env_covariates.columns = ['FID', 'IID'] + [f"cov{i+1}" for i in range(num_covars)]

# each line starts with FID and IID

env_covariates['FID'] = env_covariates['FID'].astype(int)
env_covariates['IID'] = env_covariates['IID'].astype(int)  
env_covariates = env_covariates[['FID', 'IID'] + [f"cov{i+1}" for i in range(num_covars)]]


env_covariates.to_csv(output_path, sep="\t", index=False)
print(f"Wrote simulated environmental covariates to {output_path}.")
