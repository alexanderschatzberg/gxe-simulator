import numpy as np

# file_path = "./maf.ld.txt"
# take cmd line arg
import sys
if len(sys.argv) != 2:
    print("Usage: python compute_mv_maf.py <file_path>")
    sys.exit(1)

file_path = sys.argv[1]

# check if file exists
try:
    with open(file_path, 'r') as f:
        pass
except FileNotFoundError:
    print(f"File {file_path} does not exist.")
    sys.exit(1)

data = np.loadtxt(file_path, skiprows=1)

col1 = data[:, 0]
col2 = data[:, 1]

mean1, var1 = np.mean(col1), np.var(col1, ddof=1)
mean2, var2 = np.mean(col2), np.var(col2, ddof=1)

print(f"MAF simulated - Mean: {mean1:.4f}, Variance: {var1:.4f}")
print(f"LD simulated - Mean: {mean2:.4f}, Variance: {var2:.4f}")


file_path = "../Simulator/example/maf.ld.txt"

data = np.loadtxt(file_path, skiprows=1)

col1 = data[:, 0]
col2 = data[:, 1]

mean1, var1 = np.mean(col1), np.var(col1, ddof=1)
mean2, var2 = np.mean(col2), np.var(col2, ddof=1)

print(f"MAF example - Mean: {mean1:.4f}, Variance: {var1:.4f}")
print(f"LD example - Mean: {mean2:.4f}, Variance: {var2:.4f}")


