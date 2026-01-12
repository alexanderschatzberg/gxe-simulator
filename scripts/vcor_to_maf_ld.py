#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

def read_afreq(path):
    import numpy as np
    import pandas as pd

    af = pd.read_csv(path, sep=r"\s+")
    id_col = "ID" if "ID" in af.columns else ("SNP" if "SNP" in af.columns else None)
    if id_col is None:
        raise RuntimeError(f"No variant ID column in {path}; columns: {af.columns.tolist()}")

    if "MAF" in af.columns:
        maf = af["MAF"].astype(float).to_numpy()
    elif "ALT_FREQS" in af.columns:
        def maf_from_altfreqs(s):
            parts = [p for p in str(s).split(",") if p not in ("", ".", "NA")]
            if not parts:
                return np.nan
            alts = np.array([float(p) for p in parts], dtype=float)
            alts = alts[(alts >= 0.0) & (alts <= 1.0)]
            sum_alt = float(alts.sum())
            ref = max(0.0, 1.0 - sum_alt)
            return float(np.min(np.append(alts, ref)))
        maf = af["ALT_FREQS"].map(maf_from_altfreqs).to_numpy()
    elif "ALT_FREQ" in af.columns:
        f = af["ALT_FREQ"].astype(float).to_numpy()
        maf = np.minimum(f, 1.0 - f)
    elif "A1_FREQ" in af.columns:   
        f = af["A1_FREQ"].astype(float).to_numpy()
        maf = np.minimum(f, 1.0 - f) 
    else:
        raise RuntimeError(f"{path} lacks MAF/ALT_FREQS/ALT_FREQ/A1_FREQ; columns: {af.columns.tolist()}")

    return pd.DataFrame({"SNP": af[id_col].astype(str), "MAF": maf})

def main():
    ap = argparse.ArgumentParser(description="Aggregate PLINK2 vcor to per-variant LD score and join with MAF.")
    ap.add_argument("--vcor", required=True, help="Tabular vcor file (decompressed).")
    ap.add_argument("--afreq", required=True, help="PLINK2 afreq file.")
    ap.add_argument("--out", required=True, help="Output directory for MAF and LD files.")
    args = ap.parse_args()

    hdr = pd.read_csv(args.vcor, sep=r"\s+", nrows=0)
    cols = set(hdr.columns)

    if {"ID_A","ID_B"}.issubset(cols):
        id_a, id_b = "ID_A", "ID_B"
    elif {"SNP_A","SNP_B"}.issubset(cols):
        id_a, id_b = "SNP_A", "SNP_B"
    else:
        raise RuntimeError(f"Could not find ID_A/ID_B (or SNP_A/SNP_B) in vcor; saw {list(cols)}")

    r2_candidates = ["R2", "UNPHASED_R2", "PHASED_R2"]
    r_candidates  = ["R", "UNPHASED_R", "PHASED_R"]

    ld_name, ld_kind = None, None
    for c in r2_candidates:
        if c in cols:
            ld_name, ld_kind = c, "r2"
            break
    if ld_name is None:
        for c in r_candidates:
            if c in cols:
                ld_name, ld_kind = c, "r"
                break
    if ld_name is None:
        raise RuntimeError(f"Could not find an LD column (R2/PHASED_R2/UNPHASED_R2 or R/PHASED_R/UNPHASED_R). Columns: {list(cols)}")

    # stream and accumulate scores 
    usecols = [id_a, id_b, ld_name]
    ld_sum = {}
    for chunk in pd.read_csv(args.vcor, sep=r"\s+", usecols=usecols, chunksize=1_000_000):
        vals = chunk[ld_name].astype(float).to_numpy()
        if ld_kind == "r":  
            vals = vals * vals
        a = chunk[id_a].astype(str).to_numpy()
        b = chunk[id_b].astype(str).to_numpy()
        for snp, v in zip(a, vals):
            ld_sum[snp] = ld_sum.get(snp, 0.0) + float(v)
        for snp, v in zip(b, vals):
            ld_sum[snp] = ld_sum.get(snp, 0.0) + float(v)

    # join with MAF and add self-term
    af = read_afreq(args.afreq)
    ld_scores = [ld_sum.get(snp, 0.0) + 1.0 for snp in af["SNP"]]
    out = pd.DataFrame({"MAF": af["MAF"].to_numpy(), "LD": np.asarray(ld_scores)})

    import os
    if not os.path.exists(args.out):
        os.makedirs(args.out)
        
    out.to_csv(args.out + "maf_ld.txt", sep="\t", index=False)
    print(f"[done] Wrote {args.out}maf_ld.txt with header: MAF\\tLD")

if __name__ == "__main__":
    main()

