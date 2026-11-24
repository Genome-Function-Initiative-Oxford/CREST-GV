#!/usr/bin/env python
"""Split each combined.<chr>.l2.ldscore.gz (CHR, SNP, BP + 729 per-cell-type L2
columns) into 729 per-cell-type, per-chromosome files, matching the file-prefix
layout ldsc.py's --ref-ld-chr-cts expects (one independently-readable prefix per
cell type): by_celltype/<name>/<name>.<chr>.l2.ldscore.gz, plus matching
.l2.M / .l2.M_5_50 files sliced from the combined .l2.M / .l2.M_5_50 lines
(column order in those files matches the ldscore columns' order, confirmed by
inspection).
"""
import argparse
import gzip
import os

import pandas as pd

CHROMS = [str(c) for c in range(1, 23)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--combined-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--chroms", nargs="+", default=CHROMS)
    args = parser.parse_args()

    for chrom in args.chroms:
        prefix = os.path.join(args.combined_dir, f"combined.{chrom}")
        df = pd.read_csv(f"{prefix}.l2.ldscore.gz", sep="\t", compression="gzip")
        names = [c[:-2] for c in df.columns if c.endswith("L2")]  # strip 'L2'

        m_vals = [float(x) for x in open(f"{prefix}.l2.M").readline().split()]
        m5_vals = [float(x) for x in open(f"{prefix}.l2.M_5_50").readline().split()]
        assert len(names) == len(m_vals) == len(m5_vals), (len(names), len(m_vals), len(m5_vals))

        for i, name in enumerate(names):
            ct_dir = os.path.join(args.out_dir, name)
            os.makedirs(ct_dir, exist_ok=True)
            out_prefix = os.path.join(ct_dir, name)

            ldscore_path = f"{out_prefix}.{chrom}.l2.ldscore.gz"
            if not os.path.exists(ldscore_path):
                sub = df[["CHR", "SNP", "BP", f"{name}L2"]]
                with gzip.open(ldscore_path, "wt") as f:
                    sub.to_csv(f, sep="\t", index=False)

            with open(f"{out_prefix}.{chrom}.l2.M", "w") as f:
                f.write(f"{m_vals[i]:g}\n")
            with open(f"{out_prefix}.{chrom}.l2.M_5_50", "w") as f:
                f.write(f"{m5_vals[i]:g}\n")

        print(f"chr{chrom}: split {len(names)} cell types", flush=True)


if __name__ == "__main__":
    main()
