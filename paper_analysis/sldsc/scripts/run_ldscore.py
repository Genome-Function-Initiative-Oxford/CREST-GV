#!/usr/bin/env python
"""Drive the combine-annots + ldsc.py --l2 steps across chromosomes in parallel.

For each chromosome: build the wide 729-column annot file (if missing), then run
ldsc.py --l2 against it (if missing). Both steps are idempotent so this script can
be re-run safely after a partial failure.
"""
import argparse
import multiprocessing as mp
import os
import subprocess
import sys

SLDSC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLINK_DIR = "/project/Wellcome_Discovery/adelaney/Projects/02_sLDSC/00_ref_files/hg38/1000G_EUR_plink"
PLINK_PREFIX = "1000G.EUR.hg38"
HM3_SNPS = "/project/Wellcome_Discovery/adelaney/Projects/02_sLDSC/00_ref_files/HapMap3_SNPs/list.txt"


def run_chrom(chrom, annot_dir, combined_dir):
    annot_path = os.path.join(combined_dir, f"combined.{chrom}.annot.gz")
    if not os.path.exists(annot_path):
        subprocess.run(
            [sys.executable, os.path.join(SLDSC_DIR, "scripts", "combine_annots.py"),
             "--annot-root", annot_dir, "--out-dir", combined_dir, "--chroms", chrom],
            check=True, cwd=SLDSC_DIR,
        )

    ldscore_path = os.path.join(combined_dir, f"combined.{chrom}.l2.ldscore.gz")
    if os.path.exists(ldscore_path):
        return chrom, "already done"

    bfile = os.path.join(PLINK_DIR, f"{PLINK_PREFIX}.{chrom}")
    out_prefix = os.path.join(combined_dir, f"combined.{chrom}")
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    result = subprocess.run(
        [sys.executable, os.path.join(SLDSC_DIR, "ldsc", "ldsc.py"), "--l2",
         "--bfile", bfile, "--ld-wind-cm", "1",
         "--annot", annot_path, "--thin-annot",
         "--print-snps", HM3_SNPS, "--out", out_prefix],
        cwd=SLDSC_DIR, env=env, capture_output=True, text=True,
    )
    if result.returncode != 0:
        return chrom, f"FAILED:\n{result.stdout[-3000:]}\n{result.stderr[-3000:]}"
    return chrom, "done"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annot-dir", default=os.path.join(SLDSC_DIR, "annot"))
    parser.add_argument("--combined-dir", default=os.path.join(SLDSC_DIR, "combined_annot"))
    parser.add_argument("--chroms", nargs="+", required=True)
    parser.add_argument("--nproc", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.combined_dir, exist_ok=True)
    with mp.Pool(args.nproc) as pool:
        for chrom, status in pool.starmap(
            run_chrom, [(c, args.annot_dir, args.combined_dir) for c in args.chroms]
        ):
            print(f"chr{chrom}: {status}", flush=True)


if __name__ == "__main__":
    main()
