"""Generate null (random genome-wide) variant sets for calibration.

Loads the 1000 Genomes background variant file ONCE (the slow part -- ~83M
rows), then draws random null variant sets by fast in-memory sampling (with
replacement -- collision probability is negligible even for the largest
trait, ~130k draws from an 83M population) for each of the 4 highlighted
traits, matched in size to that trait's real variant count, N_REPS times
each. SNPS column is a synthetic ID (CREST-GV's enrichment computation only
uses CHR_ID/CHR_POS, never SNPS, for peak intersection).

Run from the repo root:
    python scripts/null_calibration_generate.py
"""

import os
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "data", "null_calibration")
BG_FILE = "/project/Wellcome_Discovery/datashare/CREST-GV_collection/1000genomes/ALL_1000_genomes.variants.hg38.bed"

N_REPS = 100
SEED = 20261027  # fixed for reproducibility
TRAITS = {
    "heart_rate": "/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/heart_rate.bed",
    "red_blood_cell_traits": "/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/red_blood_cell_traits.bed",
    "type_1_diabetes_mellitus": "/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_1_diabetes_mellitus.bed",
    "type_2_diabetes_mellitus": "/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_2_diabetes_mellitus.bed",
}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)

    trait_n = {trait: sum(1 for _ in open(path)) for trait, path in TRAITS.items()}
    print("Real variant counts per trait:", trait_n, flush=True)

    print("Loading 1000 Genomes background (one-time cost) ...", flush=True)
    t0 = time.time()
    bg = pd.read_csv(
        BG_FILE, sep="\t", header=None, usecols=[0, 1], names=["chrom", "start"],
        dtype={"chrom": "category", "start": np.int64},
        engine="c",
    )
    chroms = bg["chrom"].to_numpy()
    positions = bg["start"].to_numpy()
    n_bg = len(positions)
    print(f"Loaded {n_bg} background variants in {time.time()-t0:.1f}s", flush=True)

    for trait, n in trait_n.items():
        trait_dir = os.path.join(OUT_DIR, trait)
        os.makedirs(trait_dir, exist_ok=True)
        for rep in range(N_REPS):
            out_path = os.path.join(trait_dir, f"null_rep{rep:03d}.bed")
            if os.path.exists(out_path):
                continue
            idx = rng.integers(0, n_bg, size=n)
            with open(out_path, "w") as f:
                for i, j in enumerate(idx):
                    f.write(f"null_{i}\t{chroms[j]}\t{positions[j]}\n")
        print(f"[{trait}] wrote {N_REPS} null replicates (n={n} variants each) -> {trait_dir}", flush=True)

    print("DONE.", flush=True)


if __name__ == "__main__":
    main()
