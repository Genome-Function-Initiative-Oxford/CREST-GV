"""Fig. 5 -- Analysis of user-provided data with CREST-GV.

Two worked examples, run with the refactored package
(CREST-GV_newCode/src/crestgv), not the precomputed resource:

(a) Held-out collection re-upload: the `cad` collection re-analysed via the
    in-house-collection code path (same bigwigs/peaks, accessed as if a user
    had uploaded them), compared against its precomputed reference-atlas
    result for the same trait (platelet_measurement). Demonstrates that an
    uploaded collection reproduces the atlas result and can be ranked
    alongside it.
(b) User-provided variant set: a GWAS trait not among the four highlighted in
    the main Results (osteoporosis), run fresh against all 11 reference
    collections, exactly as a first-time user upload would be processed.

Run time for each is reported (wall-clock, this run).

Requires CREST-GV_newCode's package on the path and a symlinked in-house shim
for `cad` (see paper_analysis/scripts/build_cad_upload_shim.sh).

Run from the repo root:
    python scripts/fig5_upload_example.py
"""

import os
import sys
import time

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from palette import apply_base_style, SERIES_1, SERIES_2, GRIDLINE, DIVERGING_POS, DIVERGING_NEG

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEWCODE_SRC = os.path.join(os.path.dirname(ROOT), "CREST-GV_newCode", "src")
sys.path.insert(0, NEWCODE_SRC)
from crestgv import crestgv  # noqa: E402

GENETICS_BASE = "/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection"
TRAITS_DIR = "/ceph/project/Wellcome_Discovery/sriva/Git/CREST-GV_v2_backup/ALL_GWAS_CATALOG_TRAITS"
CAD_UPLOAD_SHIM = os.path.join(ROOT, "data", "cad_upload_shim")

PANEL_A_TRAIT = "platelet_measurement"
PANEL_B_TRAIT = "osteoporosis"
ALL_COLLECTIONS = [
    "cad", "catlas_fetal", "catlas_adult", "mpal", "super_pbmc", "calderon",
    "ludwig2019", "erythoid_d7_d10_d13_d17", "h1_hescs", "immune_cell", "greenleaf_brain",
]


def panel_a():
    """Re-upload cad as an in-house collection; compare to the precomputed atlas result."""
    t0 = time.time()
    cgv = crestgv(
        genetic=os.path.join(GENETICS_BASE, f"{PANEL_A_TRAIT}.bed"),
        output="/tmp/fig5_panel_a_out",
        in_house_collection_path=CAD_UPLOAD_SHIM,
    )
    reuploaded = cgv.calculate_enrichment_score()
    runtime_a = time.time() - t0

    atlas = pd.read_csv(
        os.path.join(TRAITS_DIR, PANEL_A_TRAIT, "cad_statistics_CREST-GV.csv"),
        sep="\t", index_col=0,
    )

    common = sorted(set(reuploaded.index) & set(atlas.index))
    compare = pd.DataFrame({
        "atlas_Z": atlas.loc[common, "CREST-GV"],
        "reuploaded_Z": reuploaded.loc[common, "CREST-GV"],
    })
    return compare, runtime_a


def panel_b():
    """Run a novel (not atlas-precomputed) trait fresh against all 11 reference collections."""
    genetic = os.path.join(GENETICS_BASE, f"{PANEL_B_TRAIT}.bed")
    results = {}
    t0 = time.time()
    for coll in ALL_COLLECTIONS:
        cgv = crestgv(genetic=genetic, output="/tmp/fig5_panel_b_out", collection_name=coll)
        df = cgv.calculate_enrichment_score()
        for ct, row in df.iterrows():
            results[f"{coll}_{ct}"] = row["CREST-GV_norm"]
    runtime_b = time.time() - t0
    return pd.Series(results).sort_values(ascending=False), runtime_b


def main():
    apply_base_style()

    print(f"Panel (a): re-uploading cad, trait={PANEL_A_TRAIT} ...")
    compare, runtime_a = panel_a()
    print(f"  done in {runtime_a:.1f}s")

    print(f"Panel (b): fresh run of {PANEL_B_TRAIT} across all 11 collections ...")
    ranked_b, runtime_b = panel_b()
    print(f"  done in {runtime_b:.1f}s")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel a: concordance scatter
    ax1.scatter(compare["atlas_Z"], compare["reuploaded_Z"], s=60, color=SERIES_1, edgecolor="white", linewidth=0.5, zorder=3)
    lims = [compare.values.min() - 0.5, compare.values.max() + 0.5]
    ax1.plot(lims, lims, color=GRIDLINE, linewidth=1.5, linestyle="--", zorder=1)
    for ct in compare.index:
        ax1.annotate(ct, (compare.loc[ct, "atlas_Z"], compare.loc[ct, "reuploaded_Z"]),
                     fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax1.set_xlabel("Reference-atlas raw Z (precomputed)")
    ax1.set_ylabel("Re-uploaded raw Z (this run)")
    ax1.set_title(f"a   Held-out collection re-upload (cad, {PANEL_A_TRAIT})\nrun time: {runtime_a:.1f}s", loc="left", fontsize=10)

    # Panel b: top 10 for the fresh user-trait run
    top10 = ranked_b.head(10)[::-1]
    colors = [DIVERGING_POS if v >= 0 else DIVERGING_NEG for v in top10.values]
    ax2.barh(range(len(top10)), top10.values, color=colors, height=0.6)
    ax2.set_yticks(range(len(top10)))
    ax2.set_yticklabels(top10.index, fontsize=8)
    ax2.axvline(0, color=GRIDLINE, linewidth=1)
    ax2.set_xlabel(r"$\tilde{Z}$ (this run)")
    ax2.set_title(f"b   User-provided variant set ({PANEL_B_TRAIT})\ntop 10 of fresh 11-collection run, run time: {runtime_b:.1f}s", loc="left", fontsize=10)

    fig.tight_layout()
    out_base = os.path.join(ROOT, "figures", "Fig5_upload_example")
    fig.savefig(out_base + ".png")
    fig.savefig(out_base + ".pdf")
    print("wrote", out_base + ".png/.pdf")

    compare.to_csv(os.path.join(ROOT, "tables", "Fig5a_reupload_concordance.csv"))
    ranked_b.to_frame("CREST-GV_norm").to_csv(os.path.join(ROOT, "tables", f"Fig5b_{PANEL_B_TRAIT}_ranked.csv"))
    with open(os.path.join(ROOT, "tables", "Fig5_runtimes.txt"), "w") as f:
        f.write(f"Panel a (re-upload cad, {PANEL_A_TRAIT}): {runtime_a:.2f}s\n")
        f.write(f"Panel b (fresh {PANEL_B_TRAIT}, 11 collections): {runtime_b:.2f}s\n")


if __name__ == "__main__":
    main()
