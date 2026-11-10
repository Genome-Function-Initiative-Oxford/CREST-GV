# CREST-GV vs CHEERS comparison

The exact pipeline used to produce the CHEERS-vs-CREST-GV benchmark reported in the
paper (Results, "Comparison with CHEERS"; Supplementary Note 1): runs
[CHEERS](https://github.com/TrynkaLab/CHEERS) (Chromatin Element Enrichment Ranking by
Specificity) on the same chromatin accessibility collections used by CREST-GV, so its
enrichment statistic can be compared against CREST-GV's own, collection by collection
and trait by trait.

## What's here

- `config.py` — shared paths, collection/trait definitions.
- `scripts/`
  - `build_consensus_peaks.py` — merges each collection's per-cell-type LanceOtron
    peak calls into one consensus set (CHEERS requires a single shared peak set).
  - `quantify_bigwigs.py` — quantifies each cell type's signal over the consensus
    peaks from its RPKM-normalised bigWig track.
  - `run_collection.py` / `run_all.py` — runs CHEERS's own normalisation and
    enrichment on the quantified signal.
  - `compare_with_crestgv.py` — compares CHEERS enrichment against CREST-GV's
    standardised statistic (Spearman rank correlation; BH-FDR significance overlap).
- `docs/`
  - `pipeline_report.md` — full write-up of the pipeline design, including the
    derivation for why RPKM-normalised bigWig signal can be used in place of CHEERS's
    expected raw read counts (its normalisation is invariant to per-sample
    multiplicative constants).
  - `cheers_vs_crestgv_scatter.html` / `scatter_data.json` — an interactive scatter
    plot of the comparison.

## What's *not* here

The consensus peaks, quantified signal, and per-collection/per-trait results
(`data/`, `results/`, `logs/` in the original working directory) are not included —
together they're several GB of intermediate output, not meant for version control.
The headline numbers from these results (median Spearman rho, per-collection/per-trait
breakdowns, FDR-significance overlap) are reported in the paper itself (Results +
Supplementary Note 1), not reproduced as raw files here.

## Running this yourself

This pipeline was run on an internal HPC cluster, and `config.py` reflects that
environment directly: absolute paths to the reference chromatin accessibility
collections, GWAS variant files, and specific tool binaries (`bedtools`,
`bigWigAverageOverBed`) and Python interpreter paths particular to that cluster's
module system. It is published here for transparency and reproducibility of the exact
steps taken, not as a portable, drop-in pipeline — expect to need to adapt `config.py`'s
paths (`COLLECTIONS_BASE`, `TRAITS`, `CRESTGV_RESULTS_BASE`, the tool/interpreter
binaries at the bottom of the file) to your own environment before running it.
