# CREST-GV paper — analyses and figures

Reproducible analyses and figures behind the CREST-GV manuscript, generated from the
final, QC-filtered, 390-profile/11-collection, re-standardised atlas. These are the
same draft figures used to check numbers and layout while writing the paper, not
necessarily the final print-ready artwork submitted to the journal.

## What's here

```
data/              the score/p-value matrices these scripts read (see Sources below)
scripts/           one script per figure, runnable standalone, writes PNG+PDF to figures/
                   and the underlying numbers to tables/
figures/           rendered PNG+PDF for every figure below, plus a PowerPoint deck
                   (CREST-GV_main_figures.pptx) assembling each with its legend
tables/            CSV of the exact numbers behind every figure, and Table 1
cheers_benchmark/  CHEERS vs CREST-GV benchmark results (summary stats + per-pair data)
null_calibration/  null calibration (random genome-wide variant sets) -- see its own README
```

| Script | Manuscript item | What it shows |
|---|---|---|
| `fig2_standardisation.py` | Fig. 2 | Raw ($Z_{ij}$) vs standardised ($\tilde{Z}_{ij}$) score distributions per collection, for red_blood_cell_traits |
| `fig3_cross_trait_heatmap.py` | Fig. 3 | Standardised scores for the cell types significant (Benjamini-Hochberg FDR q<0.10) in >=1 of the 4 highlighted traits (6 cell types) |
| `fig4_per_trait_rankings.py` | Fig. 4 | Top/bottom 10 standardised scores per trait, 4 panels, FDR-significant hits marked |
| `fig5_upload_example.py` | Fig. 5 | **(a)** The `cad` collection re-analysed via the in-house-upload code path, compared to its precomputed atlas result — concordance is near-exact (differences ~1e-13, floating-point noise only). **(b)** A trait not among the four highlighted in Results (`osteoporosis`) run fresh against all 11 collections, exactly as a first-time upload would be processed. Run times reported for both. |
| `suppfig_qc_distributions.py` | Supp. Fig. 1 | QC prediction-score histogram with Good/Decent/Bad thresholds shaded |
| `suppfig_cheers_scatter.py` | Supp. Fig. 3 | CHEERS-vs-CREST-GV concordance scatter, per collection and pooled |
| `null_calibration_analyze.py` | Supp. Fig. 2 | Null calibration QQ-plot and empirical-FPR-vs-nominal summary (see `null_calibration/README.md`) |
| — | Table 1 | `tables/Table1_collections.csv` |

Colours follow a fixed-order categorical palette (`scripts/palette.py`): one hue for
magnitude, a blue/red diverging pair for signed enrichment scores, a reserved
green/amber/red status palette for the three QC classes. White background, no
gridlines, thin marks — static figures for a print PDF manuscript.

## Sources

- `data/scores_390kept.tsv` / `data/pvals_390kept.tsv`: the final, re-standardised
  390-profile x 1,125-trait matrices (`CREST-GV_norm` / `CREST-GV_pval_norm` restricted to
  QC-kept cell types, standardised with the package's own sample-std, ddof=1,
  convention).
- `data/scores_390kept_from_norm.tsv` / `data/pvals_390kept_from_norm.tsv`: a
  comparison file — the *original* `CREST-GV_norm`/`CREST-GV_pval_norm` values
  (standardised against the full, unfiltered 511-profile atlas), simply filtered down
  to the 390 kept rows, with no re-standardisation. Useful for checking how much the
  re-standardisation basis itself (kept-only vs. full atlas) changes results; not used
  by any of the figure scripts.
- `data/CREST-GV_qc_summary_all.csv`: per-track QC labels (Good/Decent/Bad), used to
  pick the 390 kept cell types.
- `cheers_benchmark/`, `null_calibration/`: outputs of the comparison and calibration
  analyses described in the paper's Results and Supplementary Note 1; each has its own
  README with full methodology.

Two scripts (`fig5_upload_example.py`, `recompute_390kept_ddof1.py`) read from this
cluster's internal paths (the live reference chromatin accessibility collections, a
GWAS variant file directory) and from `../src/crestgv` directly — they're included for
transparency of exactly how the reported numbers were produced, not as portable,
drop-in scripts; expect to need to adapt their paths to rerun them outside this
cluster. All other scripts only need the `data/` files in this folder.

## Regenerating the figures that don't need cluster-internal paths

```bash
pip install -e ..   # installs crestgv itself, used by fig5_upload_example.py
cd paper_analysis
python3 scripts/fig2_standardisation.py
python3 scripts/fig3_cross_trait_heatmap.py
python3 scripts/fig4_per_trait_rankings.py
python3 scripts/suppfig_qc_distributions.py
python3 scripts/suppfig_cheers_scatter.py
python3 scripts/null_calibration_analyze.py
```
