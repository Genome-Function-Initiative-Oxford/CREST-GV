# Null calibration

Checks whether CREST-GV's P-value is well-calibrated, using random genome-wide null
variant sets (see `../scripts/null_calibration_generate.py` and `null_calibration_run.py`).

## Method
For each of the 4 highlighted traits, 100 null variant sets were drawn (with
replacement, fixed seed) from the same 1000 Genomes background CREST-GV already uses,
matched in size to that trait's real variant count (3,796 / 130,668 / 9,673 / 48,891 for
heart rate / red blood cell traits / T1D / T2D respectively). Each null set was run
through the refactored `crestgv` package (this repo's `src/crestgv`) against all 11 atlas collections,
recording raw Z per cell type (`<trait>_null_results.tsv`, 400 files' worth of runs, 4
traits x 100 reps x 11 collections = 4,400 CREST-GV runs). The standardised P-value
(the one the manuscript's $P<\pthr$ significance criterion actually uses) was then
re-derived offline per (replicate, collection), restricted to the 390 QC-kept cell
types — exactly the same recipe as the real atlas's re-standardisation, just applied to
null data (`../scripts/null_calibration_analyze.py`).

## Results

**Per-test calibration is good.** Empirical false-positive rate at nominal P<0.01 ranges
0.0086–0.0116 across the 4 traits (expected 0.01); at P<0.05, ranges 0.049–0.053
(expected 0.05). See `SuppFig_null_calibration.png` panel (a) (QQ-plot, tracks y=x
closely) and panel (b), and `../tables/NullCalibration_fpr_by_trait.csv` /
`_by_collection.csv`.

**But across the whole 390-profile atlas, a null/random trait is almost certain to
produce at least one "significant" hit by chance.** `P(>=1 false hit in one null
replicate)` is 97–99% at P<0.01 and ~100% at P<0.05 (`../tables/NullCalibration_atleastone_hit.csv`).
This is the expected consequence of testing 390 cell types simultaneously with no
multiple-testing correction — each individual test is calibrated, but running that many
of them means *some* false positive is all but guaranteed every time, real signal or not.

**Concrete impact on the four highlighted traits**: a Bonferroni correction for a 5%
family-wise error rate across 390 tests requires P < 1.28e-4 (vs the current P<0.01).
Applying that to the real atlas (`../tables/NullCalibration_bonferroni_impact.csv`):

| Trait | Significant at P<0.01 (current) | Survives Bonferroni P<1.28e-4 |
|---|---|---|
| Heart rate | 8 | 1 (catlas_fetal_A_Cardiomyocyte) |
| Red blood cell traits | 6 | 1 (catlas_fetal_Erythroblast_2) |
| Type 1 diabetes | 7 | 2 (catlas_adult CD8+ T, catlas_fetal CD4+ T) |
| Type 2 diabetes | 2 | 0 |

So most of the current "significant" hits for 3 of the 4 traits would not survive formal
multiple-testing correction, and none would for type 2 diabetes. The single strongest
hit per trait (the one already emphasised in the Results text) does survive in 3 of 4
cases.

## What this means for the manuscript — resolved
Decision made (28 Oct 2026): apply Benjamini-Hochberg FDR correction (q<0.10) within
each trait, across the 390 kept cell types, and report that as the significance
criterion throughout — not raw P<0.01, and not Bonferroni (too conservative for this
screen; BH-FDR is the standard choice for this kind of multi-cell-type test and is
what comparable genomics papers use). Applied throughout `crest-gv_nature_methods.tex`
and `crest-gv_supplementary.tex` (the `\qthr` macro, formerly a raw P-value threshold
of 0.01, is now the FDR q threshold, 0.10) and `paper_analysis/figures/Fig3_cross_trait_heatmap.{png,pdf}`
/ `Fig4_per_trait_rankings.{png,pdf}` were regenerated against it (`stats_utils.py`'s
`bh_fdr_mask`). Net effect on the four highlighted traits: heart rate and red blood
cell traits each keep 2 significant cell types, type 1 diabetes keeps its top 2
(unchanged from the raw-threshold version), and type 2 diabetes now has zero
significant hits — reported plainly as a negative result rather than kept at the old
uncorrected count of 2.

Checked one adjacent question while deciding this: CHEERS's own enrichment P-value
(`CHEERS_computeEnrichment.py`) has no multiple-testing correction either — same issue,
never addressed by the original tool. This doesn't affect the CHEERS-vs-CREST-GV
Spearman correlation already reported (`../cheers_benchmark/`), since rank correlation
is invariant to the monotonic rescaling FDR/Bonferroni correction applies.

## Files
- `<trait>_null_results.tsv`: raw Z and (uncorrected, per-collection) P for every
  (replicate, collection, cell type) — 51,100 rows each (100 reps x 511 profiles).
- `../tables/NullCalibration_fpr_by_trait.csv`, `_by_collection.csv`: empirical FPR vs
  nominal alpha.
- `../tables/NullCalibration_atleastone_hit.csv`: P(>=1 false hit across the atlas) per
  trait per alpha.
- `../tables/NullCalibration_bonferroni_impact.csv`, `_threshold.txt`: impact of formal
  correction on the real atlas's reported hits.
- `../figures/SuppFig_null_calibration.{png,pdf}`: QQ-plot + empirical-FPR-vs-nominal bar chart.
- `../data/null_calibration/`: the 400 generated null variant `.bed` files themselves
  (~473MB total, regenerable from the fixed seed in `null_calibration_generate.py` —
  not required once `<trait>_null_results.tsv` exists, safe to delete to save space).
