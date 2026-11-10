# CREST-GV vs CHEERS benchmark — final atlas rerun

Rerun of the CHEERS comparison (`CHEERS/CREST-GV_vs_CHEERS/scripts/compare_with_crestgv.py`)
on the final, QC-filtered, 11-collection atlas. See
`../scripts/cheers_benchmark_compare.py` for the exact script and its docstring for what
changed vs the original comparison.

**Did not re-run CHEERS itself** (consensus peaks / bigwig quantification /
CHEERS_normalize.py / CHEERS_computeEnrichment.py) — those steps are unaffected by the
catlas_fetal QC-pipeline pairing bug and already cover all 11 atlas collections (verified
by a byte-for-byte diff of a sample pair's raw CREST-GV stats against the canonical copy).
Only the comparison/summary step was rerun, with: `pancreatic_pbmc` dropped, restricted to
QC-kept (Good/Decent) cell types, and using the final 390-profile re-standardised
CREST-GV score/p-value instead of the old (unfiltered-collection) standardisation.

## Results

- **44/44 collection x trait pairs** (11 collections x 4 traits).
- **Overall median Spearman rho = 0.644** (vs the manuscript's current draft text, "0.59" —
  that number is from the old, pre-fix comparison and needs updating).
- One pair (`super_pbmc` / `red_blood_cell_traits`) has an undefined correlation: CHEERS's
  own output was degenerate there (p = 0.0 for all 19 cell types), not a CREST-GV issue.
  A few more pairs have `rho_vs_raw` = NaN for the same reason on the *raw* (unstandardised)
  CREST-GV p-value specifically (it underflows to 0 for strong/large-N signals — the
  already-documented CREST-GV pmf/sf/raw-p-value underflow issue, see
  this repo's `CHANGELOG.md`, point 5); the standardised p-value (`rho_vs_norm`, the
  one actually reported) is unaffected in those cases and was used as the primary
  statistic throughout.

**Per-trait median** (highest to lowest):
| Trait | Median rho |
|---|---|
| Type 1 diabetes | 0.715 |
| Red blood cell traits | 0.674 |
| Heart rate | 0.590 |
| Type 2 diabetes | 0.521 |

**Per-collection median** (highest to lowest):
| Collection | Median rho |
|---|---|
| erythoid_d7_d10_d13_d17 | 0.871 |
| cad | 0.761 |
| greenleaf_brain | 0.753 |
| h1_hescs | 0.700 |
| catlas_fetal | 0.667 |
| immune_cell | 0.612 |
| mpal | 0.525 |
| catlas_adult | 0.473 |
| super_pbmc | 0.444 |
| calderon | 0.122 |
| ludwig2019 | 0.024 |

**This reorders the manuscript's current collection-ranking claim.** The current draft
text says concordance is "highest in catlas_adult and catlas_fetal (median rho = 0.79 and
0.75)" — on the final atlas, catlas_adult is now one of the *lower*-concordance
collections (0.47, 8th of 11), while catlas_fetal is mid-pack (0.67). The genuinely
highest-concordance collections are erythoid_d7_d10_d13_d17, cad and greenleaf_brain.
calderon and ludwig2019 remain the lowest, consistent with the existing text's explanation
(small collections / calderon's donor-replicate structure).

## Files

- `summary.tsv` — one row per (collection, trait) pair: n cell types, Spearman rho (vs raw
  and vs standardised CREST-GV p-value), significance counts.
- `medians.tsv` — overall / per-trait / per-collection medians.
- `per_pair/*.tsv` — full per-cell-type merged table for each of the 44 pairs.
- `../figures/SuppFig_CHEERS_scatter.{png,pdf}` — (a) per-collection concordance, (b) pooled
  -log10(P) scatter across all pairs, coloured by trait.
