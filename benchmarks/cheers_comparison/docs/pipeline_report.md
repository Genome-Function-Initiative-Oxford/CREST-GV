---
title: "CREST-GV vs CHEERS: Pipeline Design and RPKM Justification"
author: "CHEERS repository — CREST-GV_vs_CHEERS pipeline"
date: "2026-11-10"
geometry: margin=2.5cm
---

# 1. Purpose

CREST-GV produces, for each of 11 reference chromatin accessibility collections, a
per-cell-type enrichment statistic (`*_ps.csv`: `n_snps`, `p_success`) describing how
strongly a set of trait-associated SNPs falls into that cell type's accessible
chromatin. This pipeline runs the independent, previously published **CHEERS**
(Chromatin Element Enrichment Ranking by Specificity) method on the same underlying
data — RPKM-normalized bigWig tracks and per-cell-type peak calls — so that CHEERS's
enrichment p-values can be compared against CREST-GV's own statistic, collection by
collection and trait by trait.

CHEERS was not designed to consume bigWig tracks; it expects raw integer read counts
per peak. Section 3 explains why RPKM-normalized bigWig signal can be used directly,
without reconstructing raw counts, and why that is not a compromise but a
mathematically exact substitution given CHEERS's own normalization steps.

# 2. Input data

## 2.1 CREST-GV collections

Location: `/project/Wellcome_Discovery/datashare/CREST-GV_collection/`

Eleven collections are used, each with an identical internal layout. (An earlier
version of this comparison also included `pancreatic_pbmc`; it was subsequently
dropped from the released atlas for poor data quality and is excluded here too.)

| Internal key | On-disk directory |
|---|---|
| `cad` | `CAD` |
| `calderon` | `calderon` |
| `catlas_fetal` | `catlas_fetal` |
| `catlas_adult` | `catlas_adult` |
| `erythoid_d7_d10_d13_d17` | `Days7_10_13_17` |
| `h1_hescs` | `H1_hESCs` |
| `immune_cell` | `immune_cell` |
| `ludwig2019` | `ludwig2019` |
| `mpal` | `MPAL` |
| `super_pbmc` | `super_PBMC` |
| `greenleaf_brain` | `greenleaf_brain` |

Within each collection directory:

```
<collection>/
  bigwigs/<cell_type>.bw        # RPKM-normalized, bamCoverage --binSize 1
  peaks/<cell_type>_L-tron.bed  # per-cell-type peak calls (chr, start, end)
  <collection>_info.csv         # cell counts per type
  <collection>_qc.csv           # QC labels per type
  <collection>_ps.csv           # CREST-GV's own enrichment statistic (n_snps, p_success)
```

Peak calls are made **independently per cell type** — there is no shared consensus
peak set across cell types within a collection, which CHEERS requires (its
normalization inner-joins all samples on identical `chr/start/end` coordinates).

The comparison in Section 5 is restricted to the 390 QC-kept (Good/Decent) cell
types across these 11 collections, not the full 511-profile atlas — the same
restriction CREST-GV's own precomputed resource applies.

## 2.2 Trait SNP lists

Four GWAS traits, already in the exact format CHEERS expects
(`name<TAB>chr<TAB>position`, no header):

| Trait | Path |
|---|---|
| `type_1_diabetes_mellitus` | `/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_1_diabetes_mellitus.bed` |
| `type_2_diabetes_mellitus` | `/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_2_diabetes_mellitus.bed` |
| `red_blood_cell_traits` | `/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/red_blood_cell_traits.bed` |
| `heart_rate` | `/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/heart_rate.bed` |

No reformatting is required for these files.

# 3. Why RPKM bigWig signal can be used directly (no reversal needed)

## 3.1 The problem

The bigWigs were generated with `bamCoverage --normalizeUsing RPKM --binSize 1`. The
original BAM files and per-sample mapped-read counts (flagstat logs, QC reports) are
not available. RPKM normalization applies, per sample, a single scale factor:

$$ \text{scale}_i = \frac{10^9}{\text{mapped\_reads}_i} \quad (\text{at binSize} = 1\,\text{bp}) $$

so that at every genomic position, `RPKM_value = raw_coverage * scale_i`.

Because RPKM is defined specifically to divide out each sample's total mapped-read
count, the genome-wide sum of an RPKM track is **invariant** to that count (it
depends only on mean fragment length). This means `mapped_reads_i` cannot be
recovered from the bigWig alone — true reversal of RPKM is not merely
inconvenient without flagstat data, it is mathematically impossible.

## 3.2 Why this does not matter for CHEERS

`CHEERS_normalize.py` applies four steps, in order, to a peak x sample count matrix.
Consider feeding it RPKM-summed peak values ($y_i = \text{scale}_i \cdot x_i$, where
$x_i$ is the vector of true raw per-peak counts for sample $i$) instead of true raw
counts ($x_i$):

1. **`normalize_counts()`** rescales every sample to match the sample with the
   largest total peak signal:
   $$\text{output}_i = x_i \cdot \text{norm\_factor}_i, \quad
     \text{norm\_factor}_i = \frac{\max_j(\text{sum}(y_j))}{\text{sum}(y_i)}$$
   Substituting $y_i = \text{scale}_i \cdot x_i$ and simplifying shows that the
   RPKM-based output differs from the raw-count-based output by exactly **one
   global constant**, identical across every sample and every peak — not a
   per-sample bias. (Full algebra: letting $S_i = \text{sum}(x_i)$, the RPKM-based
   result equals $x_i \cdot \max_j(\text{scale}_j S_j) / S_i$, versus
   $x_i \cdot \max_j(S_j)/S_i$ from raw counts — the two differ only by the ratio
   $\max_j(\text{scale}_j S_j) / \max_j(S_j)$, a single number independent of $i$.)

2. **`drop_peaks()`** (bottom-10th-percentile filter by total signal) is invariant
   to any positive uniform rescaling — percentile cutoffs don't move when every
   value is scaled by the same constant. The same peaks are retained either way.

3. **`quantile_norm()`** operates on within-column ranks. Multiplying every entry
   by the same positive constant does not change any rank, so the rank-mean
   lookup table is scaled by that same constant and the output carries it through
   unchanged in relative terms.

4. **`euclidean_norm()`** divides each peak's row vector by its own Euclidean norm
   across samples. This step **exactly cancels** any remaining global constant:
   $$\frac{k \cdot z}{\lVert k \cdot z \rVert} = \frac{z}{\lVert z \rVert}$$

**Conclusion:** the final Euclidean-normalized specificity scores are
mathematically identical whether CHEERS is given true raw per-peak counts or
RPKM-normalized bigWig signal summed over the same peaks — provided the RPKM
scale factor is a pure per-sample constant applied uniformly across every
position, which is exactly what `bamCoverage --normalizeUsing RPKM --binSize 1`
produces.

## 3.3 Conditions and minor caveats

- This holds exactly only if every bigWig in a collection was generated with
  **consistent `bamCoverage` settings** (same `--binSize`, no sample-specific
  `--ignoreForNormalization`, consistent fragment-extension handling). Mixed
  settings would reintroduce a position- or sample-dependent bias that the
  derivation above does not cover.
- Raw integer counts can produce exact ties between peaks; RPKM floating-point
  values essentially never tie exactly. `rank(method='first')`'s tie-breaking
  order could therefore differ very slightly between the two inputs — a
  negligible practical effect, not a structural one.
- Output values are no longer interpretable as literal "read counts", only as
  proportional specificity scores — which is all CHEERS consumes downstream.

**Net effect: no RPKM reversal is performed anywhere in this pipeline.** Peak
signal is quantified directly from the RPKM bigWigs.

# 4. Pipeline implementation

All code lives under `CREST-GV_vs_CHEERS/` in the CHEERS repository.

```
CREST-GV_vs_CHEERS/
  config.py                        # paths, collection map, trait map, tool paths
  scripts/
    build_consensus_peaks.py       # step A, per collection
    quantify_bigwigs.py            # step B, per collection
    run_collection.py              # orchestrates A -> D for one collection
  data/
    consensus_peaks/<collection>.bed          # merged, sorted, non-overlapping
    consensus_peaks/<collection>_named.bed    # + unique 4th-column peak ID
    counts/<collection>/<cell_type>_ReadsInPeaks.txt
  results/
    normalize/<collection>/...
    enrichment/<collection>/<trait>/...
  docs/
    pipeline_report.md / .pdf      # this document
```

## Step A — Build consensus peaks (`build_consensus_peaks.py`)

Concatenates every cell type's peak BED file within a collection, sorts
(`bedtools sort`) and merges overlapping intervals (`bedtools merge`) into one
non-overlapping consensus peak set. A unique peak ID (`chrom_start_end`) is added
as a 4th column, required by `bigWigAverageOverBed`.

## Step B — Quantify bigWig signal over peaks (`quantify_bigwigs.py`)

For every cell type's bigWig in a collection, runs:

```
bigWigAverageOverBed <bigwig> <consensus_named.bed> <tmp.tab>
```

and takes the **`sum`** column (total signal over all bases in the peak — the
RPKM-domain equivalent of a raw per-peak read count, per Section 3). Output is
written as `<cell_type>_ReadsInPeaks.txt`: tab-separated `chr, start, end, value`,
no header — exactly the format `CHEERS_normalize.py` expects.

## Step C — Normalize (`run_collection.py` -> `CHEERS_normalize.py`)

Runs the repository's existing, unmodified `CHEERS_normalize.py` on all of a
collection's count files, producing library-size-normalized, quantile-normalized,
and Euclidean-normalized (cell-type specificity score) matrices.

## Step D — Enrichment (`run_collection.py` -> `CHEERS_computeEnrichment.py`)

Runs the repository's existing, unmodified `CHEERS_computeEnrichment.py` once per
trait against the collection's specificity scores, producing per-cell-type
enrichment p-values, observed mean ranks, and a log file.

## Tooling notes specific to this cluster

- `awk`, GNU `parallel` are not installed; all text processing is done in Python.
- `bedtools` and `bigWigAverageOverBed` (`ucsctools`) are available only via
  `module load`; since Lmod's `module` command is a shell function rather than an
  executable, `config.py` instead hardcodes resolved absolute binary paths
  (`/package/bedtools/2.29.2/bin/bedtools`,
  `/package/ucsctools/385/bin/bigWigAverageOverBed`) so subprocess calls work
  without depending on an interactive login shell.
- `sbatch`/`srun` are available (SLURM cluster); batch submission can be added if
  running all 11 collections sequentially proves too slow.

# 5. Results

The pipeline was run to completion on all 11 released atlas collections across
all 4 traits (44 collection/trait pairs, restricted to the 390 QC-kept cell
types). All 44 pairs produced normalized specificity scores and enrichment
p-values; none failed.

## 5.1 A real bug found and fixed during QA: corrupted `bigWigAverageOverBed` output

Before trusting the results, every collection's raw quantified counts were
scanned for implausible values (negative, or exceeding a generous sanity
ceiling). This caught a genuine, reproducible bug: for `greenleaf_brain`,
`bigWigAverageOverBed` (ucsctools v385) returned a garbage "sum" value
(`1.8783e+300`, or `-5.8e+305` in one file) for one specific peak
(`chr1:248944146-248946588`) in 5 of its 20 bigWigs, instead of erroring.

That single corrupted value was enough to silently invalidate the **entire**
collection's result: `normalize_counts()`'s library-size scaling and
`quantile_norm()`'s shared rank-value table both pool across every peak and
sample, so the contamination spread to all 20 cell types; `euclidean_norm()`
then squared those now-astronomical values, overflowing IEEE-754 double
precision to `inf`, and `finite / inf = 0.0` silently collapsed the entire
257,279-peak x 20-cell-type specificity matrix to zero - no exception was ever
raised. Every peak tied for rank 1, and every enrichment p-value for
`greenleaf_brain` was degenerate as a result, for all 4 traits.

**Fix:** `quantify_bigwigs.py` now rejects any `bigWigAverageOverBed` output
that is negative, non-finite, or exceeds `config.SANE_SUM_CEILING` (`1e9`, far
above the ~1e7 maximum legitimate value observed across all 11 collections),
replacing it with `0.0` and logging a warning. All other 10 collections were
re-scanned for the same signature and found clean; only `greenleaf_brain` was
affected. After the fix, `greenleaf_brain`'s specificity matrix shows 0%
degenerate rows and enrichment p-values spread sensibly around the theoretical
null, as expected.

## 5.2 Overall concordance with CREST-GV

CHEERS's output was compared against CREST-GV's own `CREST-GV_pval_norm`
column, not the raw `CREST-GV_pval`, via Spearman rank correlation, per
collection/trait pair. CREST-GV computes `CREST-GV_pval_norm` as
`norm.sf(Z_norm)`, where `Z_norm = (Z - mean(Z)) / std(Z)` is re-standardized
across all cell types *within one collection* (`crestgv.py`,
`__add_statistics()`). This re-standardization is what makes the statistic
comparable across collections of very different size and composition - the
raw `CREST-GV_pval` is calibrated only against a trait-independent,
genome-wide background rate and does not correct for collection-specific
scale differences, so `CREST-GV_pval_norm` is the column its own
implementation intends for this kind of comparison.

| | Value |
|---|---|
| Pairs compared | 44 (11 collections x 4 traits) |
| Median Spearman &rho; | 0.644 |
| Mean Spearman &rho; | 0.570 |

By trait (median &rho;):

| Trait | Median &rho; | Range |
|---|---|---|
| Type 1 diabetes | 0.715 | 0.400 - 0.933 |
| Red blood cell traits | 0.674 | -0.143 - 1.000 |
| Heart rate | 0.590 | -0.190 - 0.908 |
| Type 2 diabetes | 0.521 | 0.036 - 1.000 |

By collection (median &rho; across its 4 traits):

| Collection | Median &rho; | | Collection | Median &rho; |
|---|---|---|---|---|
| erythoid_d7_d10_d13_d17 | 0.871 | | catlas_adult | 0.473 |
| cad | 0.761 | | super_pbmc | 0.444 |
| greenleaf_brain | 0.753 | | calderon | 0.122 |
| h1_hescs | 0.700 | | ludwig2019 | 0.024 |
| catlas_fetal | 0.667 | | | |
| immune_cell | 0.612 | | | |
| mpal | 0.525 | | | |

Note the ranking changed substantially from an earlier version of this comparison
(run before a chromatin-track QC pairing bug in `catlas_fetal` was fixed, and before
`pancreatic_pbmc` was dropped from the atlas): `catlas_adult` was previously the
highest-concordance collection (0.785) and is now one of the lower ones (0.473);
`erythoid_d7_d10_d13_d17`, previously mid-pack (0.420), is now the highest (0.871).

A caveat on significance counts: counting cell types "significant" is lopsided and
should **not** be read as a concordance metric. Using the same Benjamini-Hochberg
FDR q<0.10 criterion the CREST-GV paper applies throughout (computed within each
trait, pooled across all 390 kept cell types), CHEERS flags 1,073 of 1,560
cell-type/trait combinations significant, CREST-GV flags only 6 (all 6 contained
within CHEERS's set). This is not CHEERS finding more real signal; CHEERS's raw
p-value has no cross-collection calibration and underflows toward 0 very easily in
large, strongly-enriched collections (see Section 5.5) - for red blood cell traits
specifically, 88% of CHEERS's raw p-values already fall below 0.01, with the median
underflowing to exactly 0. **Rank correlation (&rho;), not the significance count,
is the concordance metric used throughout this report.**

## 5.3 Biological validation: top enriched cell types per trait

As an independent face-validity check (not a comparison against CREST-GV),
CHEERS's own top hits per trait, ranked by Z-score across all 11 collections,
recover the expected biology in every case:

| Trait | Top cell types (by Z-score) |
|---|---|
| Type 1 diabetes | T lymphocytes (CD4+, CD8+), thymocytes, cytotoxic T, naive T, macrophages - consistent with T1D's T-cell-mediated autoimmune pathology |
| Type 2 diabetes | Pancreatic acinar, hepatoblast, skeletal myocytes (type I/II), adipocyte - metabolic/insulin-resistance tissues, not immune cells |
| Red blood cell traits | Erythroblasts (stages 1-4), megakaryocyte - directly the erythroid lineage (Z up to 73, the strongest signal of all 4 traits) |
| Heart rate | Atrial and ventricular cardiomyocytes, vascular smooth muscle, endothelial - cardiac tissue, as expected |

## 5.4 A statistical caveat for small collections

Three collections have few enough cell types that their individual
collection/trait Spearman &rho; should not be over-interpreted:
`h1_hescs` (n=4), `ludwig2019` (n=8), `CAD` (n=10).
Checking statistical significance of the correlation itself (not just its
point estimate) makes this concrete: **none of `ludwig2019`'s four
correlations are statistically significant** (p=0.207, 0.651, 0.736, 0.651)
- with only 8 cell types, a Spearman &rho; anywhere in a wide range around
zero is statistically indistinguishable from noise.

This also surfaces a numerical quirk worth documenting: `h1_hescs` (n=4) shows
two *perfect* correlations (&rho;=1.0 for both type 2 diabetes and red blood
cell traits). SciPy's `spearmanr` reports `pvalue=0.0` for these - but that is
a known breakdown of its t-distribution-based approximation at `|&rho;|=1`
(division by zero in the test statistic), **not** the correct value. Brute-force
enumeration of all 4! = 24 possible rank permutations for n=4 gives the true
exact two-sided p-value for a perfect correlation as 2/24 &asymp; 0.083 - not
significant at the conventional 0.05 threshold, despite being the best
possible outcome achievable with only 4 data points. **Conclusion: no
individual &rho; value from a collection with n&lt;12 or so should be read as
evidence of agreement or disagreement between the two methods on its own.**
The larger, better-powered collections (catlas_adult/fetal, calderon,
immune_cell, super_pbmc, mpal, greenleaf_brain, Days7_10_13_17; n=15-175) carry
the real statistical weight in Section 5.2's summary.

## 5.5 Case study: `calderon`'s near-zero correlation (&rho;=0.122 median)

Unlike the small-collection cases above, `calderon` (n=175, among the largest
collections) has a robust, well-powered near-zero correlation for 2 of its 4
traits that is not a statistical-power artifact, and is worth explaining
directly since it is the weakest well-powered collection in this comparison.

**Structure:** `calderon`'s 175 bigWigs are not 175 distinct biological cell
types - they are **45 distinct cell-type/stimulation conditions, each
replicated across ~4 donors** (e.g. `1001-Th1_precursors-S`,
`1002-Th1_precursors-S`, `1004-...`, `1008-...`). This is the Soskic et al.
(2019) T-cell dataset CHEERS was originally developed and validated on (the
same 45 condition names ship as `test_input/` elsewhere in this repository,
donor-pooled); CREST-GV's collection instead keeps every donor's bigWig and
peak calls unaggregated.

**Both methods are substantially donor-noise-dominated in this collection.**
Decomposing each method's own Z-score variance (across calderon's 175 points)
into within-base-cell-type (donor noise) versus between-base-cell-type (true
biological signal) gives:

| Trait | CHEERS donor-noise variance | CREST-GV donor-noise variance | Observed &rho; |
|---|---|---|---|
| Type 1 diabetes | 13.2% | 39.4% | 0.673 |
| Heart rate | 39.1% | 49.9% | 0.206 |
| Type 2 diabetes | 43.4% | 57.4% | 0.036 |
| Red blood cell traits | 45.1% | 61.7% | 0.037 |

Since correlating two independently noisy measurements is mathematically
attenuated by the product of their respective noise fractions (an upper bound
of roughly &radic;((1-w<sub>CHEERS</sub>)(1-w<sub>CREST-GV</sub>)) on the
observable correlation), donor noise alone predicts a correlation ceiling of
&asymp;0.73 for type 1 diabetes (observed: 0.673 - right at the ceiling) but
only &asymp;0.49 and &asymp;0.46 for type 2 diabetes and red blood cell traits
respectively - yet their *observed* &rho; (0.036, 0.037) sits far *below* even
that reduced ceiling. Donor noise alone does not fully explain those two.

**The remaining gap has a biological explanation.** `calderon` is a T-cell/PBMC
collection. Type 1 diabetes is a T-cell-mediated autoimmune disease - a direct
biological match, so both methods detect real, concordant signal despite the
noise. Type 2 diabetes (pancreas/liver/muscle/fat) and red blood cell traits
(erythroid lineage) are not primarily mediated through T cells at all - there
is effectively **no true enrichment signal for either method to find** in a
T-cell-only collection for those two traits. Under a near-null signal, two
different statistical tests' rankings are expected to show close to zero
correlation purely by chance, which is exactly what is observed.

**Net conclusion:** `calderon`'s low correlation is not a flaw in either tool.
It reflects (1) donor-replicate noise diluting an otherwise-real signal for
type 1 diabetes, compounded by (2) a genuine absence of true signal, for
either method, on traits this T-cell collection was never going to be
informative for. (A direct confirmatory test - re-running CHEERS on
donor-pooled calderon counts, as the method was originally validated - was
considered but not run; it was judged unnecessary to support the conclusion
above, which is already supported by the variance decomposition and the
trait-specific biological pattern.)

## 5.6 Interactive visualization

A per-collection, per-trait scatter-plot comparison (CHEERS vs CREST-GV,
toggleable between Z-score and raw p-value views, colored by significance
status in each tool) is available at
`CREST-GV_vs_CHEERS/docs/cheers_vs_crestgv_scatter.html` (self-contained,
opens in any browser, no server required). Each subplot's axes are scaled
independently (CHEERS's Z-score is unbounded and can reach into the tens or
hundreds for very strong signals - e.g. erythroblasts for red blood cell
traits, Z&asymp;73 - while CREST-GV's `_norm` statistic is standardized to
roughly [-3, 3] by construction), so panels should be read for rank/pattern,
not absolute axis position.

# 6. Conceptual differences between the two methods

Reading CREST-GV's own implementation (`crestgv/crestgv.py`) clarifies that
CHEERS and CREST-GV are not two statistics answering the same question - this
matters for interpreting both the agreement and the disagreement above.

- **CHEERS** asks whether a GWAS SNP sits in chromatin that is
  *specifically/disproportionately* accessible in one cell type relative to
  all others in the collection. It requires one shared consensus peak set
  across all cell types (built by this pipeline's `build_consensus_peaks.py`,
  since CREST-GV's collections provide independently-called peaks per cell
  type) and ranks peaks by a continuous, Euclidean-normalized specificity
  score.
- **CREST-GV** asks whether a cell type's *overall* accessible genome is
  enriched for trait SNPs relative to a genome-wide background rate (the
  fraction of 1000 Genomes common variants falling in that cell type's own
  peaks). It is a binary peak-overlap proportion (Z-)test, computed
  independently per cell type against that cell type's own native peak calls
  - it does not use bigWig signal strength in the statistic at all, and
  carries no relative-specificity weighting between cell types.

Two consequences follow. First, a locus accessible in every cell type
contributes equally to every cell type's CREST-GV statistic, but is
down-weighted by CHEERS's specificity ranking - so some disagreement between
the methods reflects this genuine conceptual difference, not error in either.
Second, because CHEERS requires a shared peak space, this pipeline's
consensus-merge step means CHEERS and CREST-GV are not even evaluated against
literally identical peak coordinates within a "matched" collection - a further
expected source of difference, independent of the statistical model.

A review of CREST-GV's source code surfaced several hardening opportunities
(a stale CLI/class parameter mismatch, division-by-zero risks in `SE` and
`Z_norm` analogous to the `bigWigAverageOverBed` failure mode found in
Section 5.1, and the exact-binomial `pmf`/`sf` columns underflowing to exactly
0.0 for strong signals). These have since been fixed in the refactored
`crestgv` package released alongside the paper (see this repository's own
`CHANGELOG.md`) and are orthogonal to the comparison results above, which use
CREST-GV's already-computed, supported output columns (`CREST-GV_pval`,
`CREST-GV_pval_norm`).

# 7. Conclusion

CHEERS and CREST-GV show positive, generally moderate-to-strong rank
concordance (median &rho;=0.644 across 44 well-defined collection/trait pairs)
on cell types they both identify as most relevant to a trait, and both
independently recover biologically expected top hits per trait (Section 5.3).
Where they diverge strongly, the divergence is explainable: either the
collection is too small to measure concordance reliably at all (Section 5.4),
or - as in `calderon` - a combination of unaggregated biological-replicate
noise and genuine absence of trait-relevant signal in that specific reference
collection (Section 5.5). Neither pattern indicates a defect in either method;
together with the conceptual differences in Section 6, they explain the full
range of agreement and disagreement observed.
