# Changelog

## 0.3.0 — refactor for GitHub/PyPI release

Rewrite of `crestgv/crestgv.py` + `crestgv.py` into an installable package
(`src/crestgv/core.py` + `src/crestgv/cli.py`). Fixes every item in the previous
`KNOWN_ISSUES.md` (numbers below refer to that file's numbering) plus two further issues
found while writing the test suite for this release. No change to the underlying
enrichment statistics for well-formed input — see "Behavioural changes" below for the
cases that do change.

### Fixed

1. **CLI/class parameter mismatch (KNOWN_ISSUES #1).** The CLI no longer passes
   `number_of_folds`, `lessNG` or `greater25k` — none of these exist on the class anymore,
   and calling the CLI as documented previously crashed immediately with
   `TypeError: __init__() got an unexpected keyword argument 'number_of_folds'`.
2. **Empty `_version.py` (KNOWN_ISSUES #2).** Populated (`0.3.0`); the package now imports
   cleanly.
3. **Division-by-zero in `SE`/`Z` when `p_succes` is 0 or 1 (KNOWN_ISSUES #3).** Cell
   types with a degenerate background rate are now excluded from the enrichment columns
   (`SE`/`Z`/`pdf` set to `NaN`) with an explicit `warnings.warn`, instead of silently
   producing `inf`/`nan`.
4. **Division-by-zero in `Z_norm` for small/flat collections (KNOWN_ISSUES #4).** If a
   collection's standard deviation of `Z` is zero or undefined (e.g. `h1_hescs`, 4 cell
   types, if they happen to score near-identically), `Z_norm` is now set to `0.0` for
   every cell type in that collection, with a `warnings.warn`, instead of silently
   producing `inf`.
5. **`pmf`/`sf` underflow documented (KNOWN_ISSUES #5).** These exact-binomial columns are
   still written (useful as a diagnostic) but the docstring and README now say explicitly
   that `CREST-GV_pval`/`CREST-GV_pval_norm` are the supported statistics.
6. **Dead code removed (KNOWN_ISSUES #6).** All commented-out code (the old pybedtools
   peak-intersection block, the old `__check_genetic_format`, `__clean_tmp`) is gone, as
   are the unused `self.tmp`/`self.mappable_bp` attributes they supported.
7. **Peak-set semantics documented (KNOWN_ISSUES #7).** Explained in the module docstring
   and the README: each cell type is scored against its own independently-called peaks,
   not a consensus peak set.
8. **`min_number_genetics` CLI/class default mismatch (KNOWN_ISSUES #8).** Both now
   default to `500`.
9. **`--run_all` silently skipped `greenleaf_brain`.** The CLI used to keep its own copy
   of the collection name → directory mapping, which was missing `greenleaf_brain` (the
   class's copy had it). There is now a single source of truth,
   `crestgv.core.COLLECTION_NAME_DICT`, imported by both the class and the CLI.
10. **Boolean CLI flags were unsettable.** `-gc/--get_coverage` (now `-wc/--with_coverage`)
    and `-ra/--run_all` used `type=bool`, a classic argparse footgun:
    `bool("False") == True` in Python, so passing e.g. `--get_coverage False` could never
    actually disable the flag. Both are now `action="store_true"`.
11. **Error message with an unfilled `%s`.** The "too few variants" error in
    `__load_genetic` literally printed the string `%s` instead of the configured
    threshold; it now interpolates `self.min_number_genetics`.
12. **Bundled example bigwigs were 0 bytes.** `example_files/in-house-data/bigwigs/*.bw`
    were empty placeholder files, so `get_coverage()` could never have been run
    successfully against the bundled example — `pyBigWig.open()` fails to open an empty
    file. Replaced with small real bigwig files (synthetic signal, same coordinate ranges
    as the example peaks and the example genetic variants).
13. **Bundled example file didn't match the documented/actual column order.**
    `example_files/genetic.vcf` (misleadingly named — it was never a VCF) had a header row
    and columns in `CHR_ID, CHR_POS, SNPS` order, matching the old README's documented
    format — but the parser (`pd.read_csv(..., names=["SNPS","CHR_ID","CHR_POS"])`) reads
    positionally with no header, expecting `SNPS, CHR_ID, CHR_POS`. Every production GWAS
    trait file (`/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/*.bed`)
    already uses the real `SNPS, CHR_ID, CHR_POS`, no-header format, so the parser itself
    was correct and battle-tested — it was the example file and the README's documented
    format that were wrong, and would have made every new user's first run fail with "no
    valid data remains after filtering". Replaced with
    `example_files/genetic.tsv` in the correct format, and corrected the documented format
    in the README and in `crestgv`'s docstring. Caught by the new test suite, which
    previously could not have existed since the CLI did not run at all (#1).

### Removed

- Unused imports that pulled in heavy, unnecessary dependencies: `matplotlib`, `seaborn`,
  `plotly` (none of their names were referenced anywhere in the class), and a duplicated,
  likewise-unused `import multiprocessing as mp`. `glob`, `subprocess`, `shutil`, `re` were
  also unused outside the removed dead code. This substantially lightens the package's
  runtime dependencies — none of these are needed to compute enrichment statistics.
- `pancreatic_pbmc` from `COLLECTION_NAME_DICT` — dropped from the published reference
  atlas for poor data quality; it no longer makes sense as an accepted `collection_name`.

### Behavioural changes

- A degenerate `p_succes` (0 or 1) or a zero-variance collection no longer silently
  produces `inf`/`nan` in the output — affected cell types get `NaN` (point 3) or `Z_norm
  = 0` (point 4) plus a `warnings.warn`. If your downstream code depended on the previous
  `inf`/`nan` values, you will need to update it.
- Coverage calculation (`get_coverage()`) is no longer run by default when using the CLI;
  pass `-wc/--with_coverage` to also compute it. (Via the Python class, nothing changes:
  call `get_coverage()` explicitly, as before.)
- `collection_name="pancreatic_pbmc"` is no longer accepted (see "Removed").
