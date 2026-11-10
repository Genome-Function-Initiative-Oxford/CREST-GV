#!/usr/bin/env python3

"""
:File: config.py

Shared configuration for the CREST-GV_vs_CHEERS pipeline: runs the CHEERS
cell-type-specificity enrichment method on the CREST-GV bigWig/peak
collections, for comparison against CREST-GV's own enrichment statistics
(the *_ps.csv files already present in each collection directory).

Peak-level signal is quantified directly from RPKM-normalized bigWigs
(bigWigAverageOverBed "sum" column over each consensus peak) rather than
raw read counts. This is intentional, not a shortcut: CHEERS_normalize.py's
own pipeline (match-to-max library size -> quantile norm -> Euclidean norm)
is invariant to any pure per-sample multiplicative constant, which is
exactly what RPKM normalization at a fixed bin size applies. See the
project discussion for the derivation. No RPKM reversal is performed.
"""

import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPELINE_ROOT = os.path.dirname(os.path.abspath(__file__))

CHEERS_NORMALIZE = os.path.join(REPO_ROOT, 'CHEERS_normalize.py')
CHEERS_ENRICHMENT = os.path.join(REPO_ROOT, 'CHEERS_computeEnrichment.py')

COLLECTIONS_BASE = '/project/Wellcome_Discovery/datashare/CREST-GV_collection'

# Internal key -> on-disk collection directory name
COLLECTIONS = {
    'cad': 'CAD',
    'calderon': 'calderon',
    'catlas_fetal': 'catlas_fetal',
    'catlas_adult': 'catlas_adult',
    'erythoid_d7_d10_d13_d17': 'Days7_10_13_17',
    'h1_hescs': 'H1_hESCs',
    'immune_cell': 'immune_cell',
    'ludwig2019': 'ludwig2019',
    'mpal': 'MPAL',
    'super_pbmc': 'super_PBMC',
    'greenleaf_brain': 'greenleaf_brain',
}

TRAITS = {
    'type_1_diabetes_mellitus': '/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_1_diabetes_mellitus.bed',
    'type_2_diabetes_mellitus': '/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/type_2_diabetes_mellitus.bed',
    'red_blood_cell_traits': '/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/red_blood_cell_traits.bed',
    'heart_rate': '/project/Wellcome_Discovery/shared/imputed_genetics_per_trait_collection/heart_rate.bed',
}

DATA_DIR = os.path.join(PIPELINE_ROOT, 'data')
CONSENSUS_PEAKS_DIR = os.path.join(DATA_DIR, 'consensus_peaks')
COUNTS_DIR = os.path.join(DATA_DIR, 'counts')

RESULTS_DIR = os.path.join(PIPELINE_ROOT, 'results')
NORMALIZE_RESULTS_DIR = os.path.join(RESULTS_DIR, 'normalize')
ENRICHMENT_RESULTS_DIR = os.path.join(RESULTS_DIR, 'enrichment')

LOG_DIR = os.path.join(PIPELINE_ROOT, 'logs')

COMPARISON_RESULTS_DIR = os.path.join(RESULTS_DIR, 'comparison')

# bigWigAverageOverBed (ucsctools/385) has a reproducible bug: for at least
# one specific peak/bigWig combination (chr1:248944146-248946588 in several
# greenleaf_brain cell types) it returns a garbage "sum" value (~1e300,
# occasionally negative) instead of erroring out. Signal sums can never be
# negative, and no real RPKM-over-peak sum approaches this ceiling (the
# largest legitimate values observed across all 12 collections are ~1e7);
# anything past it is corrupted tool output, not biology. A single such
# value otherwise silently contaminates an entire collection's results,
# since normalize_counts()'s per-sample scaling and quantile_norm()'s shared
# rank-value table both pool across peaks/samples, and euclidean_norm()
# squares values - overflowing to inf and collapsing the whole matrix to 0.
SANE_SUM_CEILING = 1e9

# CREST-GV's own per-trait enrichment results (distinct from the *_ps.csv
# background/calibration files under COLLECTIONS_BASE). Layout:
# <CRESTGV_RESULTS_BASE>/<trait_name>/<collection_key>/<collection_key>_statistics_CREST-GV.csv
# trait_name and collection_key match the keys of TRAITS and COLLECTIONS above.
CRESTGV_RESULTS_BASE = '/project/Wellcome_Discovery/sriva/Git/CREST-GV_internal/CREST_HYPER_Paper'

# External tools resolved by path rather than relying on `module load`
# inside non-interactive subprocesses (Lmod's `module` is a shell function,
# not reliably available to subprocess.run). Update if module versions change.
BEDTOOLS_BIN = '/package/bedtools/2.29.2/bin/bedtools'
BIGWIG_AVERAGE_OVER_BED_BIN = '/package/ucsctools/385/bin/bigWigAverageOverBed'

# The system default `python3` has a broken NumPy install (missing
# libblas.so.3). `module load python-cbrg/current` provides a working
# NumPy/pandas/SciPy stack; resolved here to its absolute binary path so it
# can be used from subprocess calls without depending on `module load`.
# The interpreter also needs LD_LIBRARY_PATH (to find libpython3.11.so.1.0)
# and PYTHONPATH (to find numpy/pandas/scipy, installed under python-cbrg's
# own site-packages rather than python-base's) - both normally set up by
# Lmod as part of `module load`, diffed here from `env` before/after a real
# `module load python-cbrg/current` and replicated explicitly.
PYTHON_BIN = '/package/python-base/3.11.15/bin/python3'
PYTHON_ENV = dict(os.environ)
PYTHON_ENV['LD_LIBRARY_PATH'] = (
    '/package/python-base/3.11.15/lib:'
    '/package/python-base/3.11.15/lib/python3.11/config-3.11-x86_64-linux-gnu/'
    + (':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else '')
)
# NumPy itself happens to live only under the 3.11.14 sub-tree (not 3.11.15,
# despite PYTHON_BIN being the 3.11.15 interpreter) - Lmod stacks every
# version's site-packages onto PYTHONPATH, so replicate that stack exactly
# rather than guessing which one path has which package.
PYTHON_ENV['PYTHONPATH'] = ':'.join(
    f'/package/python-cbrg/current/{v}/lib/python3.11/site-packages'
    for v in ('3.11.15', '3.11.14', '3.11.13', '3.11.11', '3.11.3')
)
