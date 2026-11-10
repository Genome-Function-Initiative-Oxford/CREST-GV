#!/usr/bin/env python3

"""
:File: run_collection.py

Runs the full CREST-GV_vs_CHEERS pipeline for one collection:

1. Build a consensus peak set from all cell types' peak calls.
2. Quantify every cell type's bigWig signal over the consensus peaks.
3. Run CHEERS_normalize.py to produce cell-type specificity scores.
4. Run CHEERS_computeEnrichment.py for every configured trait.

Usage:
    python3 run_collection.py <collection_key>

Where <collection_key> is a key of config.COLLECTIONS (e.g. "cad").
"""

import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

from quantify_bigwigs import quantify_collection


def run_normalize(collection_name, counts_dir):

    '''
    Run CHEERS_normalize.py on a collection's quantified count files.

    :param collection_name: The on-disk collection directory name.
    :param counts_dir: Directory containing per-cell-type CHEERS-format count files.
    :return: Path to the Euclidean-normalized specificity score file.
    '''

    outdir = os.path.join(config.NORMALIZE_RESULTS_DIR, collection_name)
    os.makedirs(outdir, exist_ok=True)

    count_files = sorted(glob.glob(os.path.join(counts_dir, '*_ReadsInPeaks.txt')))

    subprocess.run(
        [config.PYTHON_BIN, config.CHEERS_NORMALIZE, collection_name, outdir] + count_files,
        check=True, env=config.PYTHON_ENV,
    )

    euclidean_norm_path = os.path.join(outdir, f'{collection_name}_counts_normToMax_quantileNorm_euclideanNorm.txt')
    print(f'[{collection_name}] Normalization complete -> {euclidean_norm_path}')

    return euclidean_norm_path


def run_enrichment(collection_name, euclidean_norm_path):

    '''
    Run CHEERS_computeEnrichment.py for every configured trait against a
    collection's specificity scores.

    :param collection_name: The on-disk collection directory name.
    :param euclidean_norm_path: Path to the Euclidean-normalized specificity score file.
    '''

    for trait_name, snp_list_path in config.TRAITS.items():

        outdir = os.path.join(config.ENRICHMENT_RESULTS_DIR, collection_name, trait_name)
        os.makedirs(outdir, exist_ok=True)

        subprocess.run(
            [config.PYTHON_BIN, config.CHEERS_ENRICHMENT, trait_name, outdir, euclidean_norm_path, snp_list_path],
            check=True, env=config.PYTHON_ENV,
        )

        print(f'[{collection_name}] Enrichment ({trait_name}) complete -> {outdir}')


def run_collection(collection_key):

    '''
    Run the full pipeline for one collection.

    :param collection_key: A key of config.COLLECTIONS (e.g. "cad").
    '''

    if collection_key not in config.COLLECTIONS:
        print(f'Unknown collection key: {collection_key}', file=sys.stderr)
        print(f'Valid keys: {", ".join(config.COLLECTIONS.keys())}', file=sys.stderr)
        sys.exit(2)

    collection_name = config.COLLECTIONS[collection_key]

    print(f'=== {collection_key} ({collection_name}) ===')

    counts_dir = quantify_collection(collection_name)
    euclidean_norm_path = run_normalize(collection_name, counts_dir)
    run_enrichment(collection_name, euclidean_norm_path)

    print(f'=== {collection_key} ({collection_name}) done ===')


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(f'Usage: python3 {sys.argv[0]} <collection_key>', file=sys.stderr)
        print(f'Valid keys: {", ".join(config.COLLECTIONS.keys())}', file=sys.stderr)
        sys.exit(2)

    run_collection(sys.argv[1])
