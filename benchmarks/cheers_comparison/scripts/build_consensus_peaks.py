#!/usr/bin/env python3

"""
:File: build_consensus_peaks.py

Builds a consensus peak set for one CREST-GV collection by merging the
per-cell-type peak BED files (collection/peaks/*_L-tron.bed) into a single
sorted, non-overlapping set of intervals. CHEERS requires every sample's
read-count file to share identical peak coordinates, but CREST-GV calls
peaks independently per cell type, so a union/merge step is required first.

Usage:
    python3 build_consensus_peaks.py <collection_name>
"""

import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


def build_consensus_peaks(collection_name):

    '''
    Merge all per-cell-type peak BED files in a collection into one sorted,
    non-overlapping consensus peak set, and write it both as a plain 3-column
    BED and as a 4-column BED with a unique peak ID (required as input to
    bigWigAverageOverBed, which needs a unique "name" field).

    :param collection_name: The on-disk collection directory name (e.g. "CAD").
    :return: Tuple of (plain consensus BED path, named consensus BED path).
    '''

    peaks_dir = os.path.join(config.COLLECTIONS_BASE, collection_name, 'peaks')
    peak_files = sorted(glob.glob(os.path.join(peaks_dir, '*.bed')))

    if not peak_files:
        print(f'No peak files found in {peaks_dir}', file=sys.stderr)
        sys.exit(1)

    os.makedirs(config.CONSENSUS_PEAKS_DIR, exist_ok=True)

    combined_path = os.path.join(config.CONSENSUS_PEAKS_DIR, f'{collection_name}_combined.tmp.bed')
    sorted_path = os.path.join(config.CONSENSUS_PEAKS_DIR, f'{collection_name}_sorted.tmp.bed')
    consensus_path = os.path.join(config.CONSENSUS_PEAKS_DIR, f'{collection_name}.bed')
    consensus_named_path = os.path.join(config.CONSENSUS_PEAKS_DIR, f'{collection_name}_named.bed')

    # Concatenate every cell type's peak calls, keeping only the first 3 columns
    with open(combined_path, 'w') as f_out:
        for peak_file in peak_files:
            with open(peak_file) as f_in:
                for line in f_in:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 3:
                        continue
                    f_out.write('\t'.join(fields[:3]) + '\n')

    # Sort and merge overlapping intervals into one consensus set
    with open(sorted_path, 'w') as f_out:
        subprocess.run([config.BEDTOOLS_BIN, 'sort', '-i', combined_path], stdout=f_out, check=True)

    with open(consensus_path, 'w') as f_out:
        subprocess.run([config.BEDTOOLS_BIN, 'merge', '-i', sorted_path], stdout=f_out, check=True)

    # Add a unique peak ID as a 4th column for bigWigAverageOverBed
    with open(consensus_path) as f_in, open(consensus_named_path, 'w') as f_out:
        for line in f_in:
            chrom, start, end = line.rstrip('\n').split('\t')
            peak_id = f'{chrom}_{start}_{end}'
            f_out.write(f'{chrom}\t{start}\t{end}\t{peak_id}\n')

    os.remove(combined_path)
    os.remove(sorted_path)

    n_peaks = sum(1 for _ in open(consensus_path))
    print(f'[{collection_name}] Consensus peak set: {n_peaks} peaks from {len(peak_files)} cell types -> {consensus_named_path}')

    return consensus_path, consensus_named_path


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(f'Usage: python3 {sys.argv[0]} <collection_name>', file=sys.stderr)
        sys.exit(2)

    build_consensus_peaks(sys.argv[1])
