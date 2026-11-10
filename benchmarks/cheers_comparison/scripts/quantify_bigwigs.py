#!/usr/bin/env python3

"""
:File: quantify_bigwigs.py

Quantifies per-peak signal for every cell type's bigWig in a CREST-GV
collection, against that collection's consensus peak set, and writes one
CHEERS-format read-count file per cell type (chr, start, end, value;
tab-separated, no header - the format CHEERS_normalize.py expects).

Signal is taken directly from the RPKM-normalized bigWigs via
bigWigAverageOverBed's "sum" column (total signal over all bases in the
peak). No RPKM reversal is applied - see config.py docstring for why this
is safe for CHEERS's downstream normalization.

Usage:
    python3 quantify_bigwigs.py <collection_name>
"""

import glob
import math
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

from build_consensus_peaks import build_consensus_peaks


def load_peak_coords(consensus_named_path):

    '''
    Load the peak ID -> (chrom, start, end) mapping from the named consensus BED.

    :param consensus_named_path: Path to the 4-column named consensus BED.
    :return: A dict mapping peak ID to (chrom, start, end).
    '''

    peak_coords = dict()

    with open(consensus_named_path) as f_in:
        for line in f_in:
            chrom, start, end, peak_id = line.rstrip('\n').split('\t')
            peak_coords[peak_id] = (chrom, start, end)

    return peak_coords


def quantify_bigwig(bigwig_path, consensus_named_path, peak_coords, out_path):

    '''
    Run bigWigAverageOverBed for one bigWig against the consensus peak set,
    and write a CHEERS-format read-count file.

    :param bigwig_path: Path to the input bigWig file.
    :param consensus_named_path: Path to the 4-column named consensus BED.
    :param peak_coords: Dict mapping peak ID to (chrom, start, end).
    :param out_path: Output path for the CHEERS-format count file.
    '''

    with tempfile.NamedTemporaryFile(mode='r', suffix='.tab', delete=False) as tmp_f:
        tmp_path = tmp_f.name

    try:
        subprocess.run(
            [config.BIGWIG_AVERAGE_OVER_BED_BIN, bigwig_path, consensus_named_path, tmp_path],
            check=True,
        )

        rows = list()
        with open(tmp_path) as f_in:
            for line in f_in:
                fields = line.rstrip('\n').split('\t')
                peak_id, signal_sum = fields[0], float(fields[3])
                chrom, start, end = peak_coords[peak_id]

                if not math.isfinite(signal_sum) or signal_sum < 0 or signal_sum > config.SANE_SUM_CEILING:
                    print(
                        f'WARNING: corrupted bigWigAverageOverBed output for {bigwig_path} at '
                        f'{peak_id} (sum={signal_sum!r}, likely a tool bug) - treating as 0.0',
                        file=sys.stderr,
                    )
                    signal_sum = 0.0

                rows.append((chrom, start, int(end), int(start), signal_sum))

        # Sort by chrom, start, end to match CHEERS's expected ordering
        rows.sort(key=lambda r: (r[0], r[3], r[2]))

        with open(out_path, 'w') as f_out:
            for chrom, start, end, start_int, signal_sum in rows:
                f_out.write(f'{chrom}\t{start}\t{end}\t{signal_sum}\n')

    finally:
        os.remove(tmp_path)


def quantify_collection(collection_name):

    '''
    Build the consensus peak set (if needed) and quantify every cell type's
    bigWig in a collection against it.

    :param collection_name: The on-disk collection directory name (e.g. "CAD").
    :return: Directory containing the per-cell-type CHEERS-format count files.
    '''

    _, consensus_named_path = build_consensus_peaks(collection_name)
    peak_coords = load_peak_coords(consensus_named_path)

    bigwigs_dir = os.path.join(config.COLLECTIONS_BASE, collection_name, 'bigwigs')
    bigwig_files = sorted(glob.glob(os.path.join(bigwigs_dir, '*.bw')))

    if not bigwig_files:
        print(f'No bigWig files found in {bigwigs_dir}', file=sys.stderr)
        sys.exit(1)

    out_dir = os.path.join(config.COUNTS_DIR, collection_name)
    os.makedirs(out_dir, exist_ok=True)

    for bigwig_path in bigwig_files:
        cell_type = os.path.basename(bigwig_path)[:-len('.bw')]
        out_path = os.path.join(out_dir, f'{cell_type}_ReadsInPeaks.txt')
        quantify_bigwig(bigwig_path, consensus_named_path, peak_coords, out_path)
        print(f'[{collection_name}] {cell_type} -> {out_path}')

    return out_dir


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(f'Usage: python3 {sys.argv[0]} <collection_name>', file=sys.stderr)
        sys.exit(2)

    quantify_collection(sys.argv[1])
