#!/usr/bin/env python3

"""
:File: compare_with_crestgv.py

Compares CHEERS's per-cell-type enrichment p-values against CREST-GV's own
enrichment statistics, for every (collection, trait) pair.

CREST-GV's results live at
    <CRESTGV_RESULTS_BASE>/<trait>/<collection_key>/<collection_key>_statistics_CREST-GV.csv
with columns including CREST-GV (Z statistic), CREST-GV_pval (raw p-value),
CREST-GV_norm (normalized Z statistic), CREST-GV_pval_norm (normalized
p-value) - the *_ps.csv files under COLLECTIONS_BASE are a different thing
(the trait-agnostic background/calibration input to CREST-GV, not its
enrichment result) and are not used here.

For each (collection, trait) pair, joins on cell type and writes a merged
per-cell-type comparison table plus Spearman correlations (against both the
raw and normalized CREST-GV p-value) to results/comparison/. A final summary
table across all pairs is also written.

Usage:
    python3 compare_with_crestgv.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

import pandas as pd
from scipy.stats import spearmanr


def load_cheers_pvalues(collection_dir, trait_name):

    '''
    Load CHEERS's per-cell-type enrichment p-values for one collection/trait.

    :param collection_dir: The on-disk collection directory name (e.g. "CAD").
    :param trait_name: The trait key (e.g. "type_1_diabetes_mellitus").
    :return: A DataFrame with columns ['cellType', 'CHEERS_pval'].
    '''

    path = os.path.join(
        config.ENRICHMENT_RESULTS_DIR, collection_dir, trait_name,
        f'{trait_name}_disease_enrichment_pValues.txt',
    )

    df = pd.read_csv(path, sep='\t', header=None, names=['cellType', 'CHEERS_pval'])
    df['cellType'] = df['cellType'].str.replace('_ReadsInPeaks', '', regex=False)

    return df


def load_crestgv_stats(collection_key, trait_name):

    '''
    Load CREST-GV's per-cell-type enrichment statistics for one collection/trait.

    :param collection_key: The internal collection key (e.g. "cad").
    :param trait_name: The trait key (e.g. "type_1_diabetes_mellitus").
    :return: A DataFrame with columns ['cellType', 'CREST-GV', 'CREST-GV_pval',
        'CREST-GV_norm', 'CREST-GV_pval_norm'].
    '''

    path = os.path.join(
        config.CRESTGV_RESULTS_BASE, trait_name, collection_key,
        f'{collection_key}_statistics_CREST-GV.csv',
    )

    df = pd.read_csv(path, sep='\t')
    cols = ['cellType', 'CREST-GV', 'CREST-GV_pval', 'CREST-GV_norm', 'CREST-GV_pval_norm']

    return df[cols]


def compare_pair(collection_key, trait_name):

    '''
    Compare CHEERS vs CREST-GV for one collection/trait pair.

    :param collection_key: The internal collection key (e.g. "cad").
    :param trait_name: The trait key (e.g. "type_1_diabetes_mellitus").
    :return: A tuple of (merged per-cell-type DataFrame, summary dict) or (None, None) on failure.
    '''

    collection_dir = config.COLLECTIONS[collection_key]

    try:
        cheers_df = load_cheers_pvalues(collection_dir, trait_name)
        crestgv_df = load_crestgv_stats(collection_key, trait_name)
    except FileNotFoundError as e:
        print(f'[{collection_key}/{trait_name}] SKIPPED: {e}', file=sys.stderr)
        return None, None

    merged = pd.merge(cheers_df, crestgv_df, on='cellType', how='inner')

    n_cheers = len(cheers_df)
    n_crestgv = len(crestgv_df)
    n_merged = len(merged)

    if n_merged < n_cheers or n_merged < n_crestgv:
        print(
            f'[{collection_key}/{trait_name}] WARNING: cell type mismatch - '
            f'CHEERS={n_cheers}, CREST-GV={n_crestgv}, matched={n_merged}',
            file=sys.stderr,
        )

    summary = {'collection': collection_key, 'trait': trait_name, 'n_cell_types': n_merged}

    if n_merged >= 3:
        rho_raw, p_raw = spearmanr(merged['CHEERS_pval'], merged['CREST-GV_pval'])
        rho_norm, p_norm = spearmanr(merged['CHEERS_pval'], merged['CREST-GV_pval_norm'])
        summary['spearman_rho_vs_raw'] = rho_raw
        summary['spearman_p_vs_raw'] = p_raw
        summary['spearman_rho_vs_norm'] = rho_norm
        summary['spearman_p_vs_norm'] = p_norm
    else:
        summary['spearman_rho_vs_raw'] = float('nan')
        summary['spearman_p_vs_raw'] = float('nan')
        summary['spearman_rho_vs_norm'] = float('nan')
        summary['spearman_p_vs_norm'] = float('nan')

    for label, col in (('CHEERS', 'CHEERS_pval'), ('CREST-GV_raw', 'CREST-GV_pval'), ('CREST-GV_norm', 'CREST-GV_pval_norm')):
        summary[f'n_significant_{label}_p<0.05'] = int((merged[col] < 0.05).sum())

    summary['n_significant_both_raw_p<0.05'] = int(
        ((merged['CHEERS_pval'] < 0.05) & (merged['CREST-GV_pval'] < 0.05)).sum()
    )
    summary['n_significant_both_norm_p<0.05'] = int(
        ((merged['CHEERS_pval'] < 0.05) & (merged['CREST-GV_pval_norm'] < 0.05)).sum()
    )

    return merged, summary


def main():

    os.makedirs(config.COMPARISON_RESULTS_DIR, exist_ok=True)

    summaries = list()

    for collection_key in config.COLLECTIONS:
        for trait_name in config.TRAITS:

            merged, summary = compare_pair(collection_key, trait_name)

            if merged is None:
                continue

            out_path = os.path.join(
                config.COMPARISON_RESULTS_DIR, f'{collection_key}__{trait_name}.tsv'
            )
            merged.sort_values('CHEERS_pval').to_csv(out_path, sep='\t', index=False)

            summaries.append(summary)
            print(
                f'[{collection_key}/{trait_name}] n={summary["n_cell_types"]} '
                f'spearman_rho(raw)={summary["spearman_rho_vs_raw"]:.3f} '
                f'spearman_rho(norm)={summary["spearman_rho_vs_norm"]:.3f}'
            )

    summary_df = pd.DataFrame(summaries)
    summary_path = os.path.join(config.COMPARISON_RESULTS_DIR, 'summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)

    print(f'\nWrote per-pair comparisons and summary.tsv -> {config.COMPARISON_RESULTS_DIR}')
    print(f'\nOverall median Spearman rho (vs raw CREST-GV p-value): {summary_df["spearman_rho_vs_raw"].median():.3f}')
    print(f'Overall median Spearman rho (vs normalized CREST-GV p-value): {summary_df["spearman_rho_vs_norm"].median():.3f}')


if __name__ == '__main__':

    main()
