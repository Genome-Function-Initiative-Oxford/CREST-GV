#!/usr/bin/env python3

"""
:File: run_all.py

Runs the full CREST-GV_vs_CHEERS pipeline (build consensus peaks -> quantify
bigWigs -> CHEERS_normalize.py -> CHEERS_computeEnrichment.py per trait) for
every configured collection.

Usage:
    python3 run_all.py [collection_key ...]

With no arguments, processes all collections in config.COLLECTIONS. Each
collection is processed independently; a failure in one is logged and does
not stop the others.
"""

import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

from run_collection import run_collection


def main():

    collection_keys = sys.argv[1:] or list(config.COLLECTIONS.keys())

    results = dict()

    for collection_key in collection_keys:
        try:
            run_collection(collection_key)
            results[collection_key] = 'OK'
        except Exception:
            print(f'[{collection_key}] FAILED', file=sys.stderr)
            traceback.print_exc()
            results[collection_key] = 'FAILED'

    print('\n=== Summary ===')
    for collection_key, status in results.items():
        print(f'{collection_key}: {status}')

    if any(status == 'FAILED' for status in results.values()):
        sys.exit(1)


if __name__ == '__main__':

    main()
