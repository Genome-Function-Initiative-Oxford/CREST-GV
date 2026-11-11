"""Command-line entry point for CREST-GV.

Fixes applied relative to the pre-0.3.0 CLI:
  * removed the ``number_of_folds``/``lessNG``/``greater25k`` flags, which no
    longer correspond to anything in ``crestgv.core.crestgv`` and made every
    invocation crash with ``TypeError``.
  * the collection name list is now imported from ``crestgv.core`` instead of
    being duplicated here (the duplicate was missing ``greenleaf_brain``, so
    ``--run-all`` silently skipped it).
  * boolean flags (``--with-coverage``, ``--run-all``) now use
    ``action="store_true"`` instead of ``type=bool``, which is a well-known
    argparse footgun: ``type=bool`` calls ``bool("False")``, which is
    ``True`` — passing ``--flag False`` could never actually disable a flag.
  * ``--min-number-genetics`` default changed from 100 to 500, matching the
    class default and the published methodology.
  * coverage calculation is no longer run by default (it is a large,
    separate diagnostic output); pass ``--with-coverage`` to also compute it.
"""

from __future__ import annotations

import argparse
import os

from .core import COLLECTION_NAME_DICT, crestgv


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="crestgv",
        description="CREST-GV: Enrichment Score for selecTive Genetic Variants",
    )

    parser.add_argument(
        "-g", "--genetic", type=str, default=None,
        help="Path to the genetic variant file: tab-separated, no header, columns "
        "in this order: SNPS, CHR_ID, CHR_POS.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output",
        help="Directory where results are written.",
    )
    parser.add_argument(
        "-gb", "--genome", type=str, default="hg38",
        help="Genome build. Only 'hg38' is currently supported.",
    )
    parser.add_argument(
        "-ng", "--min_number_genetics", type=int, default=500,
        help="Minimum number of variants required after cleaning.",
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42,
        help="Seed for reproducibility.",
    )
    parser.add_argument(
        "-cn", "--collection_name", type=str, default="",
        help="Reference collection name. Available: %s. Mutually exclusive "
        "with --in_house_collection_path." % ", ".join(COLLECTION_NAME_DICT.keys()),
    )
    parser.add_argument(
        "-incp", "--in_house_collection_path", type=str, default="",
        help="Path to an in-house collection directory (<path>/<collection-name>).",
    )
    parser.add_argument(
        "-wc", "--with_coverage", action="store_true",
        help="Also compute per-variant bigwig coverage (large additional output).",
    )
    parser.add_argument(
        "-ra", "--run_all", action="store_true",
        help="Run CREST-GV for every reference collection. Ignores --collection_name "
        "and --in_house_collection_path.",
    )
    parser.add_argument(
        "--list_collections", action="store_true",
        help="Print the available reference collection names and exit.",
    )

    return parser


def _run_one(
    args: argparse.Namespace, output: str, collection_name: str = "", in_house_collection_path: str = ""
) -> None:
    cgv = crestgv(
        genetic=args.genetic,
        output=output,
        genome=args.genome,
        min_number_genetics=args.min_number_genetics,
        seed=args.seed,
        collection_name=collection_name,
        in_house_collection_path=in_house_collection_path,
    )
    cgv.calculate_enrichment_score()
    if args.with_coverage:
        cgv.get_coverage()


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.list_collections:
        for name in COLLECTION_NAME_DICT:
            print(name)
        return

    if args.genetic is None:
        parser.error("the following argument is required: -g/--genetic")

    if args.run_all:
        for collection_name in COLLECTION_NAME_DICT:
            _run_one(args, output=args.output + os.sep + collection_name, collection_name=collection_name)
    elif args.in_house_collection_path:
        _run_one(args, output=args.output, in_house_collection_path=args.in_house_collection_path)
    else:
        _run_one(args, output=args.output, collection_name=args.collection_name)


if __name__ == "__main__":
    main()
