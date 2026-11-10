"""End-to-end smoke test against the bundled in-house example data.

Does not require network access or the reference data collection: it only
exercises the in-house collection code path, using
``example_files/in-house-data`` and ``example_files/genetic.vcf``.
"""

import os
import subprocess
import sys

import pandas as pd
import pytest

from crestgv import crestgv

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENETIC_FILE = os.path.join(REPO_ROOT, "example_files", "genetic.tsv")
IN_HOUSE_PATH = os.path.join(REPO_ROOT, "example_files", "in-house-data")


def test_calculate_enrichment_score(tmp_path):
    cgv = crestgv(
        genetic=GENETIC_FILE,
        output=str(tmp_path),
        min_number_genetics=1,  # the bundled example file is small
        in_house_collection_path=IN_HOUSE_PATH,
    )
    df = cgv.calculate_enrichment_score()

    expected_cell_types = {"cell_type_1", "cell_type_2", "cell_type_3"}
    assert set(df.index) == expected_cell_types
    for col in ["GWAS_init", "p_succes", "Z", "CREST-GV", "CREST-GV_pval", "CREST-GV_norm", "CREST-GV_pval_norm"]:
        assert col in df.columns
    assert df["CREST-GV_norm"].notna().all()

    out_file = tmp_path / "in-house-data_statistics_CREST-GV.csv"
    assert out_file.exists()
    reloaded = pd.read_csv(out_file, sep="\t", index_col=0)
    assert set(reloaded.index) == expected_cell_types


def test_get_coverage(tmp_path):
    cgv = crestgv(
        genetic=GENETIC_FILE,
        output=str(tmp_path),
        min_number_genetics=1,
        in_house_collection_path=IN_HOUSE_PATH,
    )
    df = cgv.get_coverage()
    for ct in ["cell_type_1", "cell_type_2", "cell_type_3"]:
        assert ct in df.columns
    assert (tmp_path / "coverage.csv").exists()


def test_missing_collection_name_exits():
    with pytest.raises(SystemExit):
        crestgv(genetic=GENETIC_FILE)


def test_unknown_collection_name_exits():
    with pytest.raises(SystemExit):
        crestgv(genetic=GENETIC_FILE, collection_name="not_a_real_collection")


def test_cli_end_to_end(tmp_path):
    result = subprocess.run(
        [
            sys.executable, "-m", "crestgv.cli",
            "-g", GENETIC_FILE,
            "-o", str(tmp_path),
            "-ng", "1",
            "-incp", IN_HOUSE_PATH,
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "in-house-data_statistics_CREST-GV.csv").exists()


def test_cli_list_collections():
    result = subprocess.run(
        [sys.executable, "-m", "crestgv.cli", "--list_collections"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "greenleaf_brain" in result.stdout
    assert "pancreatic_pbmc" not in result.stdout
