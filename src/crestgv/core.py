"""Core CREST-GV implementation.

CREST-GV tests whether a set of trait-associated genetic variants is enriched
in the open chromatin of each cell type of a chromatin accessibility
collection, relative to a trait-agnostic background (1000 Genomes Project
variants). Enrichment is expressed as a Z-score, which is then standardised
within the collection so that cell types from different studies can be
compared on a common scale.

Peak-set semantics
-------------------
Each cell type's variants are intersected against *that cell type's own,
independently called peak set* (``<cell type>_L-tron.bed`` for reference
collections, or ``peaks/<cell type>.bed`` for in-house collections) — not a
consensus/merged peak set shared across the collection. This is a deliberate
design choice: it keeps each cell type's test self-contained and lets new
cell types be added to a collection without recomputing anything for the
existing ones. It does mean CREST-GV and consensus-peak-set methods such as
CHEERS are not scoring literally the same intervals, even within the "same"
collection.
"""

from __future__ import annotations

import os
import random
import sys
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import pyBigWig
import pyranges as pr
from scipy import stats
from tqdm import tqdm

warnings.filterwarnings("ignore")

#: Public collection name -> literal directory name under ``URL`` (the
#: reference data collection root). Kept as a single module-level constant so
#: the CLI and the class can never drift out of sync with each other.
#: ``pancreatic_pbmc`` was removed from the reference atlas (poor data
#: quality) and is intentionally absent here.
COLLECTION_NAME_DICT = {
    "cad": "CAD",
    "calderon": "calderon",
    "catlas_fetal": "catlas_fetal",
    "catlas_adult": "catlas_adult",
    "erythoid_d7_d10_d13_d17": "Days7_10_13_17",
    "greenleaf_brain": "greenleaf_brain",
    "h1_hescs": "H1_hESCs",
    "immune_cell": "immune_cell",
    "ludwig2019": "ludwig2019",
    "mpal": "MPAL",
    "super_pbmc": "super_PBMC",
}

#: Reference data collection root. Kept as a module-level default so it can
#: be overridden (e.g. in tests) without touching the class.
DEFAULT_COLLECTION_URL = "/project/Wellcome_Discovery/datashare/CREST-GV_collection"


class crestgv:
    """CREST-GV enrichment pipeline for a single (genetic file, collection) pair."""

    def __init__(
        self,
        genetic: Optional[str] = None,
        output: str = "output",
        genome: str = "hg38",
        min_number_genetics: int = 500,
        seed: Optional[int] = 42,
        collection_name: str = "",
        in_house_collection_path: str = "",
        collection_url: str = DEFAULT_COLLECTION_URL,
    ):
        """
        Parameters
        ----------
        genetic : str
            Path to the genetic variant file: tab-separated, **no header
            row**, exactly 3 columns in this fixed order: variant ID
            (``rs#####`` or similar), chromosome (``chr#`` or ``#``) and
            1-based position. For example::

                rs60692108          chr10   801748
                chr10-840700-A-C    chr10   840700
        output : str
            Directory where results are written.
        genome : str
            Genome build. Only ``"hg38"`` is currently supported.
        min_number_genetics : int
            Minimum number of variants required after cleaning; the run
            aborts below this threshold.
        seed : int or None
            Seed for reproducibility (passed to ``random`` and ``numpy``).
        collection_name : str
            Reference collection name, one of ``COLLECTION_NAME_DICT``.
            Mutually exclusive with ``in_house_collection_path``.
        in_house_collection_path : str
            Path to an in-house collection
            (``<path>/<collection-name>``), structured as::

                <path>/<collection-name>/
                ├── bigwigs/<cell-type>.bw
                ├── peaks/<cell-type>.bed
                └── <collection-name>_info.csv

            where ``<collection-name>_info.csv`` lists the cell type names
            (one per line, matching the bigwig/peak file stems).
        collection_url : str
            Root directory of the reference data collection. Override only
            if you are not running on the original CCB cluster layout.
        """

        self.seed = seed
        if self.seed is not None and isinstance(self.seed, int):
            random.seed(self.seed)
            np.random.seed(self.seed)

        if genetic is not None and os.path.isfile(genetic):
            self.genetic_file = genetic
        else:
            sys.exit("Provide correct path to genetic file.")

        self.genetic_df = pd.DataFrame()
        self.min_number_genetics = min_number_genetics

        self.output = output
        self.genome = genome

        if self.genome != "hg38":
            if self.genome == "hg19":
                sys.exit("Genome 'hg19' not functional at the moment. Please use genome 'hg38'.")
            sys.exit("Select genome between 'hg19' and 'hg38'.")

        self.URL = collection_url
        self.collection_name_dict = COLLECTION_NAME_DICT

        if not in_house_collection_path and not collection_name:
            sys.exit("Please specify a collection name or a in house collection path.")
        elif not in_house_collection_path:
            if collection_name not in self.collection_name_dict:
                sys.exit(
                    "Wrong 'collection_name' selected, please specify a collection "
                    "name from: [%s]" % (", ".join(self.collection_name_dict.keys()))
                )
            self.collection_name = collection_name
            self.collectionHouseBool = False
        else:
            self.in_house_collection_path = in_house_collection_path
            self.collection_name = in_house_collection_path.rstrip("/").split("/")[-1]
            self.collectionHouseBool = True

        if self.collectionHouseBool:
            # In-house collections have no precomputed background overlap
            # (<collection>_ps.csv), so the 1000 Genomes background variant
            # set is loaded once here and intersected with the collection's
            # own peaks in __prepare_data.
            print("Loading 1000genome snps...")
            chunk_size = 100_000
            chunks = []
            for chunk in pd.read_csv(
                self.URL + os.sep + "1000genomes/ALL_1000_genomes.variants.hg38.bed",
                sep="\t",
                chunksize=chunk_size,
                low_memory=False,
                header=None,
                dtype={0: "category", 1: "int32", 2: "int32", 3: "category"},
            ):
                chunks.append(chunk[[0, 1, 2]])

            snp_data = pd.concat(chunks, ignore_index=True)
            snp_data.columns = ["Chromosome", "Start", "End"]
            self.n_snp_data = snp_data.shape[0]
            self.snps = pr.PyRanges(snp_data)

    def __check_genetic_format(self, df_genetics: pd.DataFrame) -> pd.DataFrame:
        """Validate and clean a raw genetic variant DataFrame.

        Parameters
        ----------
        df_genetics : pandas.DataFrame
            Raw genetic variants, expected to contain at least the columns
            ``CHR_ID``, ``CHR_POS`` and ``SNPS``.

        Returns
        -------
        pandas.DataFrame
            Cleaned variants (invalid rows dropped).
        """

        if not os.path.exists(self.genetic_file):
            sys.exit("Error, genetic file does not exist.")

        print("Checking and cleaning genetic file...")

        required_columns = ["CHR_ID", "CHR_POS", "SNPS"]
        try:
            df_genetics = df_genetics[required_columns]
        except KeyError as e:
            sys.exit(
                f"Error: Missing required column {e}. Genetic file must contain: "
                "'CHR_ID' (int/float), 'CHR_POS' (int/float), 'SNPS' (string)."
            )

        df_genetics["CHR_POS"] = df_genetics["CHR_POS"].astype(str)
        df_genetics = df_genetics[
            ~df_genetics["CHR_POS"].str.contains(";")
            & ~df_genetics["CHR_POS"].str.contains("x", case=False)
        ]
        df_genetics["CHR_POS"] = pd.to_numeric(df_genetics["CHR_POS"], errors="coerce")
        df_genetics = df_genetics.dropna(subset=["CHR_POS", "CHR_ID", "SNPS"])
        df_genetics["CHR_POS"] = df_genetics["CHR_POS"].astype(int)

        if df_genetics.empty:
            sys.exit("Error: No valid data remains after filtering.")

        print(f"Processed genetic file. Retained {len(df_genetics)} rows.")
        return df_genetics

    def __load_genetic(self) -> pd.DataFrame:
        """Load and clean the genetic variant file provided at construction."""

        # No header row in the input file: names= assigns columns by
        # position (SNPS, CHR_ID, CHR_POS, in that order), it does not match
        # against a header. See the `genetic` parameter docstring above.
        df_genetics = pd.read_csv(self.genetic_file, sep="\t", names=["SNPS", "CHR_ID", "CHR_POS"])
        df_genetics = self.__check_genetic_format(df_genetics)
        df_genetics = df_genetics.dropna(subset=["CHR_ID"], axis=0)
        df_genetics["CHR_ID"] = df_genetics["CHR_ID"].astype(str)
        df_genetics["CHR_POS"] = df_genetics["CHR_POS"].astype(int)
        df_genetics["SNPS"] = df_genetics["SNPS"].astype(str)
        df_genetics["CHR_ID"] = "chr" + df_genetics["CHR_ID"].str.replace(".0", "", regex=False)
        df_genetics["CHR_ID"] = df_genetics["CHR_ID"].str.replace("chrchr", "chr", regex=False)
        df_genetics["CHR_POS+1"] = df_genetics["CHR_POS"] + 1
        df_genetics = df_genetics.sort_values(["CHR_ID", "CHR_POS"])
        df_genetics = df_genetics.reset_index(drop=True)
        df_genetics = df_genetics[["CHR_ID", "CHR_POS", "CHR_POS+1", "SNPS"]]
        df_genetics = df_genetics.drop_duplicates()

        print("Total number of used variants in CREST-GV: %s" % df_genetics.shape[0])
        if df_genetics.shape[0] < self.min_number_genetics:
            sys.exit(
                "Genetics provided after quality control contains less than "
                f"{self.min_number_genetics} entry variants."
            )
        self.genetic_df = df_genetics
        return df_genetics

    def __loading_info(self) -> pd.DataFrame:
        """Load the list of cell type names for the selected collection."""

        if self.collectionHouseBool:
            file = "%s/%s_info.csv" % (self.in_house_collection_path, self.collection_name)
        else:
            dir_name = self.collection_name_dict[self.collection_name]
            file = "%s/%s/%s_info.csv" % (self.URL, dir_name, dir_name)
        # Only the cell type name column is required; any further columns
        # (e.g. single-cell counts) in <collection>_info.csv are ignored here.
        info = pd.read_csv(file, sep="\t", names=["cellType", "scCounts"])[["cellType"]]
        return info

    def __loading_bigwigs_and_beds(self, info: pd.DataFrame) -> tuple[list[str], list[str]]:
        """Resolve bigwig/peak file paths for every cell type in ``info``."""

        bigwigs, beds = [], []
        for ct in info["cellType"].tolist():
            if self.collectionHouseBool:
                bigwigs.append("%s/bigwigs/%s.bw" % (self.in_house_collection_path, ct))
                beds.append("%s/peaks/%s.bed" % (self.in_house_collection_path, ct))
            else:
                dir_name = self.collection_name_dict[self.collection_name]
                bigwigs.append("%s/%s/bigwigs/%s.bw" % (self.URL, dir_name, ct))
                beds.append("%s/%s/peaks/%s_L-tron.bed" % (self.URL, dir_name, ct))
        return bigwigs, beds

    def __loading_collection_data(self):
        """Load cell type names and their bigwig/peak file paths."""

        info = self.__loading_info()
        bigwigs, beds = self.__loading_bigwigs_and_beds(info)
        return info, bigwigs, beds

    def __prepare_data(self, beds: list[str], df_genetics: pd.DataFrame) -> pd.DataFrame:
        """Compute per-cell-type background rate and observed variant overlap.

        Parameters
        ----------
        beds : list of str
            Peak file path per cell type (same order as the collection).
        df_genetics : pandas.DataFrame
            Cleaned trait variants, as returned by ``__load_genetic``.

        Returns
        -------
        pandas.DataFrame
            One row per cell type with ``n_snps`` (background variant pool
            size), ``p_succes`` (background overlap proportion) and
            ``GWAS_init`` (observed trait-variant overlap count).
        """

        df_data = pd.DataFrame()

        if self.collectionHouseBool:
            # No precomputed background for in-house collections: intersect
            # the 1000 Genomes background against each cell type's own peaks.
            peak_area, ps, celltypes = [], [], []
            print("Calculating (1) total number of snps from 1000genomes within peaks ...")
            for bed in tqdm(beds):
                region_data = pd.read_csv(bed, sep="\t", names=["Chromosome", "Start", "End"])
                celltype = bed.split("/")[-1].replace("_L-tron.bed", "").replace(".bed", "")
                celltypes.append(celltype)
                regions = pr.PyRanges(region_data)
                intersection = self.snps.intersect(regions)
                p = len(intersection) / self.n_snp_data
                ps.append(p)

            df_data["n_snps"] = int(self.n_snp_data)
            df_data["p_succes"] = ps
            df_data.index = celltypes
        else:
            # Reference collections ship a precomputed background overlap
            # (<collection>_ps.csv) so it does not need recomputing per run.
            dir_name = self.collection_name_dict[self.collection_name]
            info_ps = pd.read_csv(
                "%s/%s/%s_ps.csv" % (self.URL, dir_name, dir_name), index_col=0, sep="\t"
            )
            df_data["n_snps"] = info_ps["n_snps"]
            df_data["p_succes"] = info_ps["p_succes"]

        df_genetics = df_genetics.copy()
        df_genetics.columns = ["Chromosome", "Start", "End", "SNPS"]
        py_df_genetics = pr.PyRanges(df_genetics[["Chromosome", "Start", "End"]])

        xs = []
        print("Intersecting genetic variants with peak regions ...")
        for bed in tqdm(beds):
            region_data = pd.read_csv(bed, sep="\t", names=["Chromosome", "Start", "End"])
            regions = pr.PyRanges(region_data)
            intersection = py_df_genetics.intersect(regions)
            xs.append(len(intersection))
        df_data["GWAS_init"] = xs

        return df_data

    def __add_statistics(self, df_data: pd.DataFrame, number_of_genetic: int) -> pd.DataFrame:
        """Compute the enrichment Z-score, within-collection standardised
        score, and their right-tailed P-values.

        Parameters
        ----------
        df_data : pandas.DataFrame
            Output of ``__prepare_data``.
        number_of_genetic : int
            Number of trait variants used (``len(df_genetics)``).

        Returns
        -------
        pandas.DataFrame
            ``df_data`` with the enrichment statistics appended.
        """

        print("Calculating statistics ...")

        # Exact binomial pmf/sf are kept for diagnostic purposes only: they
        # underflow to exactly 0.0 for strong enrichment in large collections
        # (observed for ludwig2019 on a red-blood-cell trait) and should not
        # be used for ranking or p-value-based comparisons. The supported
        # statistics are the normal-approximation Z/Z_norm based ones below.
        df_data["pmf"] = stats.binom.pmf(df_data["GWAS_init"], number_of_genetic, df_data["p_succes"])
        df_data["sf"] = stats.binom.sf(df_data["GWAS_init"] - 1, number_of_genetic, df_data["p_succes"])

        df_data["check_np"] = number_of_genetic * df_data["p_succes"]
        df_data["check_n(1-p)"] = number_of_genetic * (1 - df_data["p_succes"])

        df_data["p_hat"] = df_data["GWAS_init"] / number_of_genetic

        # Guard against p_succes == 0 or 1 (possible for a degenerate/empty
        # or genome-spanning peak set): SE would be exactly 0 and Z would
        # silently become +-inf/nan. Cell types failing this check are
        # flagged rather than silently propagating a non-finite value.
        valid_p = (df_data["p_succes"] > 0) & (df_data["p_succes"] < 1)
        se = pd.Series(np.nan, index=df_data.index)
        se[valid_p] = np.sqrt(
            (df_data.loc[valid_p, "p_succes"] * (1 - df_data.loc[valid_p, "p_succes"])) / number_of_genetic
        )
        df_data["SE"] = se

        z = pd.Series(np.nan, index=df_data.index)
        z[valid_p] = (df_data.loc[valid_p, "p_hat"] - df_data.loc[valid_p, "p_succes"]) / se[valid_p]
        df_data["Z"] = z
        df_data["pdf"] = stats.norm.sf(df_data["Z"])

        if (~valid_p).any():
            warnings.warn(
                "%d cell type(s) had a degenerate background rate (p_succes "
                "== 0 or 1) and were excluded from the enrichment statistics: %s"
                % ((~valid_p).sum(), ", ".join(df_data.index[~valid_p].astype(str)))
            )

        # Within-collection standardisation. Guard against a zero (or
        # near-zero) standard deviation, which otherwise blows Z_norm up to
        # +-inf for very small or near-flat collections (e.g. h1_hescs, 4
        # cell types). Cell types in a collection with no usable spread are
        # reported as Z_norm = 0 (no information to rank them against each
        # other) rather than a meaningless infinite value.
        z_mean = df_data["Z"].mean()
        z_std = df_data["Z"].std()
        if not np.isfinite(z_std) or z_std == 0:
            warnings.warn(
                "Standard deviation of Z across the collection is zero or "
                "undefined (collection may be too small or all cell types "
                "have identical enrichment); Z_norm set to 0 for all cell types."
            )
            df_data["Z_norm"] = 0.0
        else:
            df_data["Z_norm"] = (df_data["Z"] - z_mean) / z_std
        df_data["pdf_cgv_norm"] = stats.norm.sf(df_data["Z_norm"])

        df_data["CREST-GV"] = df_data["Z"]
        df_data["CREST-GV_pval"] = df_data["pdf"]

        df_data["CREST-GV_norm"] = df_data["Z_norm"]
        df_data["CREST-GV_pval_norm"] = df_data["pdf_cgv_norm"]

        return df_data

    def calculate_enrichment_score(self) -> pd.DataFrame:
        """Run the full enrichment pipeline and write ``<output>/<collection>_statistics_CREST-GV.csv``.

        Returns
        -------
        pandas.DataFrame
            Per-cell-type enrichment statistics (also written to disk).
        """

        df_genetics = self.__load_genetic()
        info, _, beds = self.__loading_collection_data()
        df_collection = self.__prepare_data(beds, df_genetics)
        df_collection = self.__add_statistics(df_collection, df_genetics.shape[0])

        if not os.path.exists(self.output):
            os.makedirs(self.output)
        df_collection.to_csv(
            self.output + os.sep + "%s_statistics_CREST-GV.csv" % self.collection_name, sep="\t"
        )

        print("Processing genetics finished.")
        return df_collection

    def get_coverage(self) -> pd.DataFrame:
        """Compute per-variant bigwig coverage for every cell type.

        This is a separate, much heavier diagnostic output (one column per
        cell type, one row per input variant) and is not required for the
        enrichment statistics above; call it only if you need raw per-base
        signal at each variant.

        Returns
        -------
        pandas.DataFrame
            Input variants with one appended column of coverage per cell type
            (also written to ``<output>/coverage.csv``).
        """

        self.genetic_df = self.__load_genetic()
        self.genetic_df["CHR_POS+1"] = self.genetic_df["CHR_POS"] + 1

        info, bigwigs, _ = self.__loading_collection_data()

        print("Getting coverage ...")
        for ct in tqdm(info["cellType"].tolist()):
            if self.collectionHouseBool:
                bigwig = "%s/bigwigs/%s.bw" % (self.in_house_collection_path, ct)
            else:
                dir_name = self.collection_name_dict[self.collection_name]
                bigwig = "%s/%s/bigwigs/%s.bw" % (self.URL, dir_name, ct)

            bw = pyBigWig.open(bigwig)
            values_list = [
                bw.values(c, s, e, numpy=True).item()
                for c, s, e in zip(
                    self.genetic_df["CHR_ID"].tolist(),
                    self.genetic_df["CHR_POS"].tolist(),
                    self.genetic_df["CHR_POS+1"].tolist(),
                )
            ]
            bw.close()
            self.genetic_df[ct] = values_list

        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.genetic_df.to_csv(self.output + os.sep + "coverage.csv", sep="\t")
        return self.genetic_df
