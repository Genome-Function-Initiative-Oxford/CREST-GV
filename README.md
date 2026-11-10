# CREST-GV - Cell types Ranking and Enrichment Score for selecTive Genetic Variants

CREST-GV tests whether a set of trait-associated genetic variants is enriched in the open
chromatin of each cell type of a chromatin accessibility collection, relative to a
trait-agnostic background (1000 Genomes Project variants). It can query the
[reference data collection](our_collection/data_collection.txt) of published bulk and
single-cell ATAC-seq datasets, leveraging [LanceOtron](https://github.com/LHentges/LanceOtron)
peak calls, or your own in-house data (see [In-house data format](#in-house-data-format)).

The only mandatory input is a genetic variant file in the format described in
[Genetic variant file format](#genetic-variant-file-format).

***

## Installation

### From PyPI (recommended)

```bash
pip install crestgv
```

This installs the `crestgv` Python package and the `crestgv` command-line tool.
`pyBigWig` and `pyranges` wheels are platform-dependent; if `pip install` fails to build
them, install via conda/mamba instead (see below) or install `pyBigWig`/`pyranges` from
conda-forge/bioconda first, then `pip install crestgv` in the same environment.

### From source

```bash
git clone https://github.com/Genome-Function-Initiative-Oxford/CREST-GV.git
cd CREST-GV
pip install -e .
```

### Conda/mamba environment (alternative)

```bash
git clone https://github.com/Genome-Function-Initiative-Oxford/CREST-GV.git
cd CREST-GV
mamba env create --file=envs/environment.yml   # or: conda env create --file=envs/environment.yml
conda activate crestgv
pip install -e .
```

### Note on peak calling

CREST-GV's statistics are computed over existing peak calls; it does not call peaks
itself. Peak files for the reference collections must be generated with
[LanceOtron](https://github.com/LHentges/LanceOtron) (the full LanceOtron output columns
are not required by CREST-GV itself, but the `<cell-type>_L-tron.bed`/`.bed` naming and
3-column `chrom, start, end` content are). For in-house collections you may use any peak
caller of your choice, as long as each peak file is a plain tab-separated
`chrom, start, end` BED file.

***

## How to use CREST-GV

### As a library

```python
from crestgv import crestgv

cgv = crestgv(
    genetic="path/to/genetic.tsv",
    output="path/to/output_directory",
    collection_name="cad",   # or in_house_collection_path="path/to/my-collection"
)
df = cgv.calculate_enrichment_score()
```

For the full list of accepted `collection_name` values:

```python
from crestgv import COLLECTION_NAME_DICT
print(list(COLLECTION_NAME_DICT))
```

Per-parameter help:

```python
help(crestgv)
```

### From the command line

```bash
crestgv -g path/to/genetic.tsv -o path/to/output_directory -cn cad
```

```
usage: crestgv [-h] [-g GENETIC] [-o OUTPUT] [-gb GENOME]
                [-ng MIN_NUMBER_GENETICS] [-s SEED] [-cn COLLECTION_NAME]
                [-incp IN_HOUSE_COLLECTION_PATH] [-wc] [-ra]
                [--list_collections]

  -g, --genetic                  Path to the genetic variant file (required).
  -o, --output                   Output directory (default: "output").
  -gb, --genome                  Genome build; only "hg38" is supported.
  -ng, --min_number_genetics     Minimum variants required after cleaning (default: 500).
  -s, --seed                     Random seed (default: 42).
  -cn, --collection_name         Reference collection name.
  -incp, --in_house_collection_path   Path to an in-house collection directory.
  -wc, --with_coverage           Also compute per-variant bigwig coverage (large output).
  -ra, --run_all                 Run for every reference collection.
  --list_collections              Print available collection names and exit.
```

***

## Genetic variant file format

Tab-separated, **no header row**, exactly 3 columns in this order:

1. Variant ID (`rs#####`, or `chr#-pos-ref-alt`)
2. Chromosome, as `chr#`
3. 1-based position

```
rs60692108          chr10   801748
rs74876360          chr10   823912
chr10-840700-A-C    chr10   840700
```

See [`example_files/genetic.tsv`](example_files/genetic.tsv) for a minimal working
example.

***

## In-house data format

All the data for an in-house collection must be stored in a folder named after the
collection, structured as (see
[`example_files/in-house-data`](example_files/in-house-data)):

```
in-house-data/
├── bigwigs
│   ├── cell_type_1.bw
│   ├── cell_type_2.bw
│   └── cell_type_3.bw
├── in-house-data_info.csv
└── peaks
    ├── cell_type_1.bed
    ├── cell_type_2.bed
    └── cell_type_3.bed
```

`in-house-data_info.csv` lists the cell type names, one per line, matching the bigwig and
peak file stems. Pass the path to this folder as `in_house_collection_path`.

***

## Peak-set semantics

Each cell type's variants are intersected against that cell type's own, independently
called peak set, not a consensus/merged peak set shared across the collection. This keeps
each cell type's enrichment test self-contained and lets new cell types be added to a
collection without recomputing anything for the existing ones. If you are comparing
CREST-GV's output against a method that requires a consensus peak set (e.g. CHEERS), be
aware the two methods are not scoring literally the same intervals even within the "same"
collection.

***

## Output statistics

`calculate_enrichment_score()` writes `<output>/<collection_name>_statistics_CREST-GV.csv`
with one row per cell type. The supported statistics (used for ranking and significance)
are the normal-approximation ones:

| Column | Meaning |
|---|---|
| `CREST-GV` | Raw enrichment Z-score for this cell type, against the trait-agnostic background. |
| `CREST-GV_pval` | Right-tailed P-value for `CREST-GV`. |
| `CREST-GV_norm` | `CREST-GV` standardised within the collection (subtract the collection mean, divide by the collection standard deviation). Comparable across collections of different size/scale. |
| `CREST-GV_pval_norm` | Right-tailed P-value for `CREST-GV_norm`. |

`pmf`/`sf` (exact binomial probability/survival function) are also written but are
**diagnostic only**: they underflow to exactly `0.0` for strong enrichment in large
collections and should not be used for ranking or downstream comparison.

***

### Contact us

If you have any suggestions, spot any errors, or have any questions, please do not
hesitate to contact us, or open an
[issue](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/issues).

:email: [simone.riva@imm.ox.ac.uk](mailto:simone.riva@imm.ox.ac.uk)
