"""Builds quickstart.ipynb for the CREST-GV repo. Run from this directory
(examples/), with the package installed (pip install -e ..), then execute it
with nbconvert so the published notebook shows real output:

    cd examples
    python3 build_notebook.py
    jupyter nbconvert --to notebook --execute --inplace quickstart.ipynb
"""

import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# CREST-GV quickstart

CREST-GV tests whether a set of trait-associated genetic variants is enriched in the
open chromatin of each cell type in a chromatin accessibility collection, relative to
a trait-agnostic background (1000 Genomes Project variants).

This notebook walks through the library API end to end against the small example
dataset bundled with the package (`../example_files/`) — three synthetic cell types, so
it runs in a few seconds and needs no external data. The same calls work against the
full reference atlas; see the final section for that.
"""))

cells.append(nbf.v4.new_code_cell(
"""from crestgv import crestgv, COLLECTION_NAME_DICT

print("Reference collections available in the atlas:")
for key, on_disk_name in COLLECTION_NAME_DICT.items():
    print(f"  {key:25s} -> {on_disk_name}")
"""))

cells.append(nbf.v4.new_markdown_cell(
"""## 1. The genetic variant file

The only mandatory input is a tab-separated file with **no header row** and exactly
three columns: variant ID, chromosome (`chr#`), 1-based position.
"""))

cells.append(nbf.v4.new_code_cell(
"""with open("../example_files/genetic.tsv") as f:
    for line in f:
        print(line.rstrip())
"""))

cells.append(nbf.v4.new_markdown_cell(
"""## 2. Running against an in-house collection

An in-house collection is a folder with `bigwigs/<cell_type>.bw`,
`peaks/<cell_type>_L-tron.bed` (or any 3-column BED of peak calls), and a
`<collection>_info.csv`. See [`../example_files/in-house-data`](../../example_files/in-house-data)
for the exact layout.
"""))

cells.append(nbf.v4.new_code_cell(
"""cgv = crestgv(
    genetic="../example_files/genetic.tsv",
    output="../example_output",
    min_number_genetics=1,  # the bundled example file is tiny; real runs use the default (500)
    in_house_collection_path="../example_files/in-house-data",
)
results = cgv.calculate_enrichment_score()
results
"""))

cells.append(nbf.v4.new_markdown_cell(
"""Key columns:

- **`CREST-GV`** / **`CREST-GV_pval`** — the raw enrichment Z-score and its P-value,
  calibrated against the genome-wide background. Comparable across cell types *within*
  this collection, but not across collections of different size/composition.
- **`CREST-GV_norm`** / **`CREST-GV_pval_norm`** — the same Z-score standardised
  within this collection (subtract the mean, divide by the standard deviation across
  all of this collection's cell types). This is the statistic used throughout the
  CREST-GV paper, since it's what makes cell types from different collections
  comparable on a common scale.

A results file is also written to disk:
"""))

cells.append(nbf.v4.new_code_cell(
"""import pandas as pd

reloaded = pd.read_csv("../example_output/in-house-data_statistics_CREST-GV.csv", sep="\\t", index_col=0)
reloaded
"""))

cells.append(nbf.v4.new_markdown_cell(
"""## 3. Ranking cell types

The standardised score (`CREST-GV_norm`) is what you rank cell types by — higher means
more enriched for this trait's variants relative to the other cell types in the
collection.
"""))

cells.append(nbf.v4.new_code_cell(
"""import matplotlib.pyplot as plt

ranked = results["CREST-GV_norm"].sort_values(ascending=False)
colors = ["#2a78d6" if v >= 0 else "#e34948" for v in ranked.values]

fig, ax = plt.subplots(figsize=(5, 2.5))
ax.barh(ranked.index[::-1], ranked.values[::-1], color=colors[::-1])
ax.axvline(0, color="#c3c2b7", linewidth=1)
ax.set_xlabel(r"Standardised enrichment score ($\\tilde{Z}$)")
fig.tight_layout()
plt.show()
"""))

cells.append(nbf.v4.new_markdown_cell(
"""## 4. Coverage: per-variant signal at each cell type's bigWig

`get_coverage()` reports the raw bigWig signal under each input variant, per cell
type — useful for inspecting individual loci rather than the aggregate enrichment
statistic.
"""))

cells.append(nbf.v4.new_code_cell(
"""coverage = cgv.get_coverage()
coverage
"""))

cells.append(nbf.v4.new_markdown_cell(
"""## 5. Running against the reference atlas instead

Pass `collection_name` (one of the keys printed in Section 0) instead of
`in_house_collection_path` to run against one of the 11 precomputed reference
collections. This requires access to the reference data collection itself (not
bundled with the package — see the main README's Installation section), so it isn't
run here, but the call is otherwise identical:

```python
cgv = crestgv(
    genetic="path/to/your_trait_variants.tsv",
    output="your_output_directory",
    collection_name="cad",
)
results = cgv.calculate_enrichment_score()
```

Or equivalently from the command line:

```bash
crestgv -g path/to/your_trait_variants.tsv -o your_output_directory -cn cad
```

Pass `-ra`/`--run_all` to run against every reference collection in one call, or
`-wc`/`--with_coverage` to also compute per-variant coverage.
"""))

nb["cells"] = cells
with open("quickstart.ipynb", "w") as f:
    nbf.write(nb, f)
print("wrote quickstart.ipynb")
