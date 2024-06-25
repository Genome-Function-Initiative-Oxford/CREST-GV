# ESVAR

Enrichment Score for VARiants


### edit readme

For the time being, it is working only in CentOS (CEPH)
Module to load:
	- module load python-cbrg
	- module load bedtools



--- From terminal ---

```
import sys
import pandas as pd

sys.path.append('ensgen/')
from ensgen import ensgen

genetic_file = "<path>/<file>.<ext>"


# catlas adult

es = ensgen(genetic=genetic_file, seed=42)
df = es.process_genetics_for_catlas(origin='adult')
es.plot_genetics_for_catlas(show=False)

es.multicoverage_genetics_for_catlas(variants=genetic_file, origin='adult')

groupA = ["<cell type of interest>"]
groupB = ["other"]
es.plot_genetics_enrichments(groupA=groupA, groupB=groupB, show=False)

##########


# catlas fetal

es = ensgen(genetic=genetic_file, seed=42)
df = es.process_genetics_for_catlas(origin='fetal')
es.plot_genetics_for_catlas(show=False)

es.multicoverage_genetics_for_catlas(variants=genetic_file, origin='fetal')

groupA = ["<cell type of interest>"]
groupB = ["other"]
es.plot_genetics_enrichments(groupA=groupA, groupB=groupB, show=False)

##########


# calderon

es = ensgen(genetic=genetic_file, seed=42)
df = es.process_genetics_for_calderon()
es.plot_genetics_for_calderon(show=False)

es.multicoverage_genetics_for_calderon(variants=genetic_file)

groupA = ["<cell type of interest>"]
groupB = ["other"]
es.plot_genetics_enrichments(groupA=groupA, groupB=groupB, show=False)

##########


# ludwig2019

es = ensgen(genetic=genetic_file, seed=42)
df = es.process_genetics_for_ludwig2019()
es.plot_genetics_for_ludwig2019(show=False)

es.multicoverage_genetics_for_ludwig2019(variants=genetic_file)

groupA = ["<cell type of interest>"]
groupB = ["other"]
es.plot_genetics_enrichments(groupA=groupA, groupB=groupB, show=False)

##########


# MPAL_lowGr

es = ensgen(genetic=genetic_file, seed=42)
df = es.process_genetics_for_MPAL_lowGr()
es.plot_genetics_for_MPAL_lowGr(show=False)

es.multicoverage_genetics_for_MPAL_lowGr(variants=genetic_file)

groupA = ["<cell type of interest>"]
groupB = ["other"]
es.plot_genetics_enrichments(groupA=groupA, groupB=groupB, show=False)

##########
```


--- From API ---

```
# catlas adult
python ensgen.py --catlas-adult True -g <path>/<file>.<ext> -o api

# catlas fetal
python ensgen.py --catlas-fetal True -g <path>/<file>.<ext> -o api

# calderon
python ensgen.py --calderon True -g <path>/<file>.<ext> -o api

# ludwig2019
python ensgen.py --ludwig2019 True -g <path>/<file>.<ext> -o api

# MPAL_lowGr
python ensgen.py --MPAL_lowGr True -g <path>/<file>.<ext> -o api


python esvar.py -g /project/Wellcome_Discovery/sriva/Git/ESVAR/genetics_test/Bcell_GO_0035456.tsv -cn super_pbmc -o output_api


```
