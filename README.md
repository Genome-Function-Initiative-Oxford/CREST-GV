# ESVAR - Enrichment Score for genetic VARiants

ESVAR is a method which allows querying a data collection of ~490 cell types (we keep piling more data to add to the data collection) to determine the enrichment score of a set of genetic variants. ESVAR rely on peak properties for each cell type in the data collection leveraging [LanceOtron](https://github.com/LHentges/LanceOtron) peak caller.
ESVAR can also query your personal (in-house) data where formatted correctly (see [In-house data format](#In-house-data-format) section for a properly formatted data structure), in this case, the user can rely on the peak caller of their own choice.
The only mandatory input for the tool is the path of a genetic file, stored following the format shown in [Genetic variant file format](#genetic-variant-file-format) section.

***

## Getting started
ESVAR can be uesed and run using the __esvar__ conda envrironment. Please follow the installation instruction detailed below.

### Installation instructions for conda environment

This section rely on the assumption that any distribution of ```Conda``` (e.g., [Anaconda](https://www.anaconda.com), [Miniconda](https://docs.conda.io/en/latest/miniconda.html), ...) or [```Mamba```](https://mamba.readthedocs.io/en/latest/index.html) is already installed.

#### Clone the repository
```
git clone git@github.com:Genome-Function-Initiative-Oxford/ESVAR.git
cd ESVAR
```

 
#### Create Anaconda environment
Activate the conda 'base' environment (if not active): 
```
conda activate base
```

There are two ways to create the conda env __esvar__ environment:
1) Using mamba (if ```Mamba``` is installed), and follow the on screen instructions:
```
mamba env create --file=envs/esvar.yml
```
2) Using conda, and follow the on screen instructions.
```
conda env create --file=envs/esvar.yml
```

#### Activate the environment
Now, the __esvar__ environment is created it needs to be activated: 
```
conda activate esvar
```
You can then use ESVAR using this environment, enjoy!

### Environment installation note
ESVAR has been successfully tested for the following operating systems: Ubuntu, CentOS, macOS (Intel CPU), and Windows. Unfortunately, it is not possible to install on macOS with M CPUs at the moment. 
For any error in the installation step, please open an [issue](https://github.com/Genome-Function-Initiative-Oxford/ESVAR/issues) so we can give a general solution for users.

### Reproducibility :repeat:
If required for publication, package versions within the environment can be exported as follows:
```
conda activate esvar
conda env export > esvar_environment_versions.yml
```

***

## How to use ESVAR

There are 2 ways to use ESVAR:

#### 1) from a code editor (e.g., VS Code) or a web-based interactive computing (e.g., Jupyter Notebook):
```
import pandas as pd
import sys
sys.path.append('esvar/')
from esvar import esvar

genetic  = "genetics_test/Bcell_GO_0035456.tsv"
es       = esvar(genetic=genetic, output="Test", collection_name="cad")
df_esvar = es.calculate_enrichment_score(less100=False, greater25k=False)
df_cover = es.get_coverage()
```

You can find some helpful parameter information using:
```
es = esvar(genetic=genetic)
help(es)
```

#### 2) from a terminal:
- for a single data collection ```python esvar.py -g genetics_test/Bcell_GO_0035456.tsv -cn super_pbmc -o test_api/single_run```
- for all data collection ```python esvar.py -g genetics_test/Bcell_GO_0035456.tsv -cn super_pbmc -o test_api/all_run -ra True```

Here some usage information
```
usage: esvar.py \
	-g/--genetic [Required. Path to genetic file.] : str \
	-nof/--number_of_folds [Number of folds to create backgound using the 1000genomes.] : int \
	-o/--output [Directory where to save the scores.]  : str \
	-gb/--genome [Genome to use, available 'hg19' and 'hg38'.] : str \
	-s/--seed [Seed for reproducibility, shuffle 1000genomes excluded.] : int \
	-cn/--collection_name [Data collection name, available 'cad', 'calderon', 'catlas_fetal', 'catlas_adult', 'erythoid_d7_d10_d13_d17', 'h1_hescs', 'immune_cell', 'ludwig2019', 'mpal', 'pancreatic_pbmc', and 'super_pbmc'.] : str \
	-incp/--in_house_collection_path [In house data collection path (<path-to-directory>/<collection-name>).] : str \
	-l100/--less100 [Boolean variable to force the software to run also with less than 100 variants per file.] : bool \
	-g25k/--greater25k [Boolean variable to check if you want to run ESVAR on more than 25k variants.] : bool \
	-gc/--get_coverage [Run only coverage calculation. Enrichment score will be ignored.] : bool \
	-ra/--run_all [Run ESVAR for all data collection. If -cn or -incp are set, they will be ignored.] : bool

```

***

### Output ESVAR result

ESVAR will create the following tree-like format example result.

```
ESVAR_result
├── coverage.csv
├── folds
│	├── ALL_1000_genomes.variants.hg38.bed
│	├── round_1
│	│	├── SUB_1.bed
│	│	├── SUB_k.bed
│	│	└── SUB_nof.bed
│	└── roundk
│		├── SUB_1.bed
│		├── SUB_k.bed
│		└── SUB_nof.bed
├── rounds
│	├── statistics_intermediate_round_1.csv
│	└── statistics_intermediate_round_k.csv
└── statistics_ESVAR.csv
```


***

### Genetic variant file format

Any genetic variant file provided to ESVAR has to be tab (\t) separated and must contain at least 3 columns:
1) "CHR_ID" : chromosome in the follwoing format *chrV*
2) "CHR_POS" : chromosome position
3) "SNPS" : ID of the variant (e.g., rs#####, or chrV-pos-ref-alt)

For example like:
```
CHR_ID  CHR_POS SNPS
chr10   801748  rs60692108
chr10   823912  rs74876360
chr10   840700  chr10-840700-A-C
```

***

### In-house data format

All the data has to be stored in a folder called with your collection name __<collection-name>__, following the tree-like format example below.

```
└── <path-to-directory>/<collection-name>
 ├── bigwigs/cell-type-name*.bw
 ├── peaks/cell-type-name*.bed
 └── <collection-name>_info.csv 
```

***

### Pipeline updates :construction:
If any changes are made to ESVAR, it is possible to update the repository by entering the main folder and pulling the update using:
   ```
   # Enter the main folder
   cd ESVAR

   # Pull updates
   git pull           
   ```
Alternatively, remove the cloned repository and then re-clone the repository as described above.   
Warning: use rm carefully!

```
rm -rf ESVAR
``` 
<hr>

### :warning: Warning for University of Oxford CCB users :warning:
When using this repository, use the default terminal and __do not__ load any module in the server (if logged-in).

***


### Contact us
If you have any suggestions, spot any errors, or have any questions regarding the pipelines, please do no hesitate to contact us anytime.   

:email: &emsp; [<simone.riva@imm.ox.ac.uk>](simone.riva@imm.ox.ac.uk)
