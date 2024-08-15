# CREST-GV - Cell types Ranking and Enrichment Score for selecTive Genetic Variants

CREST-GV is a method which allows querying [our described data collection](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/blob/main/our_collection/data_collection.txt) of ~500 cell types (we keep piling more data to add to the data collection) to determine the enrichment score of a set of genetic variants. CREST-GV rely on peak properties for each cell type in the data collection leveraging [LanceOtron](https://github.com/LHentges/LanceOtron) peak caller.
CREST-GV can also query your personal (in-house) data where formatted correctly (see [In-house data format](#In-house-data-format) section for a properly formatted data structure), in this case, the user can rely on the peak caller of their own choice.
The only mandatory input for the tool is the path of a genetic file, stored following the format shown in [Genetic variant file format](#genetic-variant-file-format) section.

***

## Getting started
CREST-GV can be uesed and run using the __crestgv__ conda envrironment. Please follow the installation instruction detailed below.

### Installation instructions for conda environment

This section rely on the assumption that any distribution of ```Conda``` (e.g., [Anaconda](https://www.anaconda.com), [Miniconda](https://docs.conda.io/en/latest/miniconda.html), ...) or [```Mamba```](https://mamba.readthedocs.io/en/latest/index.html) is already installed.

#### Clone the repository
```
git clone git@github.com:Genome-Function-Initiative-Oxford/CREST-GV.git
cd CREST-GV
```

 
#### Create Anaconda environment
Activate the conda 'base' environment (if not active): 
```
conda activate base
```

There are two ways to create the conda env __crestgv__ environment:
1) Using mamba (if ```Mamba``` is installed), and follow the on screen instructions:
```
mamba env create --file=envs/crestgv.yml
```
2) Using conda, and follow the on screen instructions.
```
conda env create --file=envs/crestgv.yml
```

#### Activate the environment
Now, the __crestgv__ environment is created it needs to be activated: 
```
conda activate crestgv
```
You can then use CREST-GV using this environment, enjoy!

### Environment installation note
CREST-GV has been successfully tested for the following operating systems: Ubuntu, CentOS, macOS (Intel CPU), and Windows. Unfortunately, it is not possible to install on macOS with M CPUs at the moment. 
For any error in the installation step, please open an [issue](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/issues) so we can give a general solution for users.

### Reproducibility :repeat:
If required for publication, package versions within the environment can be exported as follows:
```
conda activate crestgv
conda env export > crestgv_environment_versions.yml
```

***

## How to use CREST-GV

There are 2 ways to use CREST-GV:

#### 1) from a code editor (e.g., VS Code) or a web-based interactive computing (e.g., Jupyter Notebook):
```
import pandas as pd
import sys
sys.path.append('crestgv/')
from crestgv import crestgv

genetic    = "genetics_test/Bcell_GO_0035456.tsv"
es         = crestgv(genetic=genetic, output="Test", collection_name="super_pbmc")
df_crestgv = es.calculate_enrichment_score(less100=False, greater25k=False)
df_cover   = es.get_coverage()
```

You can find some helpful parameter information using:
```
cgv = crestgv(genetic=genetic)
help(cgv)
```

We created a Jupyter Notebook with a heatmap plot as CREST-GV example run (see [heatmap_test.ipynb](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/blob/main/heatmap_test.ipynb)).

#### 2) from a terminal:
- for a single data collection ```python crestgv.py -g genetics_test/Bcell_GO_0035456.tsv -cn super_pbmc -o test_api/single_run```
- for all data collection ```python crestgv.py -g genetics_test/Bcell_GO_0035456.tsv -o test_api/all_run -ra True```

Here some usage information
```
usage: crestgv.py \
	-g/--genetic [Required. Path to genetic file.] : str \
	-nof/--number_of_folds [Number of folds to create backgound using the 1000genomes.] : int \
	-o/--output [Directory where to save the scores.]  : str \
	-gb/--genome [Genome to use, available 'hg19' and 'hg38'.] : str \
	-s/--seed [Seed for reproducibility, shuffle 1000genomes excluded.] : int \
	-cn/--collection_name [Data collection name, available 'cad', 'calderon', 'catlas_fetal', 'catlas_adult', 'erythoid_d7_d10_d13_d17', 'h1_hescs', 'immune_cell', 'ludwig2019', 'mpal', 'pancreatic_pbmc', and 'super_pbmc'.] : str \
	-incp/--in_house_collection_path [In house data collection path (<path-to-directory>/<collection-name>).] : str \
	-l100/--less100 [Boolean variable to force the software to run also with less than 100 variants per file.] : bool \
	-g25k/--greater25k [Boolean variable to check if you want to run CREST-GV on more than 25k variants.] : bool \
	-gc/--get_coverage [Run only coverage calculation. Enrichment score will be ignored.] : bool \
	-ra/--run_all [Run CREST-GV for all data collection. If -cn or -incp are set, they will be ignored.] : bool

```

***

### Output CREST-GV result

CREST-GV will create the following tree-like format result (see [heatmap_test folder](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/tree/main/heatmap_test) example).

```
heatmap_test
├── coverage.csv
├── folds
│   ├── round0
│   │   ├── SUB1.bed
│   │   ├── SUB2.bed
│   │   ├── SUB3.bed
│   │   ├── SUB4.bed
│   │   └── SUB5.bed
│   ├── round1
│   │   ├── SUB1.bed
│   │   ├── SUB2.bed
│   │   ├── SUB3.bed
│   │   ├── SUB4.bed
│   │   └── SUB5.bed
│   ├── round2
│   │   ├── SUB1.bed
│   │   ├── SUB2.bed
│   │   ├── SUB3.bed
│   │   ├── SUB4.bed
│   │   └── SUB5.bed
│   ├── round3
│   │   ├── SUB1.bed
│   │   ├── SUB2.bed
│   │   ├── SUB3.bed
│   │   ├── SUB4.bed
│   │   └── SUB5.bed
│   └── round4
│       ├── SUB1.bed
│       ├── SUB2.bed
│       ├── SUB3.bed
│       ├── SUB4.bed
│       └── SUB5.bed
├── rounds
│   ├── statistics_intermediate_round1.csv
│   ├── statistics_intermediate_round2.csv
│   ├── statistics_intermediate_round3.csv
│   ├── statistics_intermediate_round4.csv
│   └── statistics_intermediate_round5.csv
└── statistics_CREST-GV.csv
```


***

### Genetic variant file format

Any genetic variant file provided to CREST-GV has to be tab (\t) separated and must contain at least 3 columns:
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

All the data has to be stored in a folder called with your collection name __<collection-name>__, following the tree-like format example below (see [example_files/in-house-data folder](https://github.com/Genome-Function-Initiative-Oxford/CREST-GV/tree/main/example_files/in-house-data) example), where *in-house-data* is the *collection name*.

```
example_files/in-house-data/
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

***

### Pipeline updates :construction:
If any changes are made to CREST-GV, it is possible to update the repository by entering the main folder and pulling the update using:
   ```
   # Enter the main folder
   cd CREST-GV

   # Pull updates
   git pull           
   ```
Alternatively, remove the cloned repository and then re-clone the repository as described above.   
Warning: use rm carefully!

```
rm -rf CREST-GV
``` 
<hr>

### :warning: Warning for University of Oxford CCB users :warning:
When using this repository, use the default terminal and __do not__ load any module in the server (if logged-in).

***


### Contact us
If you have any suggestions, spot any errors, or have any questions regarding the pipelines, please do no hesitate to contact us anytime.   

:email: &emsp; [<simone.riva@imm.ox.ac.uk>](simone.riva@imm.ox.ac.uk)
