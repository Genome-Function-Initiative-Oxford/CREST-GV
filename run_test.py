import sys, glob
import pandas as pd

sys.path.append('esvar/')
from esvar import esvar

collections = [ 'cad', 
				'calderon', 
				'catlas_fetal', 
				'catlas_adult', 
				'erythoid_d7_d10_d13_d17', 
				'h1_hescs', 
				'immune_cell', 
				'ludwig2019', 
				'mpal', 
				'pancreatic_pbmc', 
				'super_pbmc'
			   ]

genetic = "/project/Wellcome_Discovery/sriva/Git/ESVAR/genetics_test/purl.obolibrary.org_obo_GO_0035456.tsv"

for collection in collections:
	es = esvar(genetic=genetic, output="test/%s"%collection, collection_name=collection)
	df = es.calculate_enrichment_score(less100=False)