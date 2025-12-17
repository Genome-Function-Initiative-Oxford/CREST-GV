import sys, glob
import pandas as pd

sys.path.append('crestgv/')
from crestgv import crestgv

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

genetic = "<genetic_path>"

for collection in collections:
	cgv = crestgv(genetic=genetic, output="test/%s"%collection, collection_name=collection, min_number_genetics=100)
	_ = cgv.calculate_enrichment_score()