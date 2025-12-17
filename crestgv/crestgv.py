import warnings
warnings.filterwarnings('ignore')
import glob, sys, os, subprocess, shutil, pyBigWig, re, random
# import pybedtools
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import multiprocessing as mp
import seaborn as sns
import pandas as pd
import numpy as np
import multiprocessing as mp
import pyranges as pr
from scipy import stats
from tqdm import tqdm


class crestgv():
	
	def __init__(self, genetic=None, output="output", genome="hg38", min_number_genetics=500, seed=42, collection_name="", in_house_collection_path=""):

		"""\
			__init__.

			Parameters
			----------
			genetic : str
				Path to genetic file.
			number_of_folds : integer
				Number of folds to create backgound using the 1000genomes.
			output : str
				Directory where to save the scores.
			genome : str
				Genome to use, available 'hg19' and 'hg38'.
			seed : int
				Seed for reproducibility, shuffle 1000genomes excluded.
			collection_name : str
				Data collection name, available 'cad', 'calderon', 'catlas_fetal', 'catlas_adult', 'erythoid_d7_d10_d13_d17', 'greenleaf_brain', 'h1_hescs', 'immune_cell', 'ludwig2019', 'mpal', 'pancreatic_pbmc', and 'super_pbmc'.
			in_house_collection_path : str
				In house data collection path (<path-to-directory>/<collection-name>) The directory has to be structured as: └── <path-to-directory>/<collection-name>
																																├── bigwigs/cell-type-name*.bw
																																├── peaks/cell-type-name*.bed
																																└── <collection-name>_info.csv # list of name present in bigwig and peak folders without extensionts (.bw and .bed)

		"""

		self.seed = seed
		if (self.seed != None) and (isinstance(self.seed , int)):
			random.seed(self.seed)
			np.random.seed(self.seed)

		if os.path.isfile(genetic):
			self.genetic_file = genetic
		else:
			sys.exit("Provide correct path to genetic file.")

		self.genetic_df = pd.DataFrame()
		self.min_number_genetics = min_number_genetics

		self.output = output
		self.tmp = self.output+'/tmp'
		
		self.genome = genome

		if self.genome == "hg38":
			self.mappable_bp = 3049315783 #https://genomewiki.ucsc.edu/index.php?title=Hg38_27-way_Genome_size_statistics
		elif self.genome == "hg19":
		# 	self.mappable_bp = 2897310462 #https://genomewiki.ucsc.edu/index.php?title=Hg19_100way_Genome_size_statistics
			sys.exit("Genome 'hg19' not functional at the moment. Please use genome 'hg38'.")
		else:
			sys.exit("Select genome between 'hg19' and 'hg38'.")
		
# 		self.URL = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/CREST-GV_collection"
		self.URL = "/project/Wellcome_Discovery/datashare/CREST-GV_collection"

		self.collection_name_dict = {'cad'                     : 'CAD',
									 'calderon'                : 'calderon',
									 'catlas_fetal'            : 'catlas_fetal',
									 'catlas_adult'            : 'catlas_adult',
									 'erythoid_d7_d10_d13_d17' : 'Days7_10_13_17',
									 'greenleaf_brain'         : 'greenleaf_brain',
									 'h1_hescs'                : 'H1_hESCs',
									 'immune_cell'             : 'immune_cell',
									 'ludwig2019'              : 'ludwig2019',
									 'mpal'                    : 'MPAL',
									 'pancreatic_pbmc'         : 'pancreatic_pbmc',
									 'super_pbmc'              : 'super_PBMC'
									 }

		if (not in_house_collection_path) & (not collection_name):
			sys.exit("Please specify a collection name or a in house collection path.")
		elif not in_house_collection_path:
			if collection_name == "":
				sys.exit("Please specify a collection name from: [%s]"%(', '.join(list(self.collection_name_dict.keys()))))
			if collection_name not in list(self.collection_name_dict.keys()):
				sys.exit("Wrong 'collection_name' selected, please specify a collection name from: [%s]"%(', '.join(list(self.collection_name_dict.keys()))))
			self.collection_name = collection_name
			self.collectionHouseBool = False
		else:
			self.in_house_collection_path = in_house_collection_path
			self.collection_name = in_house_collection_path.split("/")[-1]
			self.collectionHouseBool = True

		# Read in chunks
		if self.collectionHouseBool:
			print("Loading 1000genome snps...")
			chunk_size = 100000  # Adjust based on your memory
			chunks = []

			for chunk in pd.read_csv(self.URL+os.sep+"1000genomes/ALL_1000_genomes.variants.hg38.bed", sep='\t',
									 chunksize=chunk_size, 
									 low_memory=False,  # Handle mixed data types,
									 header=None,
									 dtype={
										 'Chromosome': 'category',  # Use category for low memory
										 'Start': 'int32',    # Reduce precision
										 'End': 'int32',    # Reduce precision
										 'ID': 'category',    # Reduce precision
									 }):
				chunks.append(chunk[[0,1,2]])

			# Concatenate chunks
			snp_data = pd.concat(chunks, ignore_index=True)
			snp_data.columns = ["Chromosome", "Start", "End"]
			self.n_snp_data = snp_data.shape[0]
			self.snps = pr.PyRanges(snp_data)


# 	def __check_genetic_format(self, df_genetics):

# 		"""\
# 			Check genetic format.

# 			Parameters
# 			----------
# 			df_genetics : pandas.DataFrame
# 				DataFrame of initial genetics.

# 			Returns
# 			-------
# 			df_genetics : pandas.DataFrame
# 				DataFrame of filtered genetics.
# 		"""

# 		if not os.path.exists(self.genetic_file):
# 			sys.exit("Error, gentic file does not exist.")

# 		print("Removing NaN rows from loaded file ...")
# 		try:
# 			df_genetics = df_genetics[["CHR_ID", "CHR_POS", "SNPS"]]
# 		except:
# 			sys.exit("Error, gentic file must contain at least 3 columns:  \
# 						collection_name_dict\t'CHR_ID': int or fload value, \
# 						collection_name_dict\t'CHR_POS': int or fload value, \
# 						collection_name_dict\t'SNPS': string value.")

# 		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(str)
# 		df_genetics = df_genetics[~df_genetics['CHR_POS'].str.contains(";")] # check when happens
# 		df_genetics = df_genetics[~df_genetics['CHR_POS'].str.contains("x")] # check when happens
		
# 		df_genetics = df_genetics[np.isfinite(df_genetics['CHR_POS'])]
		
# 		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(float).astype(int)
# 		df_genetics = df_genetics.dropna()

# 		if df_genetics.shape[-1] == 1:
# 			sys.exit("Error, gentic file must be tab delimited.")
		
# 		if df_genetics.shape[-1] < 3:
# 			sys.exit("Error, gentic file must contain at least 3 columns:  \
# 					collection_name_dict\t'CHR_ID': int or fload value, \
# 					collection_name_dict\t'CHR_POS': int or fload value, \
# 					collection_name_dict\t'SNPS': string value.")

# 		if df_genetics.shape[-1] >= 3:
# 			if ("CHR_ID" not in df_genetics.columns) | ("CHR_POS" not in df_genetics.columns) | ("SNPS" not in df_genetics.columns):
# 				sys.exit("Error, gentic file must contain at least 3 columns:  \
# 						collection_name_dict\t'CHR_ID': int or fload value, \
# 						collection_name_dict\t'CHR_POS': int or fload value, \
# 						collection_name_dict\t'SNPS': string value.")

# 		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(int)
# 		if df_genetics.dtypes["CHR_POS"] not in ['int32', 'int64', 'float32', 'float64']:
# 			sys.exit("Error, column 'CHR_POS' does not contain all numeric values.")
# 		if df_genetics.dtypes["SNPS"] not in ['str', 'object']:
# 			sys.exit("Error, column 'SNPS' does not contain all string/object values.")  
		
# 		return df_genetics

	def __check_genetic_format(self, df_genetics):
		"""
		Check genetic format.

		Parameters
		----------
		df_genetics : pandas.DataFrame
			DataFrame of initial genetics.

		Returns
		-------
		df_genetics : pandas.DataFrame
			DataFrame of filtered genetics.
		"""

		if not os.path.exists(self.genetic_file):
			sys.exit("Error, genetic file does not exist.")

		print("Checking and cleaning genetic file...")

		# Ensure required columns exist
		required_columns = ["CHR_ID", "CHR_POS", "SNPS"]
		try:
			df_genetics = df_genetics[required_columns]
		except KeyError as e:
			sys.exit(f"Error: Missing required column {e}. Genetic file must contain: " + 
					 "'CHR_ID' (int/float), 'CHR_POS' (int/float), 'SNPS' (string).")

		# Convert CHR_POS to string first to handle potential problematic entries
		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(str)
                
		# Remove rows with problematic entries
		df_genetics = df_genetics[
			~df_genetics['CHR_POS'].str.contains(';') &  # Remove entries with semicolons
			~df_genetics['CHR_POS'].str.contains('x', case=False)  # Remove entries with 'x'
		]
        
		# Convert CHR_POS to numeric, coercing errors to NaN
		df_genetics['CHR_POS'] = pd.to_numeric(df_genetics['CHR_POS'], errors='coerce')
                
		# Remove NaN values
		df_genetics = df_genetics.dropna(subset=['CHR_POS', 'CHR_ID', 'SNPS'])
                
		# Convert CHR_POS to integer
		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(int)
                
# 		# Additional type checking
# 		try:
# 			# Ensure CHR_ID is numeric
# 			df_genetics['CHR_ID'] = pd.to_numeric(df_genetics['CHR_ID'], errors='raise')

# 			# Ensure SNPS is string
# 			df_genetics['SNPS'] = df_genetics['SNPS'].astype(str)
# 		except ValueError as e:
# 			sys.exit(f"Error in data type conversion: {e}")

		# Final validation
		if df_genetics.empty:
			sys.exit("Error: No valid data remains after filtering.")

		print(f"Processed genetic file. Retained {len(df_genetics)} rows.")

		return df_genetics


	def __load_genetic(self):

		"""\
			Load genetic.

			Parameters
			----------
			lessNG : Boolean
				Boolean variable to force the software to run also with less than 100 variants per file.
			greater25k : Boolean
				Boolean variable to check if you want to run CREST-GV on more than 25k variants.

			Returns
			-------
			df_genetics : pandas.DataFrame
				DataFrame of filtered genetics.
		"""

		df_genetics = pd.read_csv(self.genetic_file, sep="\t", names=['SNPS', 'CHR_ID','CHR_POS'])#,'CHR_POS+1','SNPSa','R','A','X','Y','Z'])
		df_genetics = self.__check_genetic_format(df_genetics)
		df_genetics = df_genetics.dropna(subset=['CHR_ID'], axis=0)
		df_genetics["CHR_ID"] = df_genetics["CHR_ID"].astype(str)
		df_genetics["CHR_POS"] = df_genetics["CHR_POS"].astype(int)
		df_genetics["SNPS"] = df_genetics["SNPS"].astype(str)
		df_genetics["CHR_ID"] = "chr"+df_genetics["CHR_ID"].str.replace(".0", "", regex=False)
		df_genetics["CHR_ID"] = df_genetics["CHR_ID"].str.replace("chrchr", "chr", regex=False)
		df_genetics["CHR_POS+1"] = df_genetics["CHR_POS"]+1
		df_genetics = df_genetics.sort_values(["CHR_ID", "CHR_POS"])
		df_genetics = df_genetics.reset_index(drop=True)
		df_genetics = df_genetics[['CHR_ID', 'CHR_POS', 'CHR_POS+1', 'SNPS']]
		df_genetics = df_genetics.drop_duplicates()

		print("Total number of used variants in CREST-GV: %s"%df_genetics.shape[0])
		if df_genetics.shape[0]<self.min_number_genetics:
			sys.exit("Genetics provided after quality control contains less than %s entry variants.") 
		self.genetic_df = df_genetics
		return df_genetics


	def __loading_info(self):

		"""\
			Load cell type information per collection type.

			Returns
			-------
			info : pandas.DataFrame
				DataFrame of cell type names per data collection.
		"""
		if self.collectionHouseBool:
			file = '%s/%s_info.csv'%(self.in_house_collection_path, self.collection_name)
			info = pd.read_csv(file, sep='\t', names=['cellType','scCounts'])[['cellType']] # info has to contain at least the list of cell type name as per bigwig and bed name files
		else:
			file = '%s/%s/%s_info.csv'%(self.URL, self.collection_name_dict[self.collection_name], self.collection_name_dict[self.collection_name])
			info = pd.read_csv(file, sep='\t', names=['cellType','scCounts'])[['cellType']]
		return info


	def __loading_bigwigs_and_beds(self, info):

		"""\
			Load bigwig and bed files per data collection.

			Parameters
			----------
			info : pandas.DataFrame
				DataFrame of cell type names per data collection.

			Returns
			-------
			bigwigs : list
				List of bigwig files per datacollection ordered per cell type within data collection.
			beds : list
				List of beds files per datacollection ordered per cell type within data collection.
		"""

		bigwigs, beds = [], []
		for ct in info['cellType'].tolist():
			if self.collectionHouseBool:
				bigwigs.append("%s/bigwigs/%s.bw"%(self.in_house_collection_path, ct))
				beds.append("%s/peaks/%s.bed"%(self.in_house_collection_path, ct))
			else:
				bigwigs.append("%s/%s/bigwigs/%s.bw"%(self.URL, self.collection_name_dict[self.collection_name], ct))
				beds.append("%s/%s/peaks/%s_L-tron.bed"%(self.URL, self.collection_name_dict[self.collection_name], ct))
		return bigwigs, beds


	def __loading_collection_data(self):

		"""\
			Load selected data collection.

			Returns
			-------
			info : pandas.DataFrame
				DataFrame of cell type names per data collection.
			bigwigs : list
				List of bigwig files per datacollection ordered per cell type within data collection.
			beds : list
				List of beds files per datacollection ordered per cell type within data collection.
		"""

		info = self.__loading_info()
		bigwigs, beds = self.__loading_bigwigs_and_beds(info)
		return info, bigwigs, beds


	def __prepare_data(self, beds, df_genetics):

		"""\
			Check genetic format.

			Parameters
			----------
			beds : list
				List of beds files per datacollection ordered per cell type within data collection.
			df_genetics : pandas.DataFrame
				DataFrame of filtered genetics.

			Returns
			-------
			df_data : pandas.DataFrame
				DataFrame of initialised information for encrichment score calculation.
		"""
		
		df_data = pd.DataFrame()
# 		peak_area, ps, celltypes = [], [], []
# 		print("Calculating (1) total number of base-pairs within peaks and (2) total number of base-pairs within peaks divided by uniquely mappable base-pairs ...")
# 		for bed in tqdm(beds):
# 			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
# 			celltype = bed.split("/")[-1].replace("_L-tron.bed","")
# 			celltypes.append(celltype)
# 			# total number of base-pairs within peaks
# 			tot_bp_within_peaks = np.sum(tmp['end']-tmp['start'])
# 			peak_area.append(tot_bp_within_peaks)
# 			# total number of base-pairs within peaks divided by uniquely mappable base-pairs
# 			p = tot_bp_within_peaks/self.mappable_bp
# 			ps.append(p)

# 		df_data["Peak_area"] = peak_area
# 		df_data["p_succes"]  = ps
# 		df_data.index		 = celltypes

### new only for new collections
		if self.collectionHouseBool:
			peak_area, ps, celltypes = [], [], []
			print("Calculating (1) total number of snps from 1000genomes within peaks ...")
			for bed in tqdm(beds):
				region_data = pd.read_csv(bed, sep="\t", names=["Chromosome", "Start", "End"])
				celltype = bed.split("/")[-1].replace("_L-tron.bed","")
				celltypes.append(celltype)
				regions = pr.PyRanges(region_data)
				intersection = self.snps.intersect(regions)
				num_intersections = len(intersection)
				p = num_intersections/self.n_snp_data
				ps.append(p)

			df_data["n_snps"]    = int(self.n_snp_data)
			df_data["p_succes"]  = ps
			df_data.index        = celltypes
		else:
			info_ps = pd.read_csv('%s/%s/%s_ps.csv'%(self.URL, self.collection_name_dict[self.collection_name], self.collection_name_dict[self.collection_name]), index_col=0, sep='\t')
			df_data["n_snps"] = info_ps['n_snps']
			df_data["p_succes"] = info_ps['p_succes']
### new


# 		if not os.path.exists(self.tmp):
# 			os.makedirs(self.tmp)
# 		pybedtools.set_tempdir(self.tmp)

# 		xs = []
# 		df_bed = pybedtools.BedTool.from_dataframe(df_genetics)
# 		print("Intersecting genetic with peak regions ...")
# 		for bed in tqdm(beds):
# 			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
# 			tmp_bed = pybedtools.BedTool.from_dataframe(tmp)
# 			intersect_bed = df_bed.intersect(tmp_bed)
# 			intersect_bed = intersect_bed.to_dataframe()
# 			if intersect_bed.empty:
# 				xs.append(0.0)
# 			else:
# 				xs.append(intersect_bed.shape[0])
# 		df_data["GWAS_init"] = xs
# 		return df_data

		df_genetics.columns = ['Chromosome', 'Start', 'End', 'SNPS']   
		py_df_genetics = pr.PyRanges(df_genetics[['Chromosome', 'Start', 'End']])
		xs = []
		for bed in tqdm(beds):
			region_data = pd.read_csv(bed, sep="\t", names=['Chromosome', 'Start', 'End']  )
			regions = pr.PyRanges(region_data)
			intersection = py_df_genetics.intersect(regions)
			num_intersections = len(intersection)
			xs.append(num_intersections)
		df_data["GWAS_init"] = xs
        
		return df_data


	def __add_statistics(self, df_data, number_of_genetic):

		"""\
			Add statisticts to .

			Parameters
			----------
			df_data : pandas.DataFrame
				DataFrame of initialised and background information for encrichment score calculation.
			number_of_genetic : int
				Number of variant to extract from 1000genomes.
				
			Returns
			-------
			df_data : pandas.DataFrame
				Encrichment score DataFrame.
		"""

		print("Calculating statistics ...")

		df_data["pmf"] = stats.binom.pmf(df_data["GWAS_init"], number_of_genetic, df_data["p_succes"])
		df_data["sf"] = stats.binom.sf(df_data["GWAS_init"]-1, number_of_genetic, df_data["p_succes"])
		
		df_data["check_np"] = number_of_genetic*df_data["p_succes"]
		df_data["check_n(1-p)"] = number_of_genetic*(1-df_data["p_succes"])
		
		df_data["p_hat"] = df_data["GWAS_init"]/number_of_genetic
		
		df_data["SE"] = np.sqrt((df_data["p_succes"]*(1-df_data["p_succes"]))/number_of_genetic)
		
		df_data["Z"] = (df_data["p_hat"]-df_data["p_succes"])/(df_data["SE"])
		df_data["pdf"] = stats.norm.sf(df_data["Z"])
        
		df_data["Z_norm"] = (df_data["Z"]-np.mean(df_data["Z"]))/np.std(df_data["Z"])
		df_data["pdf_cgv_norm"] = stats.norm.sf(df_data["Z_norm"])

		df_data["CREST-GV"] = df_data["Z"]
		df_data["CREST-GV_pval"] = df_data["pdf"]
# 		df_data["CREST-GV_pval_adjust"] = stats.false_discovery_control(df_data["CREST-GV_pval"], method='by')
        
		df_data["CREST-GV_norm"] = df_data["Z_norm"]
		df_data["CREST-GV_pval_norm"] = df_data["pdf_cgv_norm"]
# 		df_data["CREST-GV_pval_norm_adjust"] = stats.false_discovery_control(df_data["CREST-GV_pval_norm"], method='by')


		return df_data


# 	def __clean_tmp(self):

# 		"""\
# 			Delete temporary files.
# 		"""
		
# 		print("Cleaning temporary files ...")
# 		shutil.rmtree(self.tmp)



	def calculate_enrichment_score(self):

		"""\
			Calculate enrichment score for provided genetics per selected data collection.

			Parameters
			----------
			lessNG : Boolean
				Boolean variable to force the software to run also with less than 100 variants per file.
			greater25k : Boolean
				Boolean variable to check if you want to run CREST-GV on more than 25k variants.

			Returns
			-------
			dfs_collection : pandas.DataFrame
				Final enrichment score DataFrame for provided genetics and selected data collection.
		"""

		if self.genetic_file is None:
			sys.exit("Error, missing genetic file.")

		df_genetics = self.__load_genetic()

		info, _, beds = self.__loading_collection_data()
        
		df_collection = self.__prepare_data(beds, df_genetics)

		df_collection = self.__add_statistics(df_collection, df_genetics.shape[0])

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_collection.to_csv(self.output+os.sep+"%s_statistics_CREST-GV.csv"%self.collection_name, sep="\t")

# 		self.__clean_tmp()
		print("Processing genetics finished.")

		return df_collection


	def get_coverage(self):

		"""\
			Get coverage for the genetic variants within per selected data collection.

			Returns
			-------
			dfs_collection : pandas.DataFrame
				Final enrichment score DataFrame for provided genetics and selected data collection.
		"""

		self.genetic_df = self.__load_genetic()
		self.genetic_df['CHR_POS+1'] = self.genetic_df['CHR_POS']+1

		info, bigwigs, _ = self.__loading_collection_data()

		print("Getting coverage ...")

		for ct in tqdm(info['cellType'].tolist()):
			values_list = []
			if self.collectionHouseBool:
				bigwig = "%s/bigwigs/%s.bw"%(self.in_house_collection_path, ct)
			else:
				bigwig = "%s/%s/bigwigs/%s.bw"%(self.URL, self.collection_name_dict[self.collection_name], ct)

			bw = pyBigWig.open(bigwig)
			for c, s, e in zip(self.genetic_df['CHR_ID'].tolist(), self.genetic_df['CHR_POS'].tolist(), self.genetic_df['CHR_POS+1'].tolist()):
				values_list.append(bw.values(c, s, e, numpy=True).item())
			bw.close()

			self.genetic_df[ct] = values_list

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		self.genetic_df.to_csv(self.output+os.sep+"coverage.csv", sep="\t")

		return self.genetic_df