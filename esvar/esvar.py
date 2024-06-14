import warnings
warnings.filterwarnings('ignore')

import glob, pybedtools, sys, os, subprocess, shutil, pyBigWig, re, random
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import multiprocessing as mp
import seaborn as sns
import pandas as pd
import numpy as np

from multiprocessing import Pool
from scipy import stats
from tqdm import tqdm


class esvar():
	
	def __init__(self, genetic=None, number_of_folds=5, output="output", genome="hg38", seed=42, collection_name=None):

		self.seed = seed
		if (self.seed != None) and (isinstance(self.seed , int)):
			random.seed(self.seed)
			np.random.seed(self.seed)

		self.df_genetics = genetic

		self.output = output
		self.tmp = self.output+'/tmp'
		self.folds = self.output+'/folds'
		
		self.number_of_folds = number_of_folds
		self.genome = genome

		if self.genome == "hg38":
			self.mappable_bp = 3049315783 #https://genomewiki.ucsc.edu/index.php?title=Hg38_27-way_Genome_size_statistics
			self.background = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/ESVAR_collection/1000genomes/ALL_1000_genomes.variants.hg38.bed"
		# elif self.genome == "hg19":
		# 	self.mappable_bp = 2897310462 #https://genomewiki.ucsc.edu/index.php?title=Hg19_100way_Genome_size_statistics
		# 	self.background = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/ESVAR_collection/1000genomes/ALL_1000_genomes.variants.hg19.bed"
		else:
			sys.exit("Select genome between 'hg19' and 'hg38'.")

		self.df_collection = pd.DataFrame()
		
		self.URL = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/ESVAR_collection"

		self.collection_name_dict = {'cad'                     : 'CAD',
									 'calderon'                : 'calderon',
									 'catlas_fetal'            : 'catlas_fetal',
									 'catlas_adult'            : 'catlas_adult',
									 'erythoid_d7_d10_d13_d17' : 'Days7_10_13_17',
									 'h1_hescs'                : 'H1_hESCs',
									 'immune_cell'             : 'immune_cell',
									 'ludwig2019'              : 'ludwig2019',
									 'mpal'                    : 'MPAL_lowGr',
									 'pancreatic_pbmc'         : 'pancreatic_pbmc',
									 'super_pbmc'              : 'super_PBMC'
									 }

		if collection_name == None:
			sys.exit("Please specify a collection name from: [%s]"%(', '.join(list(self.collection_name_dict.keys()))))
		if collection_name not in list(self.collection_name_dict.keys()):
			sys.exit("Wrong 'collection_name' selected, please specify a collection name from: [%s]"%(', '.join(list(self.collection_name_dict.keys()))))
		self.collection_name = collection_name


	def __check_genetic_format(self, df_genetics):

		if not os.path.exists(self.df_genetics):
			sys.exit("Error, gentic file does not exist.")

		print("Removing NaN rows from loaded file ...")
		try:
			df_genetics = df_genetics[["CHR_ID", "CHR_POS", "SNPS"]]
		except:
			sys.exit("Error, gentic file must contain at least 3 columns:  \
						collection_name_dict\t'CHR_ID': int or fload value, \
						collection_name_dict\t'CHR_POS': int or fload value, \
						collection_name_dict\t'SNPS': string value.")

		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(str)
		df_genetics = df_genetics[~df_genetics['CHR_POS'].str.contains(";")] # check when happens
		df_genetics = df_genetics[~df_genetics['CHR_POS'].str.contains("x")] # check when happens
		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(int)
		df_genetics = df_genetics.dropna()

		if df_genetics.shape[-1] == 1:
			sys.exit("Error, gentic file must be tab delimited.")
		
		if df_genetics.shape[-1] < 3:
			sys.exit("Error, gentic file must contain at least 3 columns:  \
					collection_name_dict\t'CHR_ID': int or fload value, \
					collection_name_dict\t'CHR_POS': int or fload value, \
					collection_name_dict\t'SNPS': string value.")

		if df_genetics.shape[-1] >= 3:
			if ("CHR_ID" not in df_genetics.columns) | ("CHR_POS" not in df_genetics.columns) | ("SNPS" not in df_genetics.columns):
				sys.exit("Error, gentic file must contain at least 3 columns:  \
						collection_name_dict\t'CHR_ID': int or fload value, \
						collection_name_dict\t'CHR_POS': int or fload value, \
						collection_name_dict\t'SNPS': string value.")

		df_genetics['CHR_POS'] = df_genetics['CHR_POS'].astype(int)
		if df_genetics.dtypes["CHR_POS"] not in ['int32', 'int64', 'float32', 'float64']:
			sys.exit("Error, column 'CHR_POS' does not contain all numeric values.")
		if df_genetics.dtypes["SNPS"] not in ['str', 'object']:
			sys.exit("Error, column 'SNPS' does not contain all string/object values.")

		return df_genetics


	def __load_genetic(self, less100=False):
		df_genetics = pd.read_csv(self.df_genetics, sep="\t")
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
		
		# ### to remove - start
		# if df_genetics.shape[0]>100:
		# 		df_genetics = df_genetics.sample(n=100, frac=None, replace=False, weights=None, random_state=42, axis=0)
		# 		df_genetics = df_genetics.reset_index(drop=True)
		# ### to remove - end

		if df_genetics.shape[0]<100:
			if less100 == False:
				sys.exit("Genetics provided after quality control contains less than the minimum number (100 variants) of entries.\nIf you want to carry on anyway with it, please set 'less100=True'.")
		return df_genetics


	def __loading_info(self):
		file = '%s/%s/%s_info.csv'%(self.URL, self.collection_name_dict[self.collection_name], self.collection_name_dict[self.collection_name])
		info = pd.read_csv(file, sep='\t', names=['cellType','scCounts'])[['cellType']]
		return info


	def __loading_bigwigs_and_beds(self, info):
		bigwigs, beds = [], []
		for ct in info['cellType'].tolist():
			bigwigs.append("%s/%s/bigwigs/%s.bw"%(self.URL, self.collection_name_dict[self.collection_name], ct))		
			beds.append("%s/%s/peaks/%s_L-tron.bed"%(self.URL, self.collection_name_dict[self.collection_name], ct))
		return bigwigs, beds


	def __loading_collection_data(self):
		info = self.__loading_info()
		bigwigs, beds = self.__loading_bigwigs_and_beds(info)
		return info, beds, bigwigs


	def __prepare_data(self, beds):
		df_data = pd.DataFrame()
		peak_area, ps, celltypes = [], [], []
		print("Calculating (1) total number of base-pairs within peaks and (2) total number of base-pairs within peaks divided by uniquely mappable base-pairs ...")
		for bed in tqdm(beds):
			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
			celltype = bed.split("/")[-1].replace("_L-tron.bed","")
			celltypes.append(celltype)
			# total number of base-pairs within peaks
			tot_bp_within_peaks = np.sum(tmp['end']-tmp['start'])
			peak_area.append(tot_bp_within_peaks)
			# total number of base-pairs within peaks divided by uniquely mappable base-pairs
			p = tot_bp_within_peaks/self.mappable_bp
			ps.append(p)

		df_data["Peak_area"] = peak_area
		df_data["p_succes"]  = ps
		df_data.index		= celltypes

		if not os.path.exists(self.tmp):
			os.makedirs(self.tmp)
		pybedtools.set_tempdir(self.tmp)

		xs = []
		df_bed = pybedtools.BedTool.from_dataframe(self.df_genetics)
		print("Intersecting genetic with peak regions ...")
		for bed in tqdm(beds):
			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
			tmp_bed = pybedtools.BedTool.from_dataframe(tmp)
			intersect_bed = df_bed.intersect(tmp_bed)
			intersect_bed = intersect_bed.to_dataframe()
			if intersect_bed.empty:
				xs.append(0.0)
			else:
				xs.append(intersect_bed.shape[0])
		df_data["GWAS_init"] = xs
		return df_data


	def __download_background(self):
		if not os.path.exists(self.folds):
			os.makedirs(self.folds)
		if not os.path.exists(self.folds+os.sep+"ALL_1000_genomes.variants.%s.bed"%self.genome):
			print("Downloading %s ALL 1000 genomes background..."%self.genome)
			_ = subprocess.run('wget -P %s/ %s'%(self.folds, self.background), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


	def shuf(self, x):
			return subprocess.run(x, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


	def __shuffle_background(self, number_of_genetic):

		print("For each fold, parallelised shuffling background ...")

		commands = []
		for f in range(1, self.number_of_folds+1):
			commands.append("shuf -n %s %s | sort -k1,1 -k2,2n > %s/SUB%s.bed"%(number_of_genetic, self.folds+os.sep+"ALL_1000_genomes.variants.%s.bed"%self.genome, self.folds, f))

		with Pool() as pool:
			_ = pool.map(self.shuf, commands)


	def __create_background(self, number_of_genetic): 
		self.__download_background()
		self.__shuffle_background(number_of_genetic)
		
	
	def __add_shuffled_background(self, df_data, beds): # it can be parallelised as __shuffle_background
		print("Adding shuffled background ...")
		for idx, fold in enumerate(tqdm(glob.glob("%s/SUB*.bed"%self.folds))):
			df_fold = pd.read_csv(fold, sep="\t", header=None)[[0,1,2]]
			xs_fold = []
			df_bed = pybedtools.BedTool.from_dataframe(df_fold)
			for bed in beds:
				try:
					tmp = pd.read_csv(bed, sep="\t")
					tmp_bed = pybedtools.BedTool.from_dataframe(tmp)
					intersect_bed = df_bed.intersect(tmp_bed)
					intersect_bed = intersect_bed.to_dataframe()
					xs_fold.append(intersect_bed.shape[0])
				except:
					xs_fold.append(0.0)
			df_data["bg_%s"%(idx+1)] = xs_fold		
		return df_data


	def __add_statistics(self, df_data, number_of_genetic):

		print("Calculating statistics ...")

		df_data["P_gwas"] = stats.binom.pmf(df_data["GWAS_init"], number_of_genetic, df_data["p_succes"])

		for f in range(1, self.number_of_folds+1):
			df_data["P_ss%s"%f]  = stats.binom.pmf(df_data["bg_%s"%f], number_of_genetic, df_data["p_succes"])
			df_data["FOLD%s"%f] = -np.log10(df_data["P_gwas"])/-np.log10(df_data["P_ss%s"%f])
		col = df_data.loc[: , "FOLD1":"FOLD%s"%f]
		df_data['ESVAR'] = col.mean(axis=1)

		return df_data


	def __clean_tmp(self):
		print("Cleaning temporary files ...")
		shutil.rmtree(self.tmp)


	def calculate_enrichment_score(self, less100=False):
		if self.df_genetics is None:
			sys.exit("Error, missing genetic file.")
		
		df_genetics = self.__load_genetic(less100=less100)
		# self.df_genetics = df_genetics
		# number_of_genetic


		if less100:
			df_genetics_list = [df_genetics]
		else:
			df_genetics_list = []
			sub_n = 100
			df_genetics_init = df_genetics.sample(sub_n, random_state=self.seed)
			df_genetics_list.append(df_genetics_init)
			for i in range(int(df_genetics.shape[0]/sub_n)-1):
			    if i == 0:
			        df_genetics_rest = df_genetics[~df_genetics.index.isin(df_genetics_init.index)]
			    else:
			        df_genetics_rest = df_genetics_rest[~df_genetics_rest.index.isin(df_genetics_round.index)]
			    df_genetics_round = df_genetics_rest.sample(sub_n, random_state=42)
			    df_genetics_list.append(df_genetics_round)

		### to parallelise - start
		dfs_collection = []
		for idx, df_genetic in enumerate(df_genetics_list):
			self.df_genetics = df_genetic
			number_of_genetic = self.df_genetics.shape[0]
			info, beds, _ = self.__loading_collection_data()
			df_collection = self.__prepare_data(beds)
			self.__create_background(number_of_genetic)
			df_collection = self.__add_shuffled_background(df_collection, beds)
			df_collection = self.__add_statistics(df_collection, number_of_genetic)
			dfs_collection.append(df_collection[['ESVAR']])

			if not os.path.exists(self.output):
				os.makedirs(self.output)
			df_collection.to_csv(self.output+os.sep+"statistics_intermediate_round%s.csv"%(idx+1), sep="\t")

			self.__clean_tmp()
			print("Processing genetics finished.")
		### to parallelise - end

		df_collection_final = pd.concat(dfs_collection, axis=1)
		df_collection_final['ESVAR_all'] = df_collection_final.mean(axis=1)
		df_collection_final = df_collection_final[['ESVAR_all']]

		df_collection_final.to_csv(self.output+os.sep+"statistics_ESVAR.csv", sep="\t")
		self.df_collection = df_collection_final

		return df_collection
