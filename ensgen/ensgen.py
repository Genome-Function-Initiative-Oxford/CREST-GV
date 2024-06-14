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

class ensgen():
	
	def __init__(self, genetic=None, tmp="tmp", folds="folds", number_of_folds=5, output="output", genome="hg38", n_subsampling=1000, seed=None):

		self.seed = seed
		if (self.seed != None) and (isinstance(self.seed , int)):
			random.seed(self.seed)
			np.random.seed(self.seed)

		self.genetic = genetic

		self.tmp = tmp
		self.folds = folds
		self.output = output
		self.number_of_folds = number_of_folds
		self.genome = genome
		self.n_subsampling = n_subsampling

		if self.genome == "hg38":
			self.mappable_bp = 3049315783 #https://genomewiki.ucsc.edu/index.php?title=Hg38_27-way_Genome_size_statistics
			self.background = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/1000genomes/ALL_1000_genomes.variants.hg38.bed"
		else:
			self.mappable_bp = 2897310462 #https://genomewiki.ucsc.edu/index.php?title=Hg19_100way_Genome_size_statistics
			self.background = "https://datashare.molbiol.ox.ac.uk/public/project/hugheslab/avocato/1000genomes/ALL_1000_genomes.variants.hg19.bed"

		### single-cell
		self.url_catlas_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/catlas/catlas_sc_info.csv"
		self.url_catlas_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/catlas/peaks"
		self.url_catlas_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/catlas/bigwigs"

		self.url_MPAL_lowGr_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/MPAL_lowGr/MPAL_lowGr_sc_info.csv"
		self.url_MPAL_lowGr_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/MPAL_lowGr/peaks"
		self.url_MPAL_lowGr_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/MPAL_lowGr/bigwigs"

		self.url_super_PBMC_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/super_PBMC/super_PBMC_sc_info.csv"
		self.url_super_PBMC_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/super_PBMC/peaks"
		self.url_super_PBMC_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/super_PBMC/bigwigs"
		self.url_super_PBMC_umap = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/super_PBMC/super_PBMC_umap.tsv"

		self.url_CAD_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/CAD/CAD_sc_info.csv"
		self.url_CAD_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/CAD/peaks"
		self.url_CAD_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/CAD/bigwigs"
		self.url_CAD_umap = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/CAD/CAD_umap.tsv"

		### bulk
		self.url_calderon_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/calderon/calderon_info.csv"
		self.url_calderon_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/calderon/peaks"
		self.url_calderon_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/calderon/bigwigs"

		self.url_ludwig2019_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/ludwig2019/ludwig2019_info.csv"
		self.url_ludwig2019_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/ludwig2019/peaks"
		self.url_ludwig2019_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/ludwig2019/bigwigs"

		self.url_Days7_10_13_17_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/Days7_10_13_17/Days7_10_13_17_info.csv"
		self.url_Days7_10_13_17_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/Days7_10_13_17/peaks"
		self.url_Days7_10_13_17_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/Days7_10_13_17/bigwigs"

		self.url_immune_cell_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/immune_cell/immune_cell_info.csv"
		self.url_immune_cell_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/immune_cell/peaks"
		self.url_immune_cell_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/immune_cell/bigwigs"

		self.url_pancreatic_pbmc_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/pancreatic_pbmc/pancreatic_pbmc_info.csv"
		self.url_pancreatic_pbmc_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/pancreatic_pbmc/peaks"
		self.url_pancreatic_pbmc_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/pancreatic_pbmc/bigwigs"

		self.url_H1_hESCs_info = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/H1_hESCs/H1_hESCs_info.csv"
		self.url_H1_hESCs_beds = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/H1_hESCs/peaks"
		self.url_H1_hESCs_bigwigs = "https://datashare.molbiol.ox.ac.uk/public/project/Wellcome_Discovery/AVOCATO/H1_hESCs/bigwigs"

		### dataframes
		self.df_catlas = pd.DataFrame()
		self.df_MPAL_lowGr = pd.DataFrame()
		self.df_super_PBMC = pd.DataFrame()
		self.df_CAD = pd.DataFrame()
		self.df_calderon = pd.DataFrame()
		self.df_ludwig2019 = pd.DataFrame()
		self.df_Days7_10_13_17 = pd.DataFrame()
		self.df_immune_cell = pd.DataFrame()
		self.df_pancreatic_pbmc = pd.DataFrame()
		self.df_H1_hESCs = pd.DataFrame()

		self.origin = ""

		self.multicoverage = pd.DataFrame()

		### custom data - start ###
		self.url_custom_info = "/project/hugheslab/shared/03_data_for_tools/avocato/custom_data/custom_data_info.csv"
		self.url_custom_beds = "/project/hugheslab/shared/03_data_for_tools/avocato/custom_data/peaks"
		self.url_custom_bigwigs = "/project/hugheslab/shared/03_data_for_tools/avocato/custom_data/bigwigs"

		self.url_multiome_ery_info = "/project/hugheslab/shared/03_data_for_tools/avocato/multiome_ery/multiome_ery_info.csv"
		self.url_multiome_ery_beds = "/project/hugheslab/shared/03_data_for_tools/avocato/multiome_ery/peaks"
		self.url_multiome_ery_bigwigs = "/project/hugheslab/shared/03_data_for_tools/avocato/multiome_ery/bigwigs"

		self.df_custom = pd.DataFrame()
		self.df_multiome_ery = pd.DataFrame()
		### custom data - end ###


	def __check_genetic_format(self):

		if not os.path.exists(self.genetic):
			sys.exit("Error, gentic file does not exist.")

		df = pd.read_csv(self.genetic, sep="\t")

		try:
			df = df[~df['CHR_POS'].str.contains(";")] # check when happens
		except:
			pass

		print("Removing NaN rows from loaded file ...")
		try:
			df = df[["CHR_ID", "CHR_POS", "SNPS"]]
		except:
			sys.exit("Error, gentic file must contain at least 3 columns:  \
						\n\t'CHR_ID': int or fload value, \
						\n\t'CHR_POS': int or fload value, \
						\n\t'SNPS': string value.")

		try:
			df = df[~df['CHR_POS'].str.contains(";")] # check when happens
		except:
			pass
		try:
			df = df[~df['CHR_POS'].str.contains("x")] # check when happens
		except:
			pass
		df = df.dropna()

		if df.shape[-1] == 1:
			sys.exit("Error, gentic file must be tab delimited.")
		
		if df.shape[-1] < 3:
			sys.exit("Error, gentic file must contain at least 3 columns:  \
					\n\t'CHR_ID': int or fload value, \
					\n\t'CHR_POS': int or fload value, \
					\n\t'SNPS': string value.")

		if df.shape[-1] >= 3:
			if ("CHR_ID" not in df.columns) | ("CHR_POS" not in df.columns) | ("SNPS" not in df.columns):
				sys.exit("Error, gentic file must contain at least 3 columns:  \
						\n\t'CHR_ID': int or fload value, \
						\n\t'CHR_POS': int or fload value, \
						\n\t'SNPS': string value.")

		try:
			df = df[~df['CHR_POS'].str.contains(";")] # check when happens
		except:
			pass
		try:
			df = df[~df['CHR_POS'].str.contains("x")] # check when happens
		except:
			pass
		df['CHR_POS'] = df['CHR_POS'].astype(int)
		# if df.dtypes["CHR_ID"] not in ['int32', 'int64', 'float32', 'float64']:
			# sys.exit("Error, column 'CHR_ID' does not contain all numeric values.")
		if df.dtypes["CHR_POS"] not in ['int32', 'int64', 'float32', 'float64']:
			sys.exit("Error, column 'CHR_POS' does not contain all numeric values.")
		if df.dtypes["SNPS"] not in ['str', 'object']:
			sys.exit("Error, column 'SNPS' does not contain all string/object values.")


	def __load_genetic(self, subsampling = False):

		self.__check_genetic_format()

		df = pd.read_csv(self.genetic, sep="\t")[["CHR_ID", "CHR_POS", "SNPS"]]

		df = df.dropna(subset=['CHR_ID'], axis=0)
		df["CHR_ID"] = df["CHR_ID"].astype(str)
		try:
			df = df[~df['CHR_POS'].str.contains(";")] # check when happens
		except:
			pass
		try:
			df = df[~df['CHR_POS'].str.contains("x")] # check when happens
		except:
			pass
		df["CHR_POS"] = df["CHR_POS"].astype(int)
		df["SNPS"] = df["SNPS"].astype(str)
		df["CHR_ID"] = "chr"+df["CHR_ID"].str.replace(".0", "", regex=False)
		df["CHR_ID"] = df["CHR_ID"].str.replace("chrchr", "chr", regex=False)
		df["CHR_POS+1"] = df["CHR_POS"]+1
		df = df.sort_values(["CHR_ID", "CHR_POS"])
		df = df.reset_index(drop=True)
		df = df[['CHR_ID', 'CHR_POS', 'CHR_POS+1', 'SNPS']]

		### random subsampling
		if subsampling:
			if df.shape[0]>self.n_subsampling:
				print("Genetic subsampled.")
				# df = df.loc[np.random.choice(df.index, 500, replace=False)]
				df = df.sample(n=self.n_subsampling, frac=None, replace=False, weights=None, random_state=42, axis=0)
				df = df.reset_index(drop=True)

		self.genetic = df

		return df.shape[0]


	def __loading_catlas_info(self, origin):

		catlas_info = pd.read_csv(self.url_catlas_info, sep='\t')
		catlas_info = catlas_info.set_index('cellType', drop=True)
		catlas_info.index.name = None
		catlas_info['scCounts'] = catlas_info['scCounts'].astype(float)
		
		if origin == 'fetal':
			catlas_info = catlas_info[catlas_info.index.str.startswith("Fetal")]
		else:
			catlas_info = catlas_info[~catlas_info.index.str.startswith("Fetal")]
		
		return catlas_info


	def __loading_catlas_beds(self, catlas_info):
		
		catlas_beds = []
		
		for ct in catlas_info.index.tolist():
			catlas_beds.append("%s/%s_L-tron.bed"%(self.url_catlas_beds, ct))
		
		return catlas_beds


	def __loading_catlas_bigwig(self, catlas_info):
		
		catlas_bigwigs = []
		
		for ct in catlas_info.index.tolist():
			catlas_bigwigs.append("%s/%s.bw"%(self.url_catlas_bigwigs, ct))
		
		return catlas_bigwigs


	def __loading_catlas_data(self, origin):

		catlas_info = self.__loading_catlas_info(origin)
		catlas_beds = self.__loading_catlas_beds(catlas_info)
		catlas_bigwigs = self.__loading_catlas_bigwig(catlas_info)

		return catlas_info, catlas_beds, catlas_bigwigs


	def __loading_calderon_beds(self):

		calderon_info = pd.read_csv(self.url_calderon_info, sep='\t', names=['cellType'])
		calderon_info = calderon_info.set_index('cellType', drop=True)

		calderon_beds = []
		for ct in calderon_info.index.tolist():
			calderon_beds.append("%s/%s_L-tron.bed"%(self.url_calderon_beds, ct))

		return calderon_beds


	def __loading_calderon_bigwig(self):

		calderon_info = pd.read_csv(self.url_calderon_info, sep='\t', names=['cellType'])
		calderon_info = calderon_info.set_index('cellType', drop=True)

		calderon_bigwigs = []
		for ct in calderon_info.index.tolist():
			calderon_bigwigs.append("%s/%s.bw"%(self.url_calderon_bigwigs, ct))

		return calderon_bigwigs


	def __loading_calderon_data(self):

		calderon_beds = self.__loading_calderon_beds()
		calderon_bigwigs = self.__loading_calderon_bigwig()

		return calderon_beds, calderon_bigwigs


	def __loading_ludwig2019_beds(self):

		ludwig2019_info = pd.read_csv(self.url_ludwig2019_info, sep='\t', names=['cellType'])
		ludwig2019_info = ludwig2019_info.set_index('cellType', drop=True)

		ludwig2019_beds = []
		for ct in ludwig2019_info.index.tolist():
			ludwig2019_beds.append("%s/%s_L-tron.bed"%(self.url_ludwig2019_beds, ct))

		return ludwig2019_beds


	def __loading_ludwig2019_bigwig(self):

		ludwig2019_info = pd.read_csv(self.url_ludwig2019_info, sep='\t', names=['cellType'])
		ludwig2019_info = ludwig2019_info.set_index('cellType', drop=True)

		ludwig2019_bigwigs = []
		for ct in ludwig2019_info.index.tolist():
			ludwig2019_bigwigs.append("%s/%s.bw"%(self.url_ludwig2019_bigwigs, ct))

		return ludwig2019_bigwigs


	def __loading_ludwig2019_data(self):

		ludwig2019_beds = self.__loading_ludwig2019_beds()
		ludwig2019_bigwigs = self.__loading_ludwig2019_bigwig()

		return ludwig2019_beds, ludwig2019_bigwigs


	def __loading_MPAL_lowGr_info(self):

		MPAL_lowGr_info = pd.read_csv(self.url_MPAL_lowGr_info, sep='\t')
		MPAL_lowGr_info = MPAL_lowGr_info.set_index('cellType', drop=True)
		MPAL_lowGr_info.index.name = None
		MPAL_lowGr_info['scCounts'] = MPAL_lowGr_info['scCounts'].astype(float)

		return MPAL_lowGr_info


	def __loading_MPAL_lowGr_beds(self, MPAL_lowGr_info):
	
		MPAL_lowGr_beds = []
		
		for ct in MPAL_lowGr_info.index.tolist():
			MPAL_lowGr_beds.append("%s/%s_L-tron.bed"%(self.url_MPAL_lowGr_beds, ct))
		
		return MPAL_lowGr_beds


	def __loading_MPAL_lowGr_bigwig(self, MPAL_lowGr_info):

		MPAL_lowGr_bigwigs = []
		for ct in MPAL_lowGr_info.index.tolist():
			MPAL_lowGr_bigwigs.append("%s/%s.bw"%(self.url_MPAL_lowGr_bigwigs, ct))

		return MPAL_lowGr_bigwigs


	def __loading_MPAL_lowGr_data(self):

		MPAL_lowGr_info = self.__loading_MPAL_lowGr_info()
		MPAL_lowGr_beds = self.__loading_MPAL_lowGr_beds(MPAL_lowGr_info)
		MPAL_lowGr_bigwigs = self.__loading_MPAL_lowGr_bigwig(MPAL_lowGr_info)

		return MPAL_lowGr_info, MPAL_lowGr_beds, MPAL_lowGr_bigwigs


	def __loading_super_PBMC_info(self):

		super_PBMC_info = pd.read_csv(self.url_super_PBMC_info, sep='\t')
		super_PBMC_info = super_PBMC_info.set_index('cellType', drop=True)
		super_PBMC_info.index.name = None
		super_PBMC_info['scCounts'] = super_PBMC_info['scCounts'].astype(float)

		return super_PBMC_info


	def __loading_super_PBMC_beds(self, super_PBMC_info):
	
		super_PBMC_beds = []
		
		for ct in super_PBMC_info.index.tolist():
			super_PBMC_beds.append("%s/%s_L-tron.bed"%(self.url_super_PBMC_beds, ct))
		
		return super_PBMC_beds


	def __loading_super_PBMC_bigwig(self, super_PBMC_info):

		super_PBMC_bigwigs = []
		for ct in super_PBMC_info.index.tolist():
			super_PBMC_bigwigs.append("%s/%s.bw"%(self.url_super_PBMC_bigwigs, ct))

		return super_PBMC_bigwigs


	def __loading_super_PBMC_data(self):

		super_PBMC_info = self.__loading_super_PBMC_info()
		super_PBMC_beds = self.__loading_super_PBMC_beds(super_PBMC_info)
		super_PBMC_bigwigs = self.__loading_super_PBMC_bigwig(super_PBMC_info)

		return super_PBMC_info, super_PBMC_beds, super_PBMC_bigwigs


	def __loading_Days7_10_13_17_beds(self):

		Days7_10_13_17_info = pd.read_csv(self.url_Days7_10_13_17_info, sep='\t', names=['cellType'])
		Days7_10_13_17_info = Days7_10_13_17_info.set_index('cellType', drop=True)

		Days7_10_13_17_beds = []
		for ct in Days7_10_13_17_info.index.tolist():
			Days7_10_13_17_beds.append("%s/%s_L-tron.bed"%(self.url_Days7_10_13_17_beds, ct))

		return Days7_10_13_17_beds


	def __loading_Days7_10_13_17_bigwig(self):

		Days7_10_13_17_info = pd.read_csv(self.url_Days7_10_13_17_info, sep='\t', names=['cellType'])
		Days7_10_13_17_info = Days7_10_13_17_info.set_index('cellType', drop=True)

		Days7_10_13_17_bigwigs = []
		for ct in Days7_10_13_17_info.index.tolist():
			Days7_10_13_17_bigwigs.append("%s/%s.bw"%(self.url_Days7_10_13_17_bigwigs, ct))

		return Days7_10_13_17_bigwigs


	def __loading_Days7_10_13_17_data(self):

		Days7_10_13_17_beds = self.__loading_Days7_10_13_17_beds()
		Days7_10_13_17_bigwigs = self.__loading_Days7_10_13_17_bigwig()

		return Days7_10_13_17_beds, Days7_10_13_17_bigwigs


	def __loading_immune_cell_beds(self):

		immune_cell_info = pd.read_csv(self.url_immune_cell_info, sep='\t', names=['cellType'])
		immune_cell_info = immune_cell_info.set_index('cellType', drop=True)

		immune_cell_beds = []
		for ct in immune_cell_info.index.tolist():
			immune_cell_beds.append("%s/%s_L-tron.bed"%(self.url_immune_cell_beds, ct))

		return immune_cell_beds


	def __loading_immune_cell_bigwig(self):

		immune_cell_info = pd.read_csv(self.url_immune_cell_info, sep='\t', names=['cellType'])
		immune_cell_info = immune_cell_info.set_index('cellType', drop=True)

		immune_cell_bigwigs = []
		for ct in immune_cell_info.index.tolist():
			immune_cell_bigwigs.append("%s/%s.bw"%(self.url_immune_cell_bigwigs, ct))

		return immune_cell_bigwigs


	def __loading_immune_cell_data(self):

		immune_cell_beds = self.__loading_immune_cell_beds()
		immune_cell_bigwigs = self.__loading_immune_cell_bigwig()

		return immune_cell_beds, immune_cell_bigwigs


	def __loading_pancreatic_pbmc_beds(self):

		pancreatic_pbmc_info = pd.read_csv(self.url_pancreatic_pbmc_info, sep='\t', names=['cellType'])
		pancreatic_pbmc_info = pancreatic_pbmc_info.set_index('cellType', drop=True)

		pancreatic_pbmc_beds = []
		for ct in pancreatic_pbmc_info.index.tolist():
			pancreatic_pbmc_beds.append("%s/%s_L-tron.bed"%(self.url_pancreatic_pbmc_beds, ct))

		return pancreatic_pbmc_beds


	def __loading_pancreatic_pbmc_bigwig(self):

		pancreatic_pbmc_info = pd.read_csv(self.url_pancreatic_pbmc_info, sep='\t', names=['cellType'])
		pancreatic_pbmc_info = pancreatic_pbmc_info.set_index('cellType', drop=True)

		pancreatic_pbmc_bigwigs = []
		for ct in pancreatic_pbmc_info.index.tolist():
			pancreatic_pbmc_bigwigs.append("%s/%s.bw"%(self.url_pancreatic_pbmc_bigwigs, ct))

		return pancreatic_pbmc_bigwigs


	def __loading_pancreatic_pbmc_data(self):

		pancreatic_pbmc_beds = self.__loading_pancreatic_pbmc_beds()
		pancreatic_pbmc_bigwigs = self.__loading_pancreatic_pbmc_bigwig()

		return pancreatic_pbmc_beds, pancreatic_pbmc_bigwigs


	def __loading_H1_hESCs_beds(self):

		H1_hESCs_info = pd.read_csv(self.url_H1_hESCs_info, sep='\t', names=['cellType'])
		H1_hESCs_info = H1_hESCs_info.set_index('cellType', drop=True)

		H1_hESCs_beds = []
		for ct in H1_hESCs_info.index.tolist():
			H1_hESCs_beds.append("%s/%s_L-tron.bed"%(self.url_H1_hESCs_beds, ct))

		return H1_hESCs_beds


	def __loading_H1_hESCs_bigwig(self):

		H1_hESCs_info = pd.read_csv(self.url_H1_hESCs_info, sep='\t', names=['cellType'])
		H1_hESCs_info = H1_hESCs_info.set_index('cellType', drop=True)

		H1_hESCs_bigwigs = []
		for ct in H1_hESCs_info.index.tolist():
			H1_hESCs_bigwigs.append("%s/%s.bw"%(self.url_H1_hESCs_bigwigs, ct))

		return H1_hESCs_bigwigs


	def __loading_H1_hESCs_data(self):

		H1_hESCs_beds = self.__loading_H1_hESCs_beds()
		H1_hESCs_bigwigs = self.__loading_H1_hESCs_bigwig()

		return H1_hESCs_beds, H1_hESCs_bigwigs


	def __loading_CAD_info(self):

		CAD_info = pd.read_csv(self.url_CAD_info, sep='\t')
		CAD_info = CAD_info.set_index('cellType', drop=True)
		CAD_info.index.name = None
		CAD_info['scCounts'] = CAD_info['scCounts'].astype(float)

		return CAD_info


	def __loading_CAD_beds(self, CAD_info):
	
		CAD_beds = []
		
		for ct in CAD_info.index.tolist():
			CAD_beds.append("%s/%s_L-tron.bed"%(self.url_CAD_beds, ct))
		
		return CAD_beds


	def __loading_CAD_bigwig(self, CAD_info):

		CAD_bigwigs = []
		for ct in CAD_info.index.tolist():
			CAD_bigwigs.append("%s/%s.bw"%(self.url_CAD_bigwigs, ct))

		return CAD_bigwigs


	def __loading_CAD_data(self):

		CAD_info = self.__loading_CAD_info()
		CAD_beds = self.__loading_CAD_beds(CAD_info)
		CAD_bigwigs = self.__loading_CAD_bigwig(CAD_info)

		return CAD_info, CAD_beds, CAD_bigwigs


	def __loading_custom_beds(self):

		custom_info = pd.read_csv(self.url_custom_info, sep='\t', names=['cellType'])
		custom_info = custom_info.set_index('cellType', drop=True)

		custom_beds = []
		for ct in custom_info.index.tolist():
			custom_beds.append("%s/%s_L-tron.bed"%(self.url_custom_beds, ct))

		return custom_beds


	def __loading_custom_bigwig(self):

		custom_info = pd.read_csv(self.url_custom_info, sep='\t', names=['cellType'])
		custom_info = custom_info.set_index('cellType', drop=True)

		custom_bigwigs = []
		for ct in custom_info.index.tolist():
			custom_bigwigs.append("%s/%s.bw"%(self.url_custom_bigwigs, ct))

		return custom_bigwigs


	def __loading_custom_data(self):

		custom_beds = self.__loading_custom_beds()
		custom_bigwigs = self.__loading_custom_bigwig()

		return custom_beds, custom_bigwigs


	def __loading_multiome_ery_info(self):

		multiome_ery_info = pd.read_csv(self.url_multiome_ery_info, sep='\t')
		multiome_ery_info = multiome_ery_info.set_index('cellType', drop=True)
		multiome_ery_info.index.name = None
		multiome_ery_info['scCounts'] = multiome_ery_info['scCounts'].astype(float)

		return multiome_ery_info


	def __loading_multiome_ery_beds(self, multiome_ery_info):
	
		multiome_ery_beds = []
		
		for ct in multiome_ery_info.index.tolist():
			multiome_ery_beds.append("%s/%s_L-tron.bed"%(self.url_multiome_ery_beds, ct))
		
		return multiome_ery_beds


	def __loading_multiome_ery_bigwig(self, multiome_ery_info):

		multiome_ery_bigwigs = []
		for ct in multiome_ery_info.index.tolist():
			multiome_ery_bigwigs.append("%s/%s.bw"%(self.url_multiome_ery_bigwigs, ct))

		return multiome_ery_bigwigs


	def __loading_multiome_ery_data(self):

		multiome_ery_info = self.__loading_multiome_ery_info()
		multiome_ery_beds = self.__loading_multiome_ery_beds(multiome_ery_info)
		multiome_ery_bigwigs = self.__loading_multiome_ery_bigwig(multiome_ery_info)

		return multiome_ery_info, multiome_ery_beds, multiome_ery_bigwigs


	def __prepare_data(self, beds):

		df_data = pd.DataFrame()

		peak_area = []
		ps = []
		celltypes = []
		print("Calculating (1) total number of base-pairs within peaks and (2) total number of base-pairs within peaks divided by uniquely mappable base-pairs ...")
		for bed in tqdm(beds):
			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
			
			celltype = bed.split("/")[-1].replace("_L-tron.bed","")
			celltypes.append(celltype)
			
			# total number of base-pairs within peaks
			tot_bp_within_peaks = np.sum(tmp['end']-tmp['start'])
			peak_area.append(tot_bp_within_peaks)
			# print(tot_bp_within_peaks)
			
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
		df_bed = pybedtools.BedTool.from_dataframe(self.genetic)
		print("Intersecting genetic with peak regions ...")
		for bed in tqdm(beds):
			# try:
			tmp = pd.read_csv(bed, sep="\t", names=["chrom", "start", "end"])
			tmp_bed = pybedtools.BedTool.from_dataframe(tmp)
			intersect_bed = df_bed.intersect(tmp_bed)
			intersect_bed = intersect_bed.to_dataframe()
			xs.append(intersect_bed.shape[0])
			# except:
				# xs.append(0.0)	
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
		df_data['MEAN'] = col.mean(axis=1)

		return df_data


	def __clean_tmp(self):

		print("Cleaning temporary files ...")
		shutil.rmtree(self.tmp)


	def process_genetics_for_catlas(self, origin="adult"):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		if origin not in ["adult", "fetal"]:
			sys.exit("Error, origin must be 'adult' or 'fetal'")

		number_of_genetic = self.__load_genetic(subsampling=True)

		self.origin = origin

		catlas_info, catlas_beds, _ = self.__loading_catlas_data(self.origin)
		df_catlas = self.__prepare_data(catlas_beds)
		self.__create_background(number_of_genetic)
		df_catlas = self.__add_shuffled_background(df_catlas, catlas_beds)
		df_catlas = self.__add_statistics(df_catlas, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_catlas.to_csv(self.output+os.sep+"catlas_%s_statistics.csv"%self.origin, sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_catlas = df_catlas

		return df_catlas


	def process_genetics_for_calderon(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		calderon_beds, _ = self.__loading_calderon_data()
		df_calderon = self.__prepare_data(calderon_beds)
		self.__create_background(number_of_genetic)
		df_calderon = self.__add_shuffled_background(df_calderon, calderon_beds)
		df_calderon = self.__add_statistics(df_calderon, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_calderon.to_csv(self.output+os.sep+"calderon_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_calderon = df_calderon

		return df_calderon


	def process_genetics_for_ludwig2019(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		ludwig2019_beds, _ = self.__loading_ludwig2019_data()
		df_ludwig2019 = self.__prepare_data(ludwig2019_beds)
		self.__create_background(number_of_genetic)
		df_ludwig2019 = self.__add_shuffled_background(df_ludwig2019, ludwig2019_beds)
		df_ludwig2019 = self.__add_statistics(df_ludwig2019, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_ludwig2019.to_csv(self.output+os.sep+"ludwig2019_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_ludwig2019 = df_ludwig2019

		return df_ludwig2019


	def process_genetics_for_MPAL_lowGr(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)

		
		MPAL_lowGr_info, MPAL_lowGr_beds, _ = self.__loading_MPAL_lowGr_data()
		df_MPAL_lowGr = self.__prepare_data(MPAL_lowGr_beds)
		self.__create_background(number_of_genetic)
		df_MPAL_lowGr = self.__add_shuffled_background(df_MPAL_lowGr, MPAL_lowGr_beds)
		df_MPAL_lowGr = self.__add_statistics(df_MPAL_lowGr, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_MPAL_lowGr.to_csv(self.output+os.sep+"MPAL_lowGr_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_MPAL_lowGr = df_MPAL_lowGr

		return df_MPAL_lowGr


	def process_genetics_for_super_PBMC(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)

		super_PBMC_info, super_PBMC_beds, _ = self.__loading_super_PBMC_data()
		df_super_PBMC = self.__prepare_data(super_PBMC_beds)
		self.__create_background(number_of_genetic)
		df_super_PBMC = self.__add_shuffled_background(df_super_PBMC, super_PBMC_beds)
		df_super_PBMC = self.__add_statistics(df_super_PBMC, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_super_PBMC.to_csv(self.output+os.sep+"super_PBMC_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_super_PBMC = df_super_PBMC

		return df_super_PBMC

	def process_genetics_for_Days7_10_13_17(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		Days7_10_13_17_beds, _ = self.__loading_Days7_10_13_17_data()
		df_Days7_10_13_17 = self.__prepare_data(Days7_10_13_17_beds)
		self.__create_background(number_of_genetic)
		df_Days7_10_13_17 = self.__add_shuffled_background(df_Days7_10_13_17, Days7_10_13_17_beds)
		df_Days7_10_13_17 = self.__add_statistics(df_Days7_10_13_17, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_Days7_10_13_17.to_csv(self.output+os.sep+"Days7_10_13_17_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_Days7_10_13_17 = df_Days7_10_13_17

		return df_Days7_10_13_17


	def process_genetics_for_immune_cell(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		immune_cell_beds, _ = self.__loading_immune_cell_data()
		df_immune_cell = self.__prepare_data(immune_cell_beds)
		self.__create_background(number_of_genetic)
		df_immune_cell = self.__add_shuffled_background(df_immune_cell, immune_cell_beds)
		df_immune_cell = self.__add_statistics(df_immune_cell, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_immune_cell.to_csv(self.output+os.sep+"immune_cell_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_immune_cell = df_immune_cell

		return df_immune_cell


	def process_genetics_for_pancreatic_pbmc(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		pancreatic_pbmc_beds, _ = self.__loading_pancreatic_pbmc_data()
		df_pancreatic_pbmc = self.__prepare_data(pancreatic_pbmc_beds)
		self.__create_background(number_of_genetic)
		df_pancreatic_pbmc = self.__add_shuffled_background(df_pancreatic_pbmc, pancreatic_pbmc_beds)
		df_pancreatic_pbmc = self.__add_statistics(df_pancreatic_pbmc, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_pancreatic_pbmc.to_csv(self.output+os.sep+"pancreatic_pbmc_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_pancreatic_pbmc = df_pancreatic_pbmc

		return df_pancreatic_pbmc


	def process_genetics_for_H1_hESCs(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		H1_hESCs_beds, _ = self.__loading_H1_hESCs_data()
		df_H1_hESCs = self.__prepare_data(H1_hESCs_beds)
		self.__create_background(number_of_genetic)
		df_H1_hESCs = self.__add_shuffled_background(df_H1_hESCs, H1_hESCs_beds)
		df_H1_hESCs = self.__add_statistics(df_H1_hESCs, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_H1_hESCs.to_csv(self.output+os.sep+"H1_hESCs_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_H1_hESCs = df_H1_hESCs

		return df_H1_hESCs


	def process_genetics_for_CAD(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)

		CAD_info, CAD_beds, _ = self.__loading_CAD_data()
		df_CAD = self.__prepare_data(CAD_beds)
		self.__create_background(number_of_genetic)
		df_CAD = self.__add_shuffled_background(df_CAD, CAD_beds)
		df_CAD = self.__add_statistics(df_CAD, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_CAD.to_csv(self.output+os.sep+"CAD_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_CAD = df_CAD

		return df_CAD


	def process_genetics_for_customdata(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		custom_beds, _ = self.__loading_custom_data()
		df_custom = self.__prepare_data(custom_beds)
		self.__create_background(number_of_genetic)
		df_custom = self.__add_shuffled_background(df_custom, custom_beds)
		df_custom = self.__add_statistics(df_custom, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_custom.to_csv(self.output+os.sep+"custom_data_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_custom = df_custom

		return df_custom


	def process_genetics_for_multiome_ery(self):

		if self.genetic is None:
			sys.exit("Error, missing genetic file.")

		number_of_genetic = self.__load_genetic(subsampling=True)
		
		multiome_ery_info, multiome_ery_beds, _ = self.__loading_multiome_ery_data()
		df_multiome_ery = self.__prepare_data(multiome_ery_beds)
		self.__create_background(number_of_genetic)
		df_multiome_ery = self.__add_shuffled_background(df_multiome_ery, multiome_ery_beds)
		df_multiome_ery = self.__add_statistics(df_multiome_ery, number_of_genetic)

		if not os.path.exists(self.output):
			os.makedirs(self.output)
		df_multiome_ery.to_csv(self.output+os.sep+"multiome_ery_statistics.csv", sep="\t")

		self.__clean_tmp()
		print("Processing genetics finished.")

		self.df_multiome_ery = df_multiome_ery

		return df_multiome_ery


	def plot_genetics_for_catlas(self, umap=False, show: bool = True):

		if self.df_catlas.empty:
			sys.exit("Error, process_genetics_for_catlas has to be run.")
				
		fetal_list = set(self.df_catlas[self.df_catlas.index.str.startswith("Fetal")].index.tolist())
		adult_list = set(self.df_catlas.index.tolist()).difference(fetal_list)

		catlas_info = pd.read_csv(self.url_catlas_info, sep='\t')
		catlas_info = catlas_info.set_index('cellType', drop=True)
		catlas_info.index.name = None
		catlas_info['scCounts'] = catlas_info['scCounts'].astype(float)

		print("Generating plots for %s CATLAS ..."%self.origin)

		if self.origin == 'fetal':
			df_fetal             = self.df_catlas[self.df_catlas.index.isin(fetal_list)]
			df_fetal             = df_fetal.sort_index()
			df_fetal['scCounts'] = catlas_info[catlas_info.index.isin(fetal_list)].sort_index()['scCounts']
			df_fetal = df_fetal.replace([np.nan, -np.inf, np.inf], 0.0)
			# df_fetal.index       = df_fetal.index.str.replace("Fetal_", "")
		else:
			df_adult             = self.df_catlas[self.df_catlas.index.isin(adult_list)]
			df_adult             = df_adult.sort_index()
			df_adult['scCounts'] = catlas_info[catlas_info.index.isin(adult_list)].sort_index()['scCounts']
			df_adult = df_adult.replace([np.nan, -np.inf, np.inf], 0.0)

		# plot fetal catlas
		if self.origin == 'fetal':
			fig_fetal, axs = plt.subplots(2, 1, sharex=True, figsize=(50, 3))

			cols = ["MEAN", "scCounts"]
			cbar_ax = fig_fetal.add_axes([.92, .53095, .02, .3475])

			for i, ax in enumerate(axs):
				if i in [0]:
					sns.heatmap(df_fetal[cols[i]].values.reshape(1, df_fetal.shape[0]), annot=True, fmt=".2f",
								ax=ax,
								cmap="coolwarm",
								annot_kws={'rotation': 90},
								cbar=i == 0,
								cbar_ax=None if i else cbar_ax,
								cbar_kws={'label': 'Enrichment score'})
				else:
					sns.heatmap(df_fetal[cols[i]].values.reshape(1, df_fetal.shape[0]), annot=True, fmt=".0f", mask=list(df_fetal[cols[i]] < 500),
								ax=ax,
								annot_kws={'rotation': 90},
								cmap=sns.color_palette(["#50C878"]),
								cbar=False)
					sns.heatmap(df_fetal[cols[i]].values.reshape(1, df_fetal.shape[0]), annot=True, fmt=".0f", mask=list(df_fetal[cols[i]] >= 500),
								ax=ax,
								annot_kws={'rotation': 90},
								cmap=sns.color_palette(["#EE4B2B"]),
								cbar=False, vmin=0, vmax=0.5)

			axs[0].set_yticklabels([])
			axs[0].set_ylabel("Catlas Fetal", rotation=0)
			axs[0].yaxis.set_label_coords(-.025, .425)

			axs[1].set_yticklabels([])
			axs[1].set_ylabel("single-cells", rotation=0)
			axs[1].yaxis.set_label_coords(-.025, .425)

			axs[1].set_xticklabels(df_fetal.index.tolist(), rotation=90)

			plt.suptitle("SNPs enrichment", y=1.1)

			if show==True:
				plt.show()
			else:
				plt.ioff()
				plt.close()
				fig_fetal.savefig("%s/catlas_fetal_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

			if show==False:
				return fig_fetal

		else:
			# plot adult catlas
			fig_adult, axs = plt.subplots(2, 1, sharex=True, figsize=(50, 3))

			cols = ["MEAN", "scCounts"]
			cbar_ax = fig_adult.add_axes([.92, .53095, .02, .3475])

			for i, ax in enumerate(axs):
				if i in [0]:
					sns.heatmap(df_adult[cols[i]].values.reshape(1, df_adult.shape[0]), annot=True, fmt=".2f",
								ax=ax,
								cmap="coolwarm",
								annot_kws={'rotation': 90},
								cbar=i == 0,
								cbar_ax=None if i else cbar_ax,
								cbar_kws={'label': 'Enrichment score'})
				else:
					sns.heatmap(df_adult[cols[i]].values.reshape(1, df_adult.shape[0]), annot=True, fmt=".0f", mask=list(df_adult[cols[i]] < 500),
								ax=ax,
								annot_kws={'rotation': 90},
								cmap=sns.color_palette(["#50C878"]),
								cbar=False)
					sns.heatmap(df_adult[cols[i]].values.reshape(1, df_adult.shape[0]), annot=True, fmt=".0f", mask=list(df_adult[cols[i]] >= 500),
								ax=ax,
								annot_kws={'rotation': 90},
								cmap=sns.color_palette(["#EE4B2B"]),
								cbar=False, vmin=0, vmax=0.5)

			axs[0].set_yticklabels([])
			axs[0].set_ylabel("Catlas Adult", rotation=0)
			axs[0].yaxis.set_label_coords(-.025, .425)

			axs[1].set_yticklabels([])
			axs[1].set_ylabel("single-cells", rotation=0)
			axs[1].yaxis.set_label_coords(-.025, .425)

			axs[1].set_xticklabels(df_adult.index.tolist(), rotation=90)

			plt.suptitle("SNPs enrichment", y=1.1)
			
			if show==True:
				plt.show()
			else:
				plt.ioff()
				plt.close()
				fig_adult.savefig("%s/catlas_adult_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

			if show==False:
				return fig_adult


	def plot_genetics_for_calderon(self, df_calderon=None, show: bool = True):

		if self.df_calderon.empty:
			sys.exit("Error, process_genetics_for_calderon has to be run.")

		self.df_calderon = self.df_calderon.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(50, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_calderon["MEAN"].values.reshape(1, self.df_calderon.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("Calderon", rotation=0)
		axs.yaxis.set_label_coords(-.015, .425)

		axs.set_xticklabels(self.df_calderon.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/calderon_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_ludwig2019(self, df_ludwig2019=None, show: bool = True):

		if self.df_ludwig2019.empty:
			sys.exit("Error, process_genetics_for_ludwig2019 has to be run.")

		self.df_ludwig2019 = self.df_ludwig2019.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_ludwig2019["MEAN"].values.reshape(1, self.df_ludwig2019.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("Ludwig2019", rotation=0)
		axs.yaxis.set_label_coords(-.085, .425)

		axs.set_xticklabels(self.df_ludwig2019.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/ludwig2019_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_MPAL_lowGr(self, show: bool = True):

		if self.df_MPAL_lowGr.empty:
			sys.exit("Error, process_genetics_for_MPAL_lowGr has to be run.")

		MPAL_lowGr_info = pd.read_csv(self.url_MPAL_lowGr_info, sep='\t')
		MPAL_lowGr_info = MPAL_lowGr_info.set_index('cellType', drop=True)
		MPAL_lowGr_info.index.name = None
		MPAL_lowGr_info['scCounts'] = MPAL_lowGr_info['scCounts'].astype(float)

		df_MPAL_lowGr = self.df_MPAL_lowGr
		df_MPAL_lowGr['scCounts'] = MPAL_lowGr_info['scCounts']
		df_MPAL_lowGr = df_MPAL_lowGr.replace([np.nan, -np.inf, np.inf], 0.0)

		print("Generating plots for MPAL_lowGr ...")

		fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 3))

		cols = ["MEAN", "scCounts"]
		cbar_ax = fig.add_axes([.92, .53095, .02, .3475])

		for i, ax in enumerate(axs):
			if i in [0]:
				sns.heatmap(df_MPAL_lowGr[cols[i]].values.reshape(1, df_MPAL_lowGr.shape[0]), annot=True, fmt=".2f",
							ax=ax,
							cmap="coolwarm",
							annot_kws={'rotation': 90},
							cbar=i == 0,
							cbar_ax=None if i else cbar_ax,
							cbar_kws={'label': 'Enrichment score'})
			else:
				sns.heatmap(df_MPAL_lowGr[cols[i]].values.reshape(1, df_MPAL_lowGr.shape[0]), annot=True, fmt=".0f", mask=list(df_MPAL_lowGr[cols[i]] < 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#50C878"]),
							cbar=False)
				sns.heatmap(df_MPAL_lowGr[cols[i]].values.reshape(1, df_MPAL_lowGr.shape[0]), annot=True, fmt=".0f", mask=list(df_MPAL_lowGr[cols[i]] >= 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#EE4B2B"]),
							cbar=False, vmin=0, vmax=0.5)

		axs[0].set_yticklabels([])
		axs[0].set_ylabel("MPAL_lowGr", rotation=0)
		axs[0].yaxis.set_label_coords(-.095, .425)

		axs[1].set_yticklabels([])
		axs[1].set_ylabel("single-cells", rotation=0)
		axs[1].yaxis.set_label_coords(-.095, .425)

		axs[1].set_xticklabels(df_MPAL_lowGr.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/MPAL_lowGr_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

		if show==False:
			return fig


	def plot_genetics_for_super_PBMC(self, umap=False, show: bool = True):

		if self.df_super_PBMC.empty:
			sys.exit("Error, process_genetics_for_super_PBMC has to be run.")

		super_PBMC_info = pd.read_csv(self.url_super_PBMC_info, sep='\t')
		super_PBMC_info = super_PBMC_info.set_index('cellType', drop=True)
		super_PBMC_info.index.name = None
		super_PBMC_info['scCounts'] = super_PBMC_info['scCounts'].astype(float)

		df_super_PBMC = self.df_super_PBMC
		df_super_PBMC['scCounts'] = super_PBMC_info['scCounts']
		df_super_PBMC = df_super_PBMC.replace([np.nan, -np.inf, np.inf], 0.0)

		print("Generating plots for super_PBMC ...")

		fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 3))

		cols = ["MEAN", "scCounts"]
		cbar_ax = fig.add_axes([.92, .53095, .02, .3475])

		for i, ax in enumerate(axs):
			if i in [0]:
				sns.heatmap(df_super_PBMC[cols[i]].values.reshape(1, df_super_PBMC.shape[0]), annot=True, fmt=".2f",
							ax=ax,
							cmap="coolwarm",
							annot_kws={'rotation': 90},
							cbar=i == 0,
							cbar_ax=None if i else cbar_ax,
							cbar_kws={'label': 'Enrichment score'})
			else:
				sns.heatmap(df_super_PBMC[cols[i]].values.reshape(1, df_super_PBMC.shape[0]), annot=True, fmt=".0f", mask=list(df_super_PBMC[cols[i]] < 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#50C878"]),
							cbar=False)
				sns.heatmap(df_super_PBMC[cols[i]].values.reshape(1, df_super_PBMC.shape[0]), annot=True, fmt=".0f", mask=list(df_super_PBMC[cols[i]] >= 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#EE4B2B"]),
							cbar=False, vmin=0, vmax=0.5)

		axs[0].set_yticklabels([])
		axs[0].set_ylabel("super_PBMC", rotation=0)
		axs[0].yaxis.set_label_coords(-.095, .425)

		axs[1].set_yticklabels([])
		axs[1].set_ylabel("single-cells", rotation=0)
		axs[1].yaxis.set_label_coords(-.095, .425)

		axs[1].set_xticklabels(df_super_PBMC.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/super_PBMC_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

		if umap:
			enScore = self.df_super_PBMC
			umap_df = pd.read_csv(self.url_super_PBMC_umap, sep='\t')

			mapCluster = []
			for c in umap_df['Annotation']:
			    mapCluster.append(enScore.loc[c]['MEAN'])
			umap_df['enScore'] = mapCluster
			umap_df['enScore'] = umap_df['enScore'].astype(np.float64)
			umap_df['Annotation'] = umap_df['Annotation'].astype('category')

			fig_umap, axs = plt.subplots(1,2,figsize=(18,6))
			# sns.set(font_scale=1.5)
			sns.set_style("white")

			axs[0] = sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="Annotation", ax=axs[0])
			axs[0].spines[['top', 'right']].set_visible(False)
			axs[0].set(xlabel='UMAP1', ylabel='UMAP2')
			axs[0].legend(frameon=False, ncol=3)

			axs[1] = sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="enScore", palette='coolwarm', ax=axs[1])
			axs[1].spines[['left', 'top', 'right']].set_visible(False)
			axs[1].tick_params(left = False, right = False , labelleft = False) 
			norm = plt.Normalize(umap_df['enScore'].min(), umap_df['enScore'].max())
			sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
			sm.set_array([])
			axs[1].get_legend().remove()

			cbaxes = fig_umap.add_axes([1.01, 0.1, 0.01, 0.85])  

			clb = axs[1].figure.colorbar(sm, cax = cbaxes)
			clb.ax.get_yaxis().labelpad = 30
			clb.ax.set_ylabel('Enrichment score', rotation=270)
			axs[1].set(xlabel='UMAP1', ylabel='')

			# sns.despine(offset=10, trim=False)
			plt.title("")
			plt.tight_layout()

			if show==True:
				plt.show()
			else:
				plt.ioff()
				plt.close()
				fig_umap.savefig("%s/super_PBMC_UMAP_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

		if show==False:
			if umap:
				return fig, fig_umap
			else:
				return fig


	def plot_genetics_for_Days7_10_13_17(self, df_Days7_10_13_17=None, show: bool = True):

		if self.df_Days7_10_13_17.empty:
			sys.exit("Error, process_genetics_for_Days7_10_13_17 has to be run.")

		self.df_Days7_10_13_17 = self.df_Days7_10_13_17.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_Days7_10_13_17["MEAN"].values.reshape(1, self.df_Days7_10_13_17.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("Days 7, 10, 13, and 17", rotation=0)
		axs.yaxis.set_label_coords(-.125, .425)

		axs.set_xticklabels(self.df_Days7_10_13_17.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/Days7_10_13_17_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_immune_cell(self, df_immune_cell=None, show: bool = True):

		if self.df_immune_cell.empty:
			sys.exit("Error, process_genetics_for_immune_cell has to be run.")

		self.df_immune_cell = self.df_immune_cell.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_immune_cell["MEAN"].values.reshape(1, self.df_immune_cell.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("Immune cells", rotation=0)
		axs.yaxis.set_label_coords(-.08, .425)

		axs.set_xticklabels(self.df_immune_cell.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/immune_cell_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_pancreatic_pbmc(self, df_pancreatic_pbmc=None, show: bool = True):

		if self.df_pancreatic_pbmc.empty:
			sys.exit("Error, process_genetics_for_pancreatic_pbmc has to be run.")

		self.df_pancreatic_pbmc = self.df_pancreatic_pbmc.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_pancreatic_pbmc["MEAN"].values.reshape(1, self.df_pancreatic_pbmc.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("Pancreatic PBMC", rotation=0)
		axs.yaxis.set_label_coords(-.1, .425)

		axs.set_xticklabels(self.df_pancreatic_pbmc.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/pancreatic_pbmc_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_H1_hESCs(self, df_H1_hESCs=None, show: bool = True):

		if self.df_H1_hESCs.empty:
			sys.exit("Error, process_genetics_for_H1_hESCs has to be run.")

		self.df_H1_hESCs = self.df_H1_hESCs.replace([np.nan, -np.inf, np.inf], 0.0)
		
		fig, axs = plt.subplots(1, 1, sharex=True, figsize=(6, 1.5))

		cbar_ax = fig.add_axes([.92, .13095, .02, .735])

		sns.heatmap(self.df_H1_hESCs["MEAN"].values.reshape(1, self.df_H1_hESCs.shape[0]), annot=True, fmt=".2f",
					ax=axs,
					cmap="coolwarm",
					annot_kws={'rotation': 90},
					cbar_ax= cbar_ax,
					cbar_kws={'label': 'Enrichment score'})

		axs.set_yticklabels([])
		axs.set_ylabel("H1 hESCs", rotation=0)
		axs.yaxis.set_label_coords(-.12, .425)

		axs.set_xticklabels(self.df_H1_hESCs.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/H1_hESCs_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)
			return fig


	def plot_genetics_for_CAD(self, umap=False, show: bool = True):

		if self.df_CAD.empty:
			sys.exit("Error, process_genetics_for_CAD has to be run.")

		CAD_info = pd.read_csv(self.url_CAD_info, sep='\t')
		CAD_info = CAD_info.set_index('cellType', drop=True)
		CAD_info.index.name = None
		CAD_info['scCounts'] = CAD_info['scCounts'].astype(float)

		df_CAD = self.df_CAD
		df_CAD['scCounts'] = CAD_info['scCounts']
		df_CAD = df_CAD.replace([np.nan, -np.inf, np.inf], 0.0)

		print("Generating plots for CAD ...")

		fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 3))

		cols = ["MEAN", "scCounts"]
		cbar_ax = fig.add_axes([.92, .53095, .02, .3475])

		for i, ax in enumerate(axs):
			if i in [0]:
				sns.heatmap(df_CAD[cols[i]].values.reshape(1, df_CAD.shape[0]), annot=True, fmt=".2f",
							ax=ax,
							cmap="coolwarm",
							annot_kws={'rotation': 90},
							cbar=i == 0,
							cbar_ax=None if i else cbar_ax,
							cbar_kws={'label': 'Enrichment score'})
			else:
				sns.heatmap(df_CAD[cols[i]].values.reshape(1, df_CAD.shape[0]), annot=True, fmt=".0f", mask=list(df_CAD[cols[i]] < 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#50C878"]),
							cbar=False)
				sns.heatmap(df_CAD[cols[i]].values.reshape(1, df_CAD.shape[0]), annot=True, fmt=".0f", mask=list(df_CAD[cols[i]] >= 500),
							ax=ax,
							annot_kws={'rotation': 90},
							cmap=sns.color_palette(["#EE4B2B"]),
							cbar=False, vmin=0, vmax=0.5)

		axs[0].set_yticklabels([])
		axs[0].set_ylabel("CAD", rotation=0)
		axs[0].yaxis.set_label_coords(-.095, .425)

		axs[1].set_yticklabels([])
		axs[1].set_ylabel("single-cells", rotation=0)
		axs[1].yaxis.set_label_coords(-.095, .425)

		axs[1].set_xticklabels(df_CAD.index.tolist(), rotation=90)

		plt.suptitle("SNPs enrichment", y=1.1)

		if show==True:
			plt.show()
		else:
			plt.ioff()
			plt.close()
			fig.savefig("%s/CAD_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

		if umap:
			enScore = self.df_CAD
			umap_df = pd.read_csv(self.url_CAD_umap, sep='\t')

			mapCluster = []
			for c in umap_df['Annotation']:
			    mapCluster.append(enScore.loc[c]['MEAN'])
			umap_df['enScore'] = mapCluster
			umap_df['enScore'] = umap_df['enScore'].astype(np.float64)
			umap_df['Annotation'] = umap_df['Annotation'].astype('category')

			fig_umap, axs = plt.subplots(1,2,figsize=(18,6))
			# sns.set(font_scale=1.5)
			sns.set_style("white")

			axs[0] = sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="Annotation", ax=axs[0])
			axs[0].spines[['top', 'right']].set_visible(False)
			axs[0].set(xlabel='UMAP1', ylabel='UMAP2')
			axs[0].legend(frameon=False, ncol=3)

			axs[1] = sns.scatterplot(data=umap_df, x="UMAP1", y="UMAP2", hue="enScore", palette='coolwarm', ax=axs[1])
			axs[1].spines[['left', 'top', 'right']].set_visible(False)
			axs[1].tick_params(left = False, right = False , labelleft = False) 
			norm = plt.Normalize(umap_df['enScore'].min(), umap_df['enScore'].max())
			sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
			sm.set_array([])
			axs[1].get_legend().remove()

			cbaxes = fig_umap.add_axes([1.01, 0.1, 0.01, 0.85])  

			clb = axs[1].figure.colorbar(sm, cax = cbaxes)
			clb.ax.get_yaxis().labelpad = 30
			clb.ax.set_ylabel('Enrichment score', rotation=270)
			axs[1].set(xlabel='UMAP1', ylabel='')

			# sns.despine(offset=10, trim=False)
			plt.title("")
			plt.tight_layout()

			if show==True:
				plt.show()
			else:
				plt.ioff()
				plt.close()
				fig_umap.savefig("%s/CAD_UMAP_enrichment_per_cell_type.pdf"%self.output, bbox_inches = 'tight', dpi=200)

		if show==False:
			if umap:
				return fig, fig_umap
			else:
				return fig


	# ### Variant prioritisation for V2

	# ### for parallelisation ###
	# ###########################
	# def _compute_coverage(self, bwi, variants):

	# 	bw = pyBigWig.open(bwi)
	# 	celltype=bwi.split("/")[-1].replace(".bw", "")
	# 	variant_values = [["%s_1"%celltype, "%s_2"%celltype, "%s_3"%celltype, "%s_4"%celltype, "%s_5"%celltype, "%s_6"%celltype, "%s_7"%celltype, "%s_8"%celltype, "%s_9"%celltype]]
		
	# 	for i in range(variants.shape[0]):
	# 		chrom=variants.iloc[i]['CHR_ID']
	# 		start=variants.iloc[i]['CHR_POS']-4
	# 		end=variants.iloc[i]['CHR_POS+1']+4
	# 		variant_value = bw.values(chrom, start, end, numpy=False)
	# 		variant_values.append(variant_value)
	# 	bw.close()

	# 	return [variant_values]
	# ###########################
	# ###########################


	# def multicoverage_genetics_for_catlas(self, variants=None, origin="adult", filename="multicoverage_catlas"):
		
	# 	if (variants is None) and (self.genetic is None):
	# 		sys.exit("Error, missing genetic file.")
	# 	elif variants is None:
	# 		variants = self.genetic
	# 	elif isinstance(variants, str):
	# 		self.__load_genetic()
	# 		variants = self.genetic

	# 	if (origin not in ["adult", "fetal"]) and (self.origin not in ["adult", "fetal"]):
	# 		sys.exit("Error, origin must be 'adult' or 'fetal'")
	# 	elif origin is None:
	# 		origin = self.origin
	# 	elif self.origin == "":
	# 		self.origin = origin
					
	# 	print("Computing catlas %s multicoverage might take time ..."%origin)

	# 	_, _, catlas_bigwigs = self.__loading_catlas_data(origin)

	# 	if not os.path.exists(self.tmp):
	# 		os.makedirs(self.tmp)

	# 	bed = self.tmp+"/multicoverage.bed"
	# 	variants.to_csv(bed, sep='\t', index=False, header=None)

	# 	if not os.path.exists(self.output):
	# 		os.makedirs(self.output)

	# 	commands = ["computeMatrix reference-point -S %s -R %s -a 2 -b 3 --missingDataAsZero --binSize 1 -o %s/multicoverage.csv --outFileNameMatrix %s/multicoverage.csv"%(" ".join(catlas_bigwigs), bed, self.output, self.output)]
	# 	subprocess.run(commands, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	# 	df = pd.read_csv("%s/multicoverage.csv"%(self.output), sep='\t', comment='#')
	# 	columns = df.columns[1:]
	# 	df.drop(columns[-1], inplace=True, axis=1)
	# 	df.columns = columns
	# 	index_id = variants['CHR_ID'].astype(str)+"-"+variants['CHR_POS'].astype(str)+"-"+variants['SNPS'].astype(str)
	# 	df.index = index_id
		
	# 	new_cols = []
	# 	for c in df.columns:
	# 		if ".1" in c:
	# 			new_cols.append(c.replace(".1", "_1"))
	# 		elif ".2" in c:
	# 			new_cols.append(c.replace(".2", "_2"))
	# 		elif ".3" in c:
	# 			new_cols.append(c.replace(".3", "_3"))
	# 		elif ".4" in c:
	# 			new_cols.append(c.replace(".4", "_4"))
	# 		else:
	# 			new_cols.append(c+"_0")
	# 	df.columns = new_cols

		
	# 	if self.origin == "fetal":
	# 		new_cols = []
	# 		for c in df.columns:
	# 			new_cols.append(c.replace("Fetal_", ""))
	# 		df.columns = new_cols

	# 	df.to_csv("%s/%s_%s.csv"%(self.output, filename, self.origin), sep='\t')

	# 	shutil.rmtree(self.tmp)

	# 	self.multicoverage = df

	# 	print("Computing catlas %s multicoverage DONE!"%origin)


	# def multicoverage_genetics_for_calderon(self, variants=None, filename="multicoverage_calderon"):
		
	# 	if (variants is None) and (self.genetic is None):
	# 		sys.exit("Error, missing genetic file.")
	# 	elif variants is None:
	# 		variants = self.genetic
	# 	elif isinstance(variants, str):
	# 		self.__load_genetic()
	# 		variants = self.genetic

	# 	print("Computing calderon multicoverage might take time ...")

	# 	_, calderon_bigwigs = self.__loading_calderon_data()

	# 	if not os.path.exists(self.tmp):
	# 		os.makedirs(self.tmp)

	# 	bed = self.tmp+"/multicoverage.bed"
	# 	variants.to_csv(bed, sep='\t', index=False, header=None)

	# 	if not os.path.exists(self.output):
	# 		os.makedirs(self.output)

	# 	commands = ["computeMatrix reference-point -S %s -R %s -a 2 -b 3 --missingDataAsZero --binSize 1 -o %s/multicoverage.csv --outFileNameMatrix %s/multicoverage.csv"%(" ".join(calderon_bigwigs), bed, self.output, self.output)]
	# 	subprocess.run(commands, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	# 	df = pd.read_csv("%s/multicoverage.csv"%(self.output), sep='\t', comment='#')
	# 	columns = df.columns[1:]
	# 	df.drop(columns[-1], inplace=True, axis=1)
	# 	df.columns = columns
	# 	index_id = variants['CHR_ID'].astype(str)+"-"+variants['CHR_POS'].astype(str)+"-"+variants['SNPS'].astype(str)
	# 	df.index = index_id
		
	# 	new_cols = []
	# 	for c in df.columns:
	# 		if ".1" in c:
	# 			new_cols.append(c.replace(".1", "_1"))
	# 		elif ".2" in c:
	# 			new_cols.append(c.replace(".2", "_2"))
	# 		elif ".3" in c:
	# 			new_cols.append(c.replace(".3", "_3"))
	# 		elif ".4" in c:
	# 			new_cols.append(c.replace(".4", "_4"))
	# 		else:
	# 			new_cols.append(c+"_0")
	# 	df.columns = new_cols

	# 	df.to_csv("%s/%s.csv"%(self.output, filename), sep='\t')

	# 	shutil.rmtree(self.tmp)

	# 	self.multicoverage = df

	# 	print("Computing calderon multicoverage DONE!")


	# def multicoverage_genetics_for_ludwig2019(self, variants=None, filename="multicoverage_ludwig2019"):
		
	# 	if (variants is None) and (self.genetic is None):
	# 		sys.exit("Error, missing genetic file.")
	# 	elif variants is None:
	# 		variants = self.genetic
	# 	elif isinstance(variants, str):
	# 		self.__load_genetic()
	# 		variants = self.genetic

	# 	print("Computing ludwig2019 multicoverage might take time ...")

	# 	_, ludwig2019_bigwigs = self.__loading_ludwig2019_data()

	# 	if not os.path.exists(self.tmp):
	# 		os.makedirs(self.tmp)

	# 	bed = self.tmp+"/multicoverage.bed"
	# 	variants.to_csv(bed, sep='\t', index=False, header=None)

	# 	if not os.path.exists(self.output):
	# 		os.makedirs(self.output)

	# 	commands = ["computeMatrix reference-point -S %s -R %s -a 2 -b 3 --missingDataAsZero --binSize 1 -o %s/multicoverage.csv --outFileNameMatrix %s/multicoverage.csv"%(" ".join(ludwig2019_bigwigs), bed, self.output, self.output)]
	# 	subprocess.run(commands, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	# 	df = pd.read_csv("%s/multicoverage.csv"%(self.output), sep='\t', comment='#')
	# 	columns = df.columns[1:]
	# 	df.drop(columns[-1], inplace=True, axis=1)
	# 	df.columns = columns
	# 	index_id = variants['CHR_ID'].astype(str)+"-"+variants['CHR_POS'].astype(str)+"-"+variants['SNPS'].astype(str)
	# 	df.index = index_id
		
	# 	new_cols = []
	# 	for c in df.columns:
	# 		if ".1" in c:
	# 			new_cols.append(c.replace(".1", "_1"))
	# 		elif ".2" in c:
	# 			new_cols.append(c.replace(".2", "_2"))
	# 		elif ".3" in c:
	# 			new_cols.append(c.replace(".3", "_3"))
	# 		elif ".4" in c:
	# 			new_cols.append(c.replace(".4", "_4"))
	# 		else:
	# 			new_cols.append(c+"_0")
	# 	df.columns = new_cols

	# 	df.to_csv("%s/%s.csv"%(self.output, filename), sep='\t')

	# 	shutil.rmtree(self.tmp)

	# 	self.multicoverage = df

	# 	print("Computing ludwig2019 multicoverage DONE!")


	# def multicoverage_genetics_for_MPAL_lowGr(self, variants=None, filename="multicoverage_MPAL_lowGr"):
		
	# 	if (variants is None) and (self.genetic is None):
	# 		sys.exit("Error, missing genetic file.")
	# 	elif variants is None:
	# 		variants = self.genetic
	# 	elif isinstance(variants, str):
	# 		self.__load_genetic()
	# 		variants = self.genetic

	# 	print("Computing MPAL_lowGr multicoverage might take time ...")

	# 	_, _, MPAL_lowGr_bigwigs = self.__loading_MPAL_lowGr_data()

	# 	if not os.path.exists(self.tmp):
	# 		os.makedirs(self.tmp)

	# 	bed = self.tmp+"/multicoverage.bed"
	# 	variants.to_csv(bed, sep='\t', index=False, header=None)

	# 	if not os.path.exists(self.output):
	# 		os.makedirs(self.output)

	# 	commands = ["computeMatrix reference-point -S %s -R %s -a 2 -b 3 --missingDataAsZero --binSize 1 -o %s/multicoverage.csv --outFileNameMatrix %s/multicoverage.csv"%(" ".join(MPAL_lowGr_bigwigs), bed, self.output, self.output)]
	# 	subprocess.run(commands, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

	# 	df = pd.read_csv("%s/multicoverage.csv"%(self.output), sep='\t', comment='#')
	# 	columns = df.columns[1:]
	# 	df.drop(columns[-1], inplace=True, axis=1)
	# 	df.columns = columns
	# 	index_id = variants['CHR_ID'].astype(str)+"-"+variants['CHR_POS'].astype(str)+"-"+variants['SNPS'].astype(str)
	# 	df.index = index_id
		
	# 	new_cols = []
	# 	for c in df.columns:
	# 		if ".1" in c:
	# 			new_cols.append(c.replace(".1", "_1"))
	# 		elif ".2" in c:
	# 			new_cols.append(c.replace(".2", "_2"))
	# 		elif ".3" in c:
	# 			new_cols.append(c.replace(".3", "_3"))
	# 		elif ".4" in c:
	# 			new_cols.append(c.replace(".4", "_4"))
	# 		else:
	# 			new_cols.append(c+"_0")
	# 	df.columns = new_cols

	# 	df.to_csv("%s/%s.csv"%(self.output, filename), sep='\t')

	# 	shutil.rmtree(self.tmp)

	# 	self.multicoverage = df

	# 	print("Computing MPAL_lowGr multicoverage DONE!")


	# def plot_genetics_enrichments(self, groupA: list = [], groupB: list = [], multicoverage: str = "", filename: str = "enrichment", show=False):
		
	# 	if multicoverage == "":
	# 		if self.multicoverage.empty:
	# 			sys.exit("Error, multicoverage_genetics_for_X has to be run first or point to the mutlicoverage file created using multicoverage_genetics_for_X.")
	# 	else:
	# 		self.multicoverage = pd.read_csv(multicoverage, sep="\t", index_col=0)

	# 	for ga in groupA:
	# 		if ga+"_2" not in self.multicoverage.columns:
	# 			sys.exit("Error, check 'groupA' entries! It must be a list of one or more cell types.")
	# 	if groupB[0] != "other":
	# 		for gb in groupB:
	# 			if gb+"_2" not in self.multicoverage.columns:
	# 				sys.exit("Error, check 'groupB' entries! It must be ['other'] or a list of one or more cell types.")


	# 	groupA_df_tmp = self.multicoverage[self.multicoverage.columns[self.multicoverage.columns.str.startswith(tuple(groupA))].tolist()]
	# 	groupA_df = pd.DataFrame()
	# 	for suf in ["_0", "_1", "_2", "_3", "_4"]:
	# 		unique_vars = groupA_df_tmp.columns[groupA_df_tmp.columns.str.endswith(suf)]
	# 		groupA_df["GroupA%s"%suf] = groupA_df_tmp[unique_vars].mean(axis=1).tolist()

	# 	if groupB[0] == "other":
	# 		groupB_df_tmp = self.multicoverage[self.multicoverage.columns[~self.multicoverage.columns.str.startswith(tuple(groupA))].tolist()]
	# 	else:
	# 		groupB_df_tmp = self.multicoverage[self.multicoverage.columns[self.multicoverage.columns.str.startswith(tuple(groupB))].tolist()] 
	# 	groupB_df = pd.DataFrame()
	# 	for suf in ["_0", "_1", "_2", "_3", "_4"]:
	# 		unique_vars = groupB_df_tmp.columns[groupB_df_tmp.columns.str.endswith(suf)]
	# 		groupB_df["GroupB%s"%suf] = groupB_df_tmp[unique_vars].mean(axis=1).tolist()

	# 	volcano_df = pd.concat([groupA_df, groupB_df], axis=1)
	# 	column_names = list(volcano_df.columns.values)
	# 	groups = {}
	# 	index = 1
	# 	for column_name in column_names:
	# 		group_name = re.sub('_[0-4]$', '', column_name)

	# 		if group_name in groups:
	# 			groups[group_name].append(index)
	# 		else:
	# 			groups[group_name] = [index];
			
	# 		index = index + 1

	# 	groupnames = list(groups.keys())
	# 	g1_col_indices = groups["GroupA"]
	# 	g2_col_indices = groups["GroupB"]

	# 	g1_values = self.multicoverage.iloc[:, g1_col_indices].to_numpy()
	# 	g2_values = self.multicoverage.iloc[:, g2_col_indices].to_numpy()

	# 	g1_means = g1_values.mean(axis = 1)
	# 	g2_means = g2_values.mean(axis = 1)

	# 	foldchanges = list(np.log2(np.divide(g2_means, g1_means)))
	# 	cov_groupA = groupA_df.sum(axis=1).tolist()

	# 	threshold_cov = np.percentile(cov_groupA, 95)

	# 	colors = []

	# 	for i in range(0, len(foldchanges)):

	# 		if cov_groupA[i] > threshold_cov:

	# 			if foldchanges[i] > 1:
	# 				colors.append('#db3232')
	# 			elif foldchanges[i] < -1:
	# 				colors.append('#3f65d4')
	# 			else:
	# 				colors.append('rgba(150,150,150,0.5)')
	# 		else:
	# 			colors.append('rgba(150,150,150,0.5)')
				
	# 	plot_title = "%s vs other"%(" ".join(groupA))
	# 	x_axis_title = "log2 fold change" 
	# 	y_axis_title = "Coverage"
	# 	point_radius = 8

	# 	fig = go.Figure()

	# 	fig.add_trace(
	# 		go.Scattergl(
	# 			x = foldchanges,
	# 			y = cov_groupA,
	# 			mode = 'markers',
	# 			text = self.multicoverage.index.tolist(),
	# 			hovertemplate ='%{text}: %{x}<br>',
	# 			marker= {
	# 				'color':colors,
	# 				'size':point_radius,
	# 			}
	# 		)
	# 	)
	# 	fig.update_layout(
	# 		title=plot_title,
	# 		xaxis_title= x_axis_title,
	# 		yaxis_title=y_axis_title,
	# 		paper_bgcolor= 'white',
	# 		plot_bgcolor='white',
	# 	)
	# 	fig.update_layout(autosize=False, width=1000, height=800)
	# 	fig.add_hline(y=threshold_cov)
	# 	fig.add_vline(x=1)
	# 	fig.add_vline(x=-1)

	# 	if show:
	# 		fig.show()
	# 	else:
	# 		fig.write_html("%s/%s.html"%(self.output, filename))
	# 		return fig