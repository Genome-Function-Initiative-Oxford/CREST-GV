import warnings
warnings.filterwarnings('ignore')

import sys, os, argparse
import pandas as pd

sys.path.append('crestgv/')
from crestgv import crestgv
from _version import __version__

collection_name_dict = {'cad'                     : 'CAD',
                        'calderon'                : 'calderon',
                        'catlas_fetal'            : 'catlas_fetal',
                        'catlas_adult'            : 'catlas_adult',
                        'erythoid_d7_d10_d13_d17' : 'Days7_10_13_17',
                        'h1_hescs'                : 'H1_hESCs',
                        'immune_cell'             : 'immune_cell',
                        'ludwig2019'              : 'ludwig2019',
                        'mpal'                    : 'MPAL',
                        'pancreatic_pbmc'         : 'pancreatic_pbmc',
                        'super_pbmc'              : 'super_PBMC'
                       }


def getArgs():

    parser = argparse.ArgumentParser(description="Enrichment Score for VARiants", add_help=False)

    # The required setting
    parser.add_argument("-g", "--genetic", help="Path to genetic file.", nargs="?", required=True, type=str)

    # Optional settings
    parser.add_argument("-nof",  "--number_of_folds", help="Number of folds to create backgound using the 1000genomes.", nargs="?", default=5, type=int)
    parser.add_argument("-o",    "--output",          help="Directory where to save the scores.", nargs="?", default="output", type=str)
    parser.add_argument("-gb",   "--genome",          help="Genome to use, available 'hg19' and 'hg38'.", nargs="?", default="hg38", type=str)    
    parser.add_argument("-ng",   "--min_number_genetics",          help="Subset number for genetic to query. Values allow in range(100, 1000).", nargs="?", default=100, type=int)
    parser.add_argument("-s",    "--seed",            help="Seed for reproducibility, shuffle 1000genomes excluded.", nargs="?", default=42, type=int)

    # Optional settings (choose one or the other)
    parser.add_argument("-cn",   "--collection_name",          help="Data collection name, available 'cad', 'calderon', 'catlas_fetal', 'catlas_adult', 'erythoid_d7_d10_d13_d17', 'h1_hescs', 'immune_cell', 'ludwig2019', 'mpal', 'pancreatic_pbmc', and 'super_pbmc'.", nargs="?", default="", type=str)
    parser.add_argument("-incp", "--in_house_collection_path", help="In house data collection path (<path-to-directory>/<collection-name>).", nargs="?", default="", type=str)

    # Optional settings when list of variants is less than 100
    parser.add_argument("-lng", "--lessNG", help="Boolean variable to force the software to run also with less than 100 variants per file.", nargs="?", default=False, type=bool)

    # Optional settings when list of variants is greater than 25,000
    parser.add_argument("-g25k", "--greater25k", help="Boolean variable to check if you want to run CREST-GV on more than 25k variants.", nargs="?", default=False, type=bool)

    # Optional settings for running only coverage calculation.
    parser.add_argument("-gc", "--get_coverage", help="Run only coverage calculation. Enrichment score will be ignored.", nargs="?", default=True, type=bool)

    # Optional settings for running CREST-GV for all data collection.
    parser.add_argument("-ra", "--run_all", help="Run CREST-GV for all data collection. If -cn or -incp are set, they will be ignored.", nargs="?", default=False, type=bool)

    args = vars(parser.parse_args())

    return args


def main():
    args = getArgs()

    if args["run_all"]:
        for collection_i in list(collection_name_dict.keys()):
            cgv = crestgv(genetic=args["genetic"], 
                          number_of_folds=args["number_of_folds"], 
                          output=args["output"]+os.sep+collection_i, 
                          genome=args["genome"], 
                          min_number_genetics=args["min_number_genetics"],
                          seed=args["seed"], 
                          collection_name=collection_i,
                          in_house_collection_path=args["in_house_collection_path"]
                        )
            if args["get_coverage"]:
                _ = cgv.calculate_enrichment_score(lessNG=args["lessNG"], greater25k=args["greater25k"])
            _ = cgv.get_coverage()
    else:
        cgv = crestgv(genetic=args["genetic"], 
                      number_of_folds=args["number_of_folds"], 
                      output=args["output"], 
                      genome=args["genome"], 
                      min_number_genetics=args["min_number_genetics"],
                      seed=args["seed"], 
                      collection_name=args["collection_name"],
                      in_house_collection_path=args["in_house_collection_path"]
                  )
        if args["get_coverage"]:
            _ = cgv.calculate_enrichment_score(lessNG=args["lessNG"], greater25k=args["greater25k"])
        _ = cgv.get_coverage()


if __name__ == "__main__":

    main()