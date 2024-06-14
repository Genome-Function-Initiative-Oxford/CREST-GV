import warnings
warnings.filterwarnings('ignore')

import sys, os, argparse
import pandas as pd

sys.path.append('ensgen/')
from ensgen import ensgen
from _version import __version__


def getArgs():

    parser = argparse.ArgumentParser(description="Enrichment score for genetics", add_help=False)

    # The required setting
    parser.add_argument("-g", "--genetic", help="Genetics file", nargs="?", required=True, type=str)

    # Optional settings
    parser.add_argument("-t",  "--tmp",             help="Temporary folder", nargs="?", default="tmp", type=str)
    parser.add_argument("-nf", "--number_of_folds", help="Number of folds for background generation", nargs="?", default=5, type=int)
    parser.add_argument("-f",  "--folds",           help="Folder where to store background", nargs="?", default="folds", type=str)
    parser.add_argument("-o",  "--output",          help="Output folder", nargs="?", default="output", type=str)
    parser.add_argument("-gb", "--genome",          help="Genome build version", nargs="?", default="hg38", type=str)
    parser.add_argument("-s",  "--seed",            help="Seed", nargs="?", default=42, type=int)

    # Settings without a parameter value
    parser.add_argument("--catlas-fetal",    help='Run genetic on CATLAS fetal', nargs="?", default=False, type=bool)
    parser.add_argument("--catlas-adult",    help='Run genetic on CATLAS adult', nargs="?", default=False, type=bool)
    parser.add_argument("--calderon",        help='Run genetic on Calderon', nargs="?", default=False, type=bool)
    parser.add_argument("--ludwig2019",      help='Run genetic on Calderon', nargs="?", default=False, type=bool)
    parser.add_argument("--MPAL_lowGr",      help='Run genetic on Calderon', nargs="?", default=False, type=bool)
    parser.add_argument("--super_PBMC",      help='Run genetic on super PBMC', nargs="?", default=False, type=bool)
    parser.add_argument("--Days7_10_13_17",  help='Run genetic on Day 7, 10, 13, and 17', nargs="?", default=False, type=bool)
    parser.add_argument("--immune_cell",     help='Run genetic on Immune Cell', nargs="?", default=False, type=bool)
    parser.add_argument("--pancreatic_pbmc", help='Run genetic on Pancreatic PBMC', nargs="?", default=False, type=bool)
    parser.add_argument("--H1_hESCs",        help='Run genetic on H1 hESCs', nargs="?", default=False, type=bool)
    parser.add_argument("--CAD",             help='Run genetic on CAD', nargs="?", default=False, type=bool)

    args = vars(parser.parse_args())

    return args



def main():
    args = getArgs()

    if args["catlas-adult"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_catlas(origin='adult')#, umap=True
        es.plot_genetics_for_catlas(show=False)
        # es.multicoverage_genetics_for_catlas(origin='adult')

    if args["catlas-fetal"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_catlas(origin='fetal')#, umap=True
        es.plot_genetics_for_catlas(show=False)
        # es.multicoverage_genetics_for_catlas(origin='fetal')

    if args["calderon"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_calderon()
        es.plot_genetics_for_calderon(show=False)
        # es.multicoverage_genetics_for_calderon()

    if args["ludwig2019"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_ludwig2019()
        es.plot_genetics_for_ludwig2019(show=False)
        # es.multicoverage_genetics_for_ludwig2019()

    if args["MPAL_lowGr"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_MPAL_lowGr()
        es.plot_genetics_for_MPAL_lowGr(show=False)#, umap=True
        # es.multicoverage_genetics_for_MPAL_lowGr()

    if args["super_PBMC"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_super_PBMC()
        es.plot_genetics_for_super_PBMC(show=False, umap=True)

    if args["Days7_10_13_17"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_Days7_10_13_17()
        es.plot_genetics_for_Days7_10_13_17(show=False)

    if args["immune_cell"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_immune_cell()
        es.plot_genetics_for_immune_cell(show=False)

    if args["pancreatic_pbmc"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_pancreatic_pbmc()
        es.plot_genetics_for_pancreatic_pbmc(show=False)

    if args["H1_hESCs"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_H1_hESCs()
        es.plot_genetics_for_H1_hESCs(show=False)

    if args["CAD"]:
        es = ensgen(
                    genetic=args["genetic"], 
                    tmp=args['tmp'], 
                    number_of_folds=args['number_of_folds'],
                    folds=args['folds'],
                    output=args['output'],
                    genome=args['genome'],
                    seed=args['seed'])
        df = es.process_genetics_for_CAD()
        es.plot_genetics_for_CAD(show=False, umap=True)

if __name__ == "__main__":

    main()