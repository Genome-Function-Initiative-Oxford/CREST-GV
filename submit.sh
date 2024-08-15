#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=CREST-GV
#SBATCH --ntasks=4
#SBATCH --mem=90G
#SBATCH --mail-user=simone.riva@imm.ox.ac.uk
#SBATCH --time=06-23:59:59
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

source /project/Wellcome_Discovery/sriva/mambaforge/bin/activate crestgv

python crestgv.py -g /project/Wellcome_Discovery/sriva/data_and_metadata/231126_GENOMICCupdated.allLD.bed -o CREST-GV_test_all_run -ra True -g25k True -ra True