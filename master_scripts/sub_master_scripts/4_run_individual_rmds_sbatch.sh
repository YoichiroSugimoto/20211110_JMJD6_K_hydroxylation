#!/bin/bash
#SBATCH --job-name=4_JMJD6_hydroxyK
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/4_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate hydroxylation_by_JMJD6

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/R/j4-de-novo-data-analysis

Rscript ../run_rmd.R j4-de-novo-data-analysis.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
