#!/bin/bash
#SBATCH --job-name=2_JMJD6_IP
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:0
#SBATCH --mem=80G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20210408_JMJD6_IP/master_scripts/sub_master_scripts/2_run_individual_rmds_sbatch.sh

source /camp/home/sugimoy/.bashrc
conda activate HP5_extension

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210408_JMJD6_IP/R/j-a-data-preprocessing

Rscript ../run_rmd.R j-a-1-protein-feature-extraction.rmd
Rscript ../run_rmd.R j-a-2-additional-protein-feature-extraction.rmd
Rscript ../run_rmd.R j-a-3-merge-all-data.rmd

cd /camp/home/sugimoy/home/CAMP_HPC/projects/20210408_JMJD6_IP/R/j-d-de-novo-data-analysis

Rscript ../run_rmd.R jd1-de-novo-data-analysis.rmd

## Template for sbatch
## Rscript ../run_rmd.R xx.rmd
