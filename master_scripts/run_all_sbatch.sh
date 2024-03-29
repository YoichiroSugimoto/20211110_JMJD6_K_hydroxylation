#!/bin/bash
#SBATCH --job-name=HP5_run_all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:0
#SBATCH --mem=8G
#SBATCH --partition=cpu

## Run with the following commands
# sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/run_all_sbatch.sh

jid00=$(sbatch /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/0_run_individual_rmds_sbatch.sh)
jid0=$(echo ${jid00##* })

jid01=$(sbatch --dependency=afterok:$jid0 /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/1_run_individual_rmds_sbatch.sh)
jid1=$(echo ${jid01##* })

jid02=$(sbatch --dependency=afterok:$jid1 /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/2_run_individual_rmds_sbatch.sh)
jid2=$(echo ${jid02##* })

jid03=$(sbatch --dependency=afterok:$jid2 /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/3_run_individual_rmds_sbatch.sh)
jid3=$(echo ${jid03##* })

jid04=$(sbatch --dependency=afterok:$jid3 /camp/home/sugimoy/home/CAMP_HPC/projects/20211110_JMJD6_K_hydroxylation/master_scripts/sub_master_scripts/4_run_individual_rmds_sbatch.sh)
jid3=$(echo ${jid04##* })