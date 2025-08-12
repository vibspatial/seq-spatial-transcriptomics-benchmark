#!/bin/bash
#SBATCH --job-name=exp66_BMK
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=128G

# Commands to execute
module load Perl/5.38.2-GCCcore-13.3.0
module load SAMtools

/home/VIB.LOCAL/julien.mortier/BMK_software/S1000_pipeline2/BSTMatrix_v2.3.r_English_new/BSTMatrix -c /data/projects/s00/2024-07-UTBenchmark-SpC/config_exp66.txt -s 0
