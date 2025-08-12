#!/bin/bash
#SBATCH --job-name=exp74_BMK
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00

# Commands to execute
module load Perl/5.38.2-GCCcore-13.3.0
module load SAMtools

/home/VIB.LOCAL/julien.mortier/BMK_software/S1000_pipeline2/BSTMatrix_v2.3.r_English_new/BSTMatrix -c /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S1000/config_exp74.txt -s 0
