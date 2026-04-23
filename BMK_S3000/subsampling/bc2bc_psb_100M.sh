#!/bin/bash
#SBATCH --job-name=BMK_bc
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

# Commands to execute
source ~/.bashrc
conda activate BSTMAtrix_env

module load Perl/5.38.2-GCCcore-13.3.0

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/subsampled_100M/BSTMAtrix/psb_100M_18Mar2026/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/subsampled_100M/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/bc2bc/psb_100M.bc