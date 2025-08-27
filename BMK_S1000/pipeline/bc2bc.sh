#!/bin/bash
#SBATCH --job-name=BMK_bc
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

# Commands to execute
module load Perl/5.38.2-GCCcore-13.3.0

perl /home/VIB.LOCAL/julien.mortier/BMK_software/S1000_pipeline2/BSTMatrix_v2.3.r_English_new/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/BSTMAtrix/exp74_5Jun2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/bc2bc/exp74.bc