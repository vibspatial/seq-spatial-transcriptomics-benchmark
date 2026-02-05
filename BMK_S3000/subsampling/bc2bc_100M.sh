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

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/subsampled_100M/BSTMAtrix/exp99_100M_11Aug2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/subsampled_100M/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bc2bc/exp99_100M.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/subsampled_100M/BSTMAtrix/exp100_100M_11Aug2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/subsampled_100M/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bc2bc/exp100_100M.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/subsampled_100M/BSTMAtrix/exp101_100M_11Aug2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/subsampled_100M/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bc2bc/exp101_100M.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/subsampled_100M/BSTMAtrix/exp102_100M_11Aug2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/subsampled_100M/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bc2bc/exp102_100M.bc