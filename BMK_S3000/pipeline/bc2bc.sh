#!/bin/bash
#SBATCH --job-name=BMK_bc
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

# Commands to execute
module load Perl/5.38.2-GCCcore-13.3.0

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/exp99_28Feb2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bc2bc/exp99.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/BSTMAtrix/exp100_28Feb2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bc2bc/exp100.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/BSTMAtrix/exp101_28Feb2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bc2bc/exp101.bc

perl /home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/tools/level_bc2bc.pl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/BSTMAtrix/exp102_28Feb2025/barcode.tsv /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz L1 > /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bc2bc/exp102.bc