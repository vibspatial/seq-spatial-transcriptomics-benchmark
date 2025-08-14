#!/bin/bash
#SBATCH --job-name=subsample_exp101
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=48:00:00

# Load modules
module load seqtk
module load Perl/5.38.2-GCCcore-13.3.0
module load SAMtools

# Specify arguments
SAMPLE_NAME=exp101
INPUT_FILE_R1=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/BMK_S3000/sequencing/4941_03_Exp101_SC58_SPC_001_R1.fastq.gz
INPUT_FILE_R2=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/BMK_S3000/sequencing/4941_03_Exp101_SC58_SPC_001_R2.fastq.gz
OUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/subsampled_100M/FASTQ

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

# Subsample FASTQ files
echo "Subsampling merged FASTQs to $TARGET_READS reads..."
seqtk sample -s$SEED "$INPUT_FILE_R1" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R1_subsampled.fastq.gz"
seqtk sample -s$SEED "$INPUT_FILE_R2" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R2_subsampled.fastq.gz"

# Run BSTMAtrix
CONFIG=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S3000/subsampling/config_${SAMPLE_NAME}_100M.txt
echo "Running BSTMatrix..."
/home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix -c ${CONFIG} -s 0

echo "DONE!"