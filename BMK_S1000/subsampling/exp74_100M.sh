#!/bin/bash
#SBATCH --job-name=subsample_exp74
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00

# Load modules
module load seqtk
module load Perl/5.38.2-GCCcore-13.3.0
module load SAMtools

# Specify arguments
SAMPLE_NAME=exp74
INPUT_FILE_R1=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/BMK_S1000/exp74/4789_Exp74_SC0010_S6_R1.fastq.gz
INPUT_FILE_R2=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/BMK_S1000/exp74/4789_Exp74_SC0010_S6_R2.fastq.gz
OUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/subsampled_100M/FASTQ

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

# Subsample FASTQ files
echo "Subsampling merged FASTQs to $TARGET_READS reads..."
seqtk sample -s$SEED "$INPUT_FILE_R1" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R1_subsampled.fastq.gz"
seqtk sample -s$SEED "$INPUT_FILE_R2" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R2_subsampled.fastq.gz"

# Run BSTMAtrix
CONFIG=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S1000/subsampling/config_${EXP}_100M.txt
echo "Running BSTMatrix..."
/home/VIB.LOCAL/julien.mortier/BMK_software/S1000_pipeline2/BSTMatrix_v2.3.r_English_new/BSTMatrix -c ${CONFIG} -s 0

echo "DONE!"