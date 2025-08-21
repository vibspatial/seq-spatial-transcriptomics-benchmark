#!/bin/bash
#SBATCH --job-name=subsample_SC11
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=256G
#SBATCH --time=6:00:00

# Load modules
module load Nextflow/22.10.1
module load seqtk

# Specify arguments
export SAMPLE_NAME=SC11
INPUT_FILE_R1=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Curio_seeker/sequencing/${SAMPLE_NAME}_R1.fastq.gz
INPUT_FILE_R2=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Curio_seeker/sequencing/${SAMPLE_NAME}_R2.fastq.gz
OUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/${SAMPLE_NAME}/subsampled_100M/FASTQ

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

# Subsample FASTQ files
echo "Subsampling merged FASTQs to $TARGET_READS reads..."
seqtk sample -s$SEED "$INPUT_FILE_R1" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R1_subsampled.fastq.gz"
seqtk sample -s$SEED "$INPUT_FILE_R2" $TARGET_READS | gzip > "${OUT_DIR}/${SAMPLE_NAME}_R2_subsampled.fastq.gz"

# Run curioseeker

export CURIO_SEEKER=/data/groups/technologies/spatial.catalyst/Julien/curioseeker/curioseeker-v3.0.0
export OUTPUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/${SAMPLE_NAME}/subsampled_100M/curioseeker
export SAMPLESHEET_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker/subsampling
export GENOMES_DIR=/data/groups/technologies/spatial.catalyst/Julien/curioseeker_references

mkdir -p ${OUTPUT_DIR}

nextflow run ${CURIO_SEEKER}/main.nf \
        --input ${SAMPLESHEET_DIR}/samplesheet_${SAMPLE_NAME}_100M.csv \
        --outdir ${OUTPUT_DIR}/results/ \
        -work-dir ${OUTPUT_DIR}/work/ \
        --igenomes_base ${GENOMES_DIR} \
        -resume \
        -profile singularity