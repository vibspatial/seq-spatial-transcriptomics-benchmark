#!/bin/bash
#SBATCH --job-name=curioseeker-SC11
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=75G
#SBATCH --time=2:00:00

module load Nextflow/22.10.1

SAMPLE=SC11

CURIO_SEEKER=/home/VIB.LOCAL/julien.mortier/curioseeker/curioseeker-v3.0.0
OUTPUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/${SAMPLE}/curioseeker
SAMPLESHEET_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker
GENOMES_DIR=/data/groups/technologies/spatial.catalyst/Julien/curioseeker_references

mkdir -p ${OUTPUT_DIR}

nextflow run ${CURIO_SEEKER}/main.nf \
        --input ${SAMPLESHEET_DIR}/samplesheet_${SAMPLE}.csv \
        --outdir ${OUTPUT_DIR} \
        -work-dir ${OUTPUT_DIR}/work/ \
        --igenomes_base ${GENOMES_DIR} \
        -resume \
        -profile singularity