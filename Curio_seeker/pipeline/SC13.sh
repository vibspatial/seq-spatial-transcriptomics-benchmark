#!/bin/bash
#SBATCH --job-name=curioseeker-SC13
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=256G
#SBATCH --time=6:00:00

module load Nextflow/22.10.1

export SAMPLE=SC13

export CURIO_SEEKER=/data/groups/technologies/spatial.catalyst/Julien/curioseeker/curioseeker-v3.0.0
export OUTPUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/${SAMPLE}/curioseeker
export SAMPLESHEET_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker/pipeline
export GENOMES_DIR=/data/groups/technologies/spatial.catalyst/Julien/curioseeker_references

mkdir -p ${OUTPUT_DIR}

nextflow run ${CURIO_SEEKER}/main.nf \
        --input ${SAMPLESHEET_DIR}/samplesheet_${SAMPLE}.csv \
        --outdir ${OUTPUT_DIR}/results/ \
        -work-dir ${OUTPUT_DIR}/work/ \
        --igenomes_base ${GENOMES_DIR} \
        -resume \
        -profile singularity