#!/bin/bash
#SBATCH --job-name=curioseeker-SC11
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=75G
#SBATCH --time=2:00:00

# Load the Nextflow module on the VIB datacore cluster
module load Nextflow/22.10.1

export SAMPLE=SC11

export CURIO_ROOT=/data/projects/s02/Frank/curio
export CURIO_SEEKER=/home/VIB.LOCAL/julien.mortier/curioseeker/curioseeker-v3.0.0
export OUTPUT_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/${SAMPLE}/curioseeker
export SAMPLESHEET_DIR=`realpath ../samplesheets`

mkdir -p ${OUTPUT_DIR}

nextflow run ${CURIO_SEEKER}/main.nf \
        --input ${SAMPLESHEET_DIR}/samplesheet_${SAMPLE}.csv \
        --outdir ${OUTPUT_DIR}/results/ \
        -work-dir ${OUTPUT_DIR}/work/ \
        --igenomes_base ${CURIO_ROOT}/genomes/ \
        -resume \
        -profile singularity