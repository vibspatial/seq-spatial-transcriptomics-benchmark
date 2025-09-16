#!/bin/bash
#SBATCH --job-name=curioseeker-SC12-merge
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=256G
#SBATCH --time=6:00:00

module load Nextflow/22.10.1

export SAMPLE=SC12
export CURIO_SEEKER_SIF=/data/groups/technologies/spatial.catalyst/Julien/curioseeker/seeker-v3.1.0/containers/seeker-v3.1.0.sif
export DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/curioseeker/work/fe/ba101d0c2b4d46e10d3b7006ac96f0/

cd ${DIR}

singularity exec -B /data/groups ${CURIO_SEEKER_SIF} \
seeker \
gen-merged-r1-r2-db \
--sample=${SAMPLE}