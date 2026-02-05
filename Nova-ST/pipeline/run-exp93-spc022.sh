#!/bin/bash
#SBATCH --account=lp_edu_ont2024
#SBATCH --clusters=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7:00:00           # hh:mm:ss
#SBATCH --mem=200G
#SBATCH --job-name=novascope93_22

novast_dir=/staging/leuven/stg_00143/julienm/nova-st
container=${novast_dir}/apptainer/novascope_latest.sif

seq1st_dir=/lustre1/project/stg_00003/cbd-bioinf/CBD__Bioinf_Unit__Nova-ST/NovaSeq6000_20241107
job_dir=${novast_dir}/data/2024-07-UTBenchmark-SpC/Exp93-sampleSPC022

singularity exec \
    --bind ${novast_dir}:${novast_dir} \
    --bind ${seq1st_dir}:${seq1st_dir} \
    --bind ${job_dir}:/data \
    $container \
    snakemake \
        -s /app/novascope/NovaScope.smk \
        -d /data \
        --configfile /data/config_job.yaml \
        --cores $SLURM_CPUS_PER_TASK \
        --rerun-incomplete \
        --latency-wait 20 \
        -p

