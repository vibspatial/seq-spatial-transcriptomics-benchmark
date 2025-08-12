#!/bin/bash
#SBATCH --job-name=SCA010
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256G
#SBATCH --time=48:00:00

spaceranger count --id=SCA010 \
    --transcriptome=/data/projects/s00/spaceranger/references/refdata-gex-GRCh38-2020-A \
    --probe-set=/home/VIB.LOCAL/julien.mortier/spaceranger/spaceranger-3.1.2/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --fastqs=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/FASTQ/0-0-Exp4825-POOL_SCC077 \
    --sample=SCA010_GEX \
    --cytaimage=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/IMAGES/SCA010_CytAssist.tif \
    --image=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/IMAGES/SCA010_HE_hires.tif \
    --slide=H1-3WGDYRZ \
    --area=D1 \
    --create-bam=true \
    --output-dir=/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010 \
    --custom-bin-size=32