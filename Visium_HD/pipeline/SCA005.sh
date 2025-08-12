#!/bin/bash
#SBATCH --job-name=SCA005
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256G
#SBATCH --time=48:00:00

spaceranger count --id=SCA005 \
    --transcriptome=/data/projects/s00/spaceranger/references/refdata-gex-mm10-2020-A \
    --probe-set=/home/VIB.LOCAL/julien.mortier/spaceranger/spaceranger-3.1.2/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv \
    --fastqs=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/2/FASTQ/0-0-Exp4799-POOL_SCC071,/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/2/FASTQ/0-0-Exp4823-POOL_SCC076 \
    --sample=SCA005_GEX \
    --cytaimage=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/2/IMAGES/SCA005_CytAssist.tif \
    --image=/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/2/IMAGES/SCA005_HE_hires.tif \
    --slide=H1-NV2FPB2 \
    --area=A1 \
    --create-bam=true \
    --output-dir=/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005 \
    --custom-bin-size=32