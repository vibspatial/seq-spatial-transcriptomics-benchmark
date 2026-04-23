#!/bin/bash
#SBATCH --job-name=SCA023
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256G
#SBATCH --time=16:00:00

spaceranger count --id=SCA023 \
    --transcriptome=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A \
    --fastqs=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/7/FASTQ/0-0-Exp5236-POOL_SCC199 \
    --sample=SCA023_cDNA \
    --cytaimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/7/IMAGES/SCA023_CytAssist.tif \
    --image=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/7/IMAGES/SCA023_HE_hires.tif \
    --slide=H1-YX2GQTZ \
    --area=A1 \
    --create-bam=true \
    --output-dir=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA023/SpaceRanger \
    --custom-bin-size=32