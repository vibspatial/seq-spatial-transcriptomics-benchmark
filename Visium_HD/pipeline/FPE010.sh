#!/bin/bash
#SBATCH --job-name=FPE010
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256G
#SBATCH --time=16:00:00

spaceranger count --id=FPE010 \
    --transcriptome=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/IRGSP.49 \
    --fastqs=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/6/FASTQ/0-0-Exp5193-POOL_SCC190,/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/6/FASTQ/0-0-Exp5210-POOL_SCC193-194 \
    --sample=FPE010_cDNA \
    --cytaimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/6/IMAGES/FPE010_CytAssist.tif \
    --darkimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/6/IMAGES/FPE010_fluorescense_grey_hires.tif \
    --loupe-alignment=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/6/IMAGES/FPE010_alignment_JM.json \
    --slide=H1-DWTNKV8 \
    --area=D1 \
    --create-bam=true \
    --output-dir=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/FPE010/SpaceRanger \
    --custom-bin-size=32