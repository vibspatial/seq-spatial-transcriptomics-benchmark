#!/bin/bash
#SBATCH --job-name=SCA013
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=8:00:00

# Specify arguments
SAMPLE_NAME="SCA013"

## Run Space Ranger
spaceranger count --id="$SAMPLE_NAME" \
    --transcriptome=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A \
    --probe-set=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/probe_sets/mouse_transcriptome_v2/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv \
    --fastqs=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA013/FASTQ/0-0-Exp4863_2-SCA,/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA013/FASTQ/0-0-Exp4972-POOL_SCC133-135,/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA013/FASTQ/0-0-Exp5098-POOL_SCC166-167 \
    --sample="${SAMPLE_NAME}_GEX" \
    --cytaimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA013/IMAGES/SCA013_CytAssist.tif \
    --image=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA013/IMAGES/SCA013_HE_hires.tif \
    --slide=V44Y13-296 \
    --area=D1 \
    --create-bam=true \
    --output-dir=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/SpaceRanger
