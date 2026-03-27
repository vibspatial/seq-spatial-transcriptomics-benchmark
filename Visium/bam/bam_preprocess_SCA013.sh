#!/bin/bash
#SBATCH --job-name=SCA013_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=SCA013_bam_%j.out
#SBATCH --error=SCA013_bam_%j.err

BAM_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/SpaceRanger/outs/possorted_genome_bam.bam"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/harpy/barcodes_in_tissue.txt"
PROBE_SET_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/SpaceRanger/outs/probe_set.csv"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output/possorted_genome_bam_preprocessed.csv"

#SIZE=000000

python bam_preprocess.py \
  --bam_path $BAM_PATH \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --probe_set_path $PROBE_SET_PATH \
  --output_path $OUTPUT_PATH
#  --size $SIZE
