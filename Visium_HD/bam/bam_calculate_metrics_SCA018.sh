#!/bin/bash
#SBATCH --job-name=SCA018_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=400G
#SBATCH --time=70:00:00
#SBATCH --output=SCA018_metrics_%j.out
#SBATCH --error=SCA018_metrics_%j.err


INPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output/possorted_genome_bam_preprocessed.csv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output/output_metrics" # directory will be created
TISSUE_AREA=21.19
STEP_SIZE=10000000 # put to 0 for no sampling

python /home/VIB.LOCAL/arne.defauw/code/2024-07-UTBenchmark-SpC/Visium/bam/bam_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE
