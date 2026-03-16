#!/bin/bash
#SBATCH --job-name=SCA014_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=400G
#SBATCH --time=200:00:00
#SBATCH --output=SCA014_metrics_%j.out
#SBATCH --error=SCA014_metrics_%j.err

source /home/VIB.LOCAL/arne.defauw/code/2024-07-UTBenchmark-SpC/.venv_bam_env/bin/activate

INPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output/possorted_genome_bam_preprocessed.csv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output/output_metrics" # directory will be created
TISSUE_AREA=28.03
STEP_SIZE=10000000 # put to 0 for no sampling

python bam_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE
