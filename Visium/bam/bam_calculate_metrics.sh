#!/bin/bash

INPUT_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/testing_arne/output.csv"
TISSUE_BARCODES_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/testing_arne/output_metrics"
TISSUE_AREA=1
STEP_SIZE=10000000 # put to 0 for no sampling


python bam_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE