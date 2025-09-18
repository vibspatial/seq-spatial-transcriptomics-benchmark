#!/bin/bash

INPUT_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/testing_arne/output_all.csv"
TISSUE_BARCODES_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/testing_arne/output_metrics"
TISSUE_AREA=1
SAMPLE_SIZE=100000000 # put to 0 for no sampling


python bam_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --sample_size $SAMPLE_SIZE