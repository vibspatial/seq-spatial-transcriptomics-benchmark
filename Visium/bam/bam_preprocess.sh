#!/bin/bash

BAM_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/SpaceRanger/outs/possorted_genome_bam.bam"
TISSUE_BARCODES_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/harpy/barcodes_in_tissue.txt"
PROBE_SET_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/SpaceRanger/outs/probe_set.csv"
OUTPUT_PATH="/Users/arnedf/VIB/DATA/untargeted_benchmark/SCA013/testing_arne/output.csv"
SIZE=50000000

python bam_preprocess.py \
  --bam_path $BAM_PATH \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --probe_set_path $PROBE_SET_PATH \
  --output_path $OUTPUT_PATH \
  --size $SIZE
