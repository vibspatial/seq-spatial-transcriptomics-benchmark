#!/bin/bash
#SBATCH --job-name=exp99_tags
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500G
#SBATCH --time=8:00:00

source ~/.bashrc
conda activate bamtools
module load SAMtools

module load Perl/5.38.2-GCCcore-13.3.0

BC_UMI_READ="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/01.fastq2BcUmi/exp99_28Feb2025.bc_umi_read.tsv.id"
SUBDATA_BARCODES="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/05.AllheStat/BSTViewer_project/subdata/L1_heAuto/barcodes.tsv.gz"
BARCODE="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/exp99_28Feb2025/barcode.tsv"
BARCODES_CLUSTER="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix/level_1/barcodes_cluster.tsv.gz"
INPUT_BAM="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/02.Umi2Gene/exp99_28Feb2025Aligned.out.bam"
OUTPUT_BAM="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output/exp99_tagged.bam"
PREFIX="L1"

# check inputs files
for f in "$BC_UMI_READ" "$SUBDATA_BARCODES" "$BARCODE" "$BARCODES_CLUSTER" "$INPUT_BAM"; do
  if [[ ! -e "$f" ]]; then
    echo "ERROR: required file not found: $f" >&2
    exit 1
  fi
done

perl /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S3000/bam/bam_add_flag.pl \
    $BC_UMI_READ \
    $SUBDATA_BARCODES \
    $BARCODE \
    $BARCODES_CLUSTER \
    $INPUT_BAM \
    $OUTPUT_BAM \
    $PREFIX