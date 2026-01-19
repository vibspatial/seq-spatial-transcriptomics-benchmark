#!/bin/bash
#SBATCH --job-name=spc014_100M
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=julien.mortier@vib.be

set -euo pipefail

ID=spc014
SN=Y40048L4
TAR=Y40048L4_SC_20251210_115137_4.2.1.tar.gz

# Load modules
module load seqtk

## Subsample FASTQ files

# Specify arguments
INPUT_DIRS=("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/fastqs")
OUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/subsampled_100M/FASTQ"

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

echo "Subsampling FASTQ files for sample $ID..."

# Gather all R1 and R2 files
R1_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${ID}*_R1_*.fastq.gz; done)
R2_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${ID}*_R2_*.fastq.gz; done)

# Count reads in each R1 file (4 lines per read)
TOTAL_READS=0
declare -A FILE_READS
for f in $R1_FILES; do
    echo "Counting reads for R1 file: $f"
    COUNT=$(zcat "$f" | wc -l)
    COUNT=$((COUNT / 4))
    FILE_READS["$f"]=$COUNT
    echo "File: $f, Reads: $COUNT"
    TOTAL_READS=$((TOTAL_READS + COUNT))
done

echo "Total reads across all R1 files: $TOTAL_READS"

# Subsample each file proportionally
for f in $R1_FILES; do
    R2_FILE=${f/_R1_/_R2_}
    COUNT=${FILE_READS["$f"]}
    KEEP=$((COUNT * TARGET_READS / TOTAL_READS))

    RUN_NAME=$(basename "$(dirname "$f")")
    RUN_DIR="${OUT_DIR}/${RUN_NAME}"
    mkdir -p "$RUN_DIR"

    OUT_R1="${RUN_DIR}/$(basename "$f")"
    OUT_R2="${RUN_DIR}/$(basename "$R2_FILE")"

    echo "Subsampling $f and $R2_FILE to $KEEP reads each..."
    seqtk sample -s$SEED "$f" $KEEP | gzip > "$OUT_R1"
    seqtk sample -s$SEED "$R2_FILE" $KEEP | gzip > "$OUT_R2"
done

# Count reads after subsampling
READS_KEPT=$(find "$OUT_DIR" -type f -name "*_R1_*.fastq.gz" -exec zcat {} + | wc -l)
READS_KEPT=$((READS_KEPT / 4))
echo "Total reads after subsampling: $READS_KEPT"

## Run SAW
cd /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW

saw count \
    --id=${ID} \
    --sn=${SN} \
    --omics=transcriptomics \
    --kit-version="Stereo-seq T FF V1.3" \
    --sequencing-type="PE75_50+100" \
    --chip-mask=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/mask/${SN}.barcodeToPos.h5 \
    --organism=mouse \
    --tissue=brain \
    --fastqs=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/FASTQ \
    --reference=/data/groups/technologies/spatial.catalyst/Julien/SAW/reference/refdata-gex-mm10-2020-A \
    --image-tar=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/image/${TAR}

echo "Converting .raw.gef ..."
saw convert gef2gem \
   --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.raw.gef \
   --bin-size=1 \
   --gem=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.raw.gem

saw convert gef2h5ad \
   --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.raw.gef \
   --bin-size=1 \
   --h5ad=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.raw.h5ad

BIN_SIZES=( 1 5 10 20 50 100 150 200 )
for BIN_SIZE in "${BIN_SIZES[@]}"
do
    echo "Converting .tissue.gef for bin size ${BIN_SIZE}..."
    saw convert gef2gem \
       --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.tissue.gef \
       --bin-size=${BIN_SIZE} \
       --gem=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.bin${BIN_SIZE}.gem

    saw convert gef2h5ad \
       --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.tissue.gef \
       --bin-size=${BIN_SIZE} \
       --h5ad=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/subsampled_100M/SAW/${ID}/outs/feature_expression/${SN}.bin${BIN_SIZE}.h5ad
done