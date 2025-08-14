#!/bin/bash
#SBATCH --job-name=subsample_SCA020
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=12:00:00

# Load modules
module load seqtk

## Subsample FASTQ files

# Specify arguments
SAMPLE_NAME="SCA020"
INPUT_DIRS=("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/4/FASTQ/0-0-Exp4961-SCC125-126" "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/4/FASTQ/0-0-Exp4988-SCC137-138")
OUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/subsampled_100M/FASTQ"

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

echo "Subsampling FASTQ files for sample $SAMPLE_NAME..."

# Gather all R1 and R2 files
R1_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${SAMPLE_NAME}*_R1_*.fastq.gz; done)
R2_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${SAMPLE_NAME}*_R2_*.fastq.gz; done)

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

## Run Space Ranger
spaceranger count --id="$SAMPLE_NAME" \
    --transcriptome=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A \
    --probe-set=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/probe_sets/human_transcriptome_v2/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --fastqs="$OUT_DIR" \
    --sample="${SAMPLE_NAME}_GEX" \
    --cytaimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/4/IMAGES/SCA020_CytAssist.tif \
    --image=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/4/IMAGES/SCA020_HE_hires.tif \
    --slide=H1-MQY8KTW \
    --area=D1 \
    --create-bam=true \
    --output-dir=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/subsampled_100M/SpaceRanger \
    --custom-bin-size=32
