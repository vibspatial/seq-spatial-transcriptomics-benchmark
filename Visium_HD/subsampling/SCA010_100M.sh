#!/bin/bash
#SBATCH --job-name=subsample_SCA010
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=256G
#SBATCH --time=12:00:00

# Load modules
module load seqtk

## Subsample FASTQ files

# Specify arguments
SAMPLE_NAME="SCA010"
INPUT_DIRS=("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/FASTQ/0-0-Exp4825-POOL_SCC077")
OUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/subsampled_100M/FASTQ"

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

# Gather all R1 and R2 files
R1_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${SAMPLE_NAME}*_R1_*.fastq.gz; done)
R2_FILES=$(for d in "${INPUT_DIRS[@]}"; do ls ${d}/${SAMPLE_NAME}*_R2_*.fastq.gz; done)

# Count reads in each R1 file (4 lines per read)
echo "Counting reads per R1 file..."
TOTAL_READS=0
declare -A FILE_READS
for f in $R1_FILES; do
    COUNT=$(zcat "$f" | wc -l)
    COUNT=$((COUNT / 4))
    FILE_READS["$f"]=$COUNT
    TOTAL_READS=$((TOTAL_READS + COUNT))
done

echo "Total reads across all R1 files: $TOTAL_READS"

# Subsample each file proportionally
for f in $R1_FILES; do
    R2_FILE=${f/_R1_/_R2_}
    COUNT=${FILE_READS["$f"]}
    KEEP=$((COUNT * TARGET_READS / TOTAL_READS))

    OUT_R1="${OUT_DIR}/$(basename "$f" .fastq.gz)_subsampled.fastq.gz"
    OUT_R2="${OUT_DIR}/$(basename "$R2_FILE" .fastq.gz)_subsampled.fastq.gz"

    echo "Subsampling $f and $R2_FILE to $KEEP reads each..."
    seqtk sample -s$SEED "$f" $KEEP | gzip > "$OUT_R1"
    seqtk sample -s$SEED "$R2_FILE" $KEEP | gzip > "$OUT_R2"
done

## Run Space Ranger
spaceranger count --id="$SAMPLE_NAME" \
    --transcriptome=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A \
    --probe-set=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/probe_sets/human_transcriptome_v2/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
    --fastqs="$OUT_DIR" \
    --sample="$SAMPLE_NAME" \
    --cytaimage=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/IMAGES/SCA010_CytAssist.tif \
    --image=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/3/IMAGES/SCA010_HE_hires.tif \
    --slide=H1-3WGDYRZ \
    --area=D1 \
    --create-bam=true \
    --output-dir=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/subsampled_100M/SpaceRanger \
    --custom-bin-size=32
