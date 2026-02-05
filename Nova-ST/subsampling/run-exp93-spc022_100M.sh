#!/bin/bash
#SBATCH --account=lp_edu_ont2024
#SBATCH --clusters=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00           # hh:mm:ss
#SBATCH --mem=200G
#SBATCH --job-name=novascope93_22_100M

export LMOD_CACHED_LOADS=yes
export LMOD_RC=/data/leuven/software/biomed/${VSC_INSTITUTE_CLUSTER}_${VSC_ARCH_LOCAL}_${VSC_OS_LOCAL}/lmodrc.lua

echo LMOD_RC=$LMOD_RC
echo LMOD_CACHED_LOADS=$LMOD_CACHED_LOADS

module use /data/leuven/software/biomed/${VSC_INSTITUTE_CLUSTER}_${VSC_ARCH_LOCAL}_${VSC_OS_LOCAL}/20*/modules/all

module load cluster/wice/batch_icelake

module load seqtk

## Subsample FASTQ files
SAMPLE_NAME="Exp93_SPC022"
INPUT_DIRS=("/staging/leuven/stg_00143/julienm/nova-st/data/2024-07-UTBenchmark-SpC/Exp93-sampleSPC022/input/5147" "/staging/leuven/stg_00143/julienm/nova-st/data/2024-07-UTBenchmark-SpC/Exp93-sampleSPC022/input/5187")
OUT_DIR="/staging/leuven/stg_00143/julienm/nova-st/data/2024-07-UTBenchmark-SpC/Exp93-sampleSPC022_100M/input/"

TARGET_READS=100000000
SEED=42

mkdir -p "$OUT_DIR"

echo "Subsampling FASTQ files for sample $SAMPLE_NAME..."

# Gather all R1 and R2 files
R1_FILES=$(for d in "${INPUT_DIRS[@]}"; do
    ls ${d}/*${SAMPLE_NAME}*_R1_*.fastq.gz | grep -v "_subset.fastq.gz"
done)

R2_FILES=$(for d in "${INPUT_DIRS[@]}"; do
    ls ${d}/*${SAMPLE_NAME}*_R2_*.fastq.gz | grep -v "_subset.fastq.gz"
done)

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

## Run NovaScope
novast_dir=/staging/leuven/stg_00143/julienm/nova-st
container=${novast_dir}/apptainer/novascope_latest.sif

seq1st_dir=/lustre1/project/stg_00003/cbd-bioinf/CBD__Bioinf_Unit__Nova-ST/NovaSeq6000_20241107
job_dir=${novast_dir}/data/2024-07-UTBenchmark-SpC/Exp93-sampleSPC022_100M

singularity exec \
    --bind ${novast_dir}:${novast_dir} \
    --bind ${seq1st_dir}:${seq1st_dir} \
    --bind ${job_dir}:/data \
    $container \
    snakemake \
        -s /app/novascope/NovaScope.smk \
        -d /data \
        --configfile /data/config_job.yaml \
        --cores $SLURM_CPUS_PER_TASK \
        --rerun-incomplete \
        --latency-wait 20 \
        -p
