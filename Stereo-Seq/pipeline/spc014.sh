#!/bin/bash
#SBATCH --job-name=spc014
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

cd /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW

saw count \
    --id=${ID} \
    --sn=${SN} \
    --omics=transcriptomics \
    --kit-version="Stereo-seq T FF V1.3" \
    --sequencing-type="PE75_50+100" \
    --chip-mask=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/mask/${SN}.barcodeToPos.h5 \
    --organism=mouse \
    --tissue=brain \
    --fastqs=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/fastqs \
    --reference=/data/groups/technologies/spatial.catalyst/Julien/SAW/reference/refdata-gex-mm10-2020-A \
    --image-tar=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Stereo-Seq/${ID}/image/${TAR}

echo "Converting .raw.gef ..."
saw convert gef2gem \
   --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.raw.gef \
   --bin-size=1 \
   --gem=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.raw.gem

saw convert gef2h5ad \
   --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.raw.gef \
   --bin-size=1 \
   --h5ad=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.raw.h5ad

BIN_SIZES=( 1 5 10 20 50 100 150 200 )
for BIN_SIZE in "${BIN_SIZES[@]}"
do
    echo "Converting .tissue.gef for bin size ${BIN_SIZE}..."
    saw convert gef2gem \
       --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.tissue.gef \
       --bin-size=${BIN_SIZE} \
       --gem=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.bin${BIN_SIZE}.gem

    saw convert gef2h5ad \
       --gef=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.tissue.gef \
       --bin-size=${BIN_SIZE} \
       --h5ad=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/${ID}/outs/feature_expression/${SN}.bin${BIN_SIZE}.h5ad
done