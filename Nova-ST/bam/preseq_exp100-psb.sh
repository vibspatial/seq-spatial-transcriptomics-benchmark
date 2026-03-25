#!/bin/bash
#SBATCH --job-name=novast_psb_preseq
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

source ~/.bashrc
conda activate preseq
set -euo pipefail

BAM_OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-samplePSB/bam_output"

echo -e "Running create_histogram.py...\n"
python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/preseq/create_histogram.py \
    --input_path "${BAM_OUTPUT_PATH}/preprocessed_reads.csv" \
    --output_dir "${BAM_OUTPUT_PATH}/preseq"

echo -e "\nc_curve - all UMIs\n"
preseq c_curve -o "${BAM_OUTPUT_PATH}/preseq/c_curve_all.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_all.txt" -s 1000000 -v -r 42
python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/preseq/correct_reads.py \
  --preseq_file "${BAM_OUTPUT_PATH}/preseq/c_curve_all.txt" \
  --reads_col "total_reads" \
  --fractions_json "${BAM_OUTPUT_PATH}/preseq/fractions.json" \
  --key "included"

echo -e "\nlc_extrap - all UMIs\n"
preseq lc_extrap -o "${BAM_OUTPUT_PATH}/preseq/lc_extrap_all.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_all.txt" -e 5000000000 -s 1000000 -n 10000 -c 0.95 -v -r 42
python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/preseq/correct_reads.py \
  --preseq_file "${BAM_OUTPUT_PATH}/preseq/lc_extrap_all.txt" \
  --reads_col "TOTAL_READS" \
  --fractions_json "${BAM_OUTPUT_PATH}/preseq/fractions.json" \
  --key "included"

echo -e "\nbound_pop - all UMIs\n"
preseq bound_pop -o "${BAM_OUTPUT_PATH}/preseq/bound_pop_all.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_all.txt" -n 10000 -c 0.95 -v -r 42

echo -e "\nc_curve - tissue UMIs\n"
preseq c_curve -o "${BAM_OUTPUT_PATH}/preseq/c_curve_tissue.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_tissue.txt" -s 1000000 -v -r 42
python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/preseq/correct_reads.py \
  --preseq_file "${BAM_OUTPUT_PATH}/preseq/c_curve_tissue.txt" \
  --reads_col "total_reads" \
  --fractions_json "${BAM_OUTPUT_PATH}/preseq/fractions.json" \
  --key "included_tissue"

echo -e "\nlc_extrap - tissue UMIs\n"
preseq lc_extrap -o "${BAM_OUTPUT_PATH}/preseq/lc_extrap_tissue.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_tissue.txt" -e 5000000000 -s 1000000 -n 10000 -c 0.95 -v -r 42
python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/preseq/correct_reads.py \
  --preseq_file "${BAM_OUTPUT_PATH}/preseq/lc_extrap_tissue.txt" \
  --reads_col "TOTAL_READS" \
  --fractions_json "${BAM_OUTPUT_PATH}/preseq/fractions.json" \
  --key "included_tissue"

echo -e "\nbound_pop - tissue UMIs\n"
preseq bound_pop -o "${BAM_OUTPUT_PATH}/preseq/bound_pop_tissue.txt" -H "${BAM_OUTPUT_PATH}/preseq/histogram_tissue.txt" -n 10000 -c 0.95 -v -r 42