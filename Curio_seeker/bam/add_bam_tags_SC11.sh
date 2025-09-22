#!/bin/bash
#SBATCH --job-name=SC11_tags
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500G
#SBATCH --time=8:00:00

ml SAMtools
source ~/.bashrc
conda activate bamtools

PARQUET_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/curioseeker/work/8c/9a59f340d1856a5ae8376c55e9bf50/SC11-r1-r2-merged/"
BAM_IN="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/curioseeker/results/samtools/SC11.bam"
BAM_OUT="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output/SC11_tagged.bam"

python3 /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker/bam/add_bam_tags.py \
  --parquet_dir "$PARQUET_DIR" \
  --bam_in "$BAM_IN" \
  --bam_out "$BAM_OUT"