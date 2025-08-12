#!/bin/bash
#SBATCH --job-name=curioseeker-merge
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=4:00:00

SEQ_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Curio_seeker/sequencing
cd "$SEQ_DIR"

# SC11
pigz -dc -p 24 *_SC11_*_L00{1,2}_R1_001.fastq.gz | pigz -p 24 > SC11_R1.fastq.gz
pigz -dc -p 24 *_SC11_*_L00{1,2}_R2_001.fastq.gz | pigz -p 24 > SC11_R2.fastq.gz

# SC12
pigz -dc -p 24 *_SC12_*_L00{1,2}_R1_001.fastq.gz | pigz -p 24 > SC12_R1.fastq.gz
pigz -dc -p 24 *_SC12_*_L00{1,2}_R2_001.fastq.gz | pigz -p 24 > SC12_R2.fastq.gz

# SC13
pigz -dc -p 24 *_SC13_*_L00{1,2}_R1_001.fastq.gz | pigz -p 24 > SC13_R1.fastq.gz
pigz -dc -p 24 *_SC13_*_L00{1,2}_R2_001.fastq.gz | pigz -p 24 > SC13_R2.fastq.gz

# SC14
pigz -dc -p 24 *_SC14_*_L00{1,2}_R1_001.fastq.gz | pigz -p 24 > SC14_R1.fastq.gz
pigz -dc -p 24 *_SC14_*_L00{1,2}_R2_001.fastq.gz | pigz -p 24 > SC14_R2.fastq.gz

echo "Merging complete."