#!/bin/bash
#SBATCH --job-name=BMK_STAR
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=10:00:00

# Commands to execute
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star/ --genomeFastaFiles fasta/genome.fa --sjdbGTFfile genes/genes.gtf
