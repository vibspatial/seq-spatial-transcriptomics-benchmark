#!/bin/bash
#SBATCH --job-name=exp66_downsampling
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=128G

# Commands to execute
python exp66_downsampling.py
