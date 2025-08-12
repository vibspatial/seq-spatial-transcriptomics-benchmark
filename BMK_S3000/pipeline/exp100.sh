#!/bin/bash
#SBATCH --job-name=BMK_S3000
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500G
#SBATCH --time=120:00:00

# Commands to execute
module load pigz
module load Perl/5.38.2-GCCcore-13.3.0
module load SAMtools

SEQ_DIR=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/BMK_S3000/sequencing
FILE_NAME_1=1_4941_02_Exp100_SC57_SPC_004_S6
FILE_NAME_2=2_4941_02_Exp100_SC57_SPC_004_S5
FILE_NAME_3=3_4941_02_Exp100_SC57_SPC_004_S5
FILE_NAME=4941_02_Exp100_SC57_SPC_004
EXP=exp100
CONFIG=/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S3000/config_${EXP}.txt

# pigz -dc -p 24 ${SEQ_DIR}/${FILE_NAME_1}_L001_R1_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_1}_L002_R1_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_2}_L001_R1_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_2}_L002_R1_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_3}_L001_R1_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_3}_L002_R1_001.fastq.gz | pigz -p 24 > ${SEQ_DIR}/${FILE_NAME}_R1.fastq.gz
# pigz -dc -p 24 ${SEQ_DIR}/${FILE_NAME_1}_L001_R2_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_1}_L002_R2_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_2}_L001_R2_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_2}_L002_R2_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_3}_L001_R2_001.fastq.gz ${SEQ_DIR}/${FILE_NAME_3}_L002_R2_001.fastq.gz | pigz -p 24 > ${SEQ_DIR}/${FILE_NAME}_R2.fastq.gz

/home/VIB.LOCAL/julien.mortier/BMK_software/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix_v2.4.f.2_English/BSTMatrix -c ${CONFIG} -s 2,3,4,5,6,8
