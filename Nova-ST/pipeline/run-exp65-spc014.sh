#!/bin/bash
#SBATCH --account=lp_edu_ont2024
#SBATCH --clusters=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7:00:00           # hh:mm:ss
#SBATCH --mem=200G
#SBATCH --job-name=novascope65
#SBATCH --mail-user=frank.vernaillen@vib.be
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT

export LMOD_CACHED_LOADS=yes
export LMOD_RC=/data/leuven/software/biomed/${VSC_INSTITUTE_CLUSTER}_${VSC_ARCH_LOCAL}_${VSC_OS_LOCAL}/lmodrc.lua

echo LMOD_RC=$LMOD_RC
echo LMOD_CACHED_LOADS=$LMOD_CACHED_LOADS

module use /data/leuven/software/biomed/${VSC_INSTITUTE_CLUSTER}_${VSC_ARCH_LOCAL}_${VSC_OS_LOCAL}/20*/modules/all

module load cluster/wice/batch_icelake
module load GCC/12.3.0 GDAL/3.7.1-foss-2023a_SciPy_2311 ImageMagick/7.1.1-15-GCCcore-12.3.0 spatula/1.0.0-GCC-12.3.0 SAMtools/1.21-GCC-12.3.0 STAR/2.7.11a-GCC-12.3.0

smk_dir=$VSC_DATA_VO/spatial/nova-st/git/NovaScope
job_dir=$VSC_DATA_VO/spatial/nova-st/data/2024-07-UTBenchmark-SpC/Exp65-sampleSPC014

source /data/gent/426/vsc42627/miniconda3/bin/activate snakemake
snakemake --latency-wait 20 -s ${smk_dir}/NovaScope.smk -d $job_dir --cores 32 --rerun-incomplete
