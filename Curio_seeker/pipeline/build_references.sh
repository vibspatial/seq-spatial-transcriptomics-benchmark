#!/bin/bash
#SBATCH --job-name=curioseeker-ref
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00

# Environment: star_2.6.1d

bash /home/VIB.LOCAL/julien.mortier/curioseeker/generate_seeker_reference/generate_seeker_reference.sh \
  /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A/fasta/genome.fa \
  /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A/genes/genes.gtf \
  chrM \
  refdata-gex-mm10-2020-A

bash /home/VIB.LOCAL/julien.mortier/curioseeker/generate_seeker_reference/generate_seeker_reference.sh \
  /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  chrM \
  refdata-gex-GRCh38-2020-A