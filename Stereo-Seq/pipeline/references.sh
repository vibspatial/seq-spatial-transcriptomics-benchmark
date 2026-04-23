## Mouse
saw makeRef \
    --mode=STAR \
    --fasta=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A/fasta/genome.fa \
    --gtf=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-mm10-2020-A/genes/genes.gtf \
    --genome=./refdata-gex-mm10-2020-A

## Human
saw makeRef \
    --mode=STAR \
    --fasta=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    --gtf=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    --genome=./refdata-gex-GRCh38-2020-A

## Rice
# Identify and remove duplicate loci
awk '$3=="gene"{print $1,$4,$5,$7}' /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/IRGSP.49/genes/genes.gtf | sort | uniq -d > rice/duplicated_coords.txt

awk 'NR==FNR{dup[$1":"$2":"$3":"$4]; next} {
  key=$1":"$4":"$5":"$7;
  if(!(key in dup)) print
}' rice/duplicated_coords.txt /data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/IRGSP.49/genes/genes.gtf > rice/genes.no_dups.gtf

# Create reference for SAW
saw makeRef \
    --mode=STAR \
    --fasta=/data/groups/technologies/spatial.catalyst/Julien/spaceranger/references/IRGSP.49/fasta/genome.fa \
    --gtf=./rice/genes.no_dups.gtf \
    --genome=./rice/IRGSP.49 \
    --params-config='--genomeSAindexNbases=13'
