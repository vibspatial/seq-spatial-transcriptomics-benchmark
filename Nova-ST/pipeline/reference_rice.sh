# Instructions to built a reference: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr

# We used the rice reference that was used for Space Ranger by the Single-Cell Core.

## index the reference genome
cd /staging/leuven/stg_00143/julienm/nova-st/ref/oryza_sativa_japonica/IRGSP.49
gunzip -c ./genes/genes.gtf.gz > ./genes/genes.gtf
rm ./genes/genes.gtf.gz

ml STAR/2.7.10a-GCC-10.3.
STAR --runMode genomeGenerate \
    --runThreadN 1 \
    --genomeDir ./indexes/STAR/2.7.10a \
    --genomeFastaFiles ./fasta/genome.fa \
    --genomeSAindexNbases 13 \
    --genomeChrBinNbits 18 \
    --genomeSAsparseD 3 \
    --limitGenomeGenerateRAM 17179869184 \
    --sjdbGTFfile ./genes/genes.gtf

# Copy genome + GTF into STAR index
cp fasta/genome.fa indexes/STAR/2.7.10a/
cp genes/genes.gtf indexes/STAR/2.7.10a/

# Create genelists 
# I am not sure whether this is entirely correct. According to the NovaScope documentation, these are just to create some PNG plots
# and these do not work in the docker version anyway. I should be able to turn of these visualizations by setting action=False for drawsge, but
# this still results in an path not found error for the genelist so I quickly created these files.

mkdir genelists

awk '$3 == "gene" {
  match($0, /gene_id "([^"]+)"/, arr);
  if (arr[1] != "") print arr[1];

}' genes/genes.gtf > genelists/all.genes.tsv

awk '$3 == "gene" && $1 == "Mt" {
  match($0, /gene_id "([^"]+)"/, arr);
  if (arr[1] != "") print arr[1];
}' genes/genes.gtf > genelists/MT.genes.tsv

grep -v -f genelists/MT.genes.tsv genelists/all.genes.tsv > genelists/nonMT.genes.tsv

awk '$3 == "gene" && $0 ~ /gene_biotype "rRNA"/ {
  match($0, /gene_id "([^"]+)"/, arr);
  if (arr[1] != "") print arr[1];
}' genes/genes.gtf > genelists/ribosomal.genes.tsv

