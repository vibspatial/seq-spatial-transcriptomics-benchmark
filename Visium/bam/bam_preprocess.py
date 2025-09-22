import os
import time
import argparse
from pathlib import Path

import pandas as pd
import pysam

import random
random.seed(42)

def main(
        bam_path : str | Path,
        tissue_barcodes_path: str | Path,
        probe_set_path: str| Path,
        output_path: str| Path,
        size: int = 0,
):
    if os.path.exists( output_path ):
        raise FileExistsError( f"File {output_path} exists." )
    os.makedirs( os.path.dirname( output_path ), exist_ok=True )

    # Load barcodes in tissue
    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)

    # Get included genes from probe set
    probe_set_df = pd.read_csv(probe_set_path, comment='#')
    included_genes = set(probe_set_df[probe_set_df['included']]['gene_id'])

    start = time.time()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with open( output_path, "w" ) as f_out:
            f_out.write( "index,included,included_tissue,barcode,umi,gene,"
                        "q30_read_bases,q30_read_total,q30_barcode_bases,q30_barcode_total,q30_umi_bases,q30_umi_total,"
                        "reads_mapped,reads_confident,reads_confident_filtered,reads_half_mapped,reads_split_mapped,random\n" )
        # i is the total number of reads
            for i, read in enumerate(bam):
                r = random.random()
                step = 1000000
                if i%step == 0:
                    print( f"{i} reads processed.", flush=True )
                    print(time.time() - start, flush=True)
                if size:
                    if i > size:
                        break

                barcode = 0
                umi = 0
                gene = 0
                included = 0
                included_tissue = 0

                mapq = read.mapping_quality
                if mapq == 255 and read.has_tag("CB") and read.has_tag("UB") and read.has_tag("GN"): #  28 s full, scales lineary, so ok

                    if read.has_tag("GX") and read.get_tag("GX") in included_genes:
                        barcode = read.get_tag("CB")
                        umi = read.get_tag("UB")
                        gene = read.get_tag("GN")


                if barcode and umi and gene:
                    included = 1 # for easier subsampling
                    if barcode in tissue_barcodes:
                        included_tissue=1

                # if read does not have query_qualities, set it to zero for this read
                q30_read_bases = 0
                q30_read_total = 0
                if read.query_qualities:
                    q30_read_bases = sum(q >= 30 for q in read.query_qualities if q > 2)
                    q30_read_total = sum(q > 2 for q in read.query_qualities)

                # if read does not has tag "CY", set it to zero
                q30_barcode_bases=0
                q30_barcode_total=0
                if read.has_tag("CY"):
                    q_barcode_scores = [ord(c) - 33 for c in read.get_tag("CY")]
                    q30_barcode_bases = sum(q >= 30 for q in q_barcode_scores if q > 2)
                    q30_barcode_total = sum(q > 2 for q in q_barcode_scores)

                q30_umi_bases=0
                q30_umi_total =0
                if read.has_tag("UY"):
                    q_umi_scores = [ord(c) - 33 for c in read.get_tag("UY")]
                    q30_umi_bases = sum(q >= 30 for q in q_umi_scores if q > 2)
                    q30_umi_total = sum(q > 2 for q in q_umi_scores)
                
                reads_mapped = 0
                if mapq >0:
                    reads_mapped =1

                reads_confident = 0
                reads_confident_filtered = 0
                reads_half_mapped = 0
                reads_split_mapped = 0
                if mapq == 255:
                    reads_confident= 1
                    if read.has_tag("GX") and read.get_tag("GX") in included_genes:
                        reads_confident_filtered = 1
                if mapq == 1:
                    reads_half_mapped = 1
                if mapq == 3:
                    reads_split_mapped = 1

                f_out.write( f"{i},{included},{included_tissue},{barcode},{umi},{gene},"
                            f"{q30_read_bases},{q30_read_total},{q30_barcode_bases},{q30_barcode_total},{q30_umi_bases},{q30_umi_total},"
                            f"{reads_mapped},{reads_confident},{reads_confident_filtered},{reads_half_mapped},{reads_split_mapped},{r}\n" )
                
    print( f"Finished, took {time.time()-start}", flush=True )
                

if __name__ == "__main__":

    # User inputs
    parser = argparse.ArgumentParser(description="BAM preprocess")
    parser.add_argument("--bam_path", required=True, help="Path to BAM file")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--probe_set_path", required=True, help="Path to probe_set.csv")
    parser.add_argument("--size", type=int, required=False, default=0, help="If specified, will only preprocess the first 'size' reads.")
    parser.add_argument("--output_path", required=True, help="Output path. Preprocessed csv file will be stored there.")

    args = parser.parse_args()

    main(
        bam_path=args.bam_path,
        tissue_barcodes_path=args.tissue_barcodes_path,
        probe_set_path=args.probe_set_path,
        size=args.size,
        output_path = args.output_path,
        )


