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
        output_path: str| Path,
):
    if os.path.exists( output_path ):
        raise FileExistsError( f"File {output_path} exists." )
    os.makedirs( os.path.dirname( output_path ), exist_ok=True )

    # Load barcodes in tissue
    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)

    start = time.time()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with open( output_path, "w" ) as f_out:
            f_out.write("index,included,included_tissue,barcode,umi,gene,random\n")
            # i is the total number of reads
            for i, read in enumerate(bam):
                r = random.random()
                step = 1000000
                if i%step == 0:
                    print( f"{i} reads processed.", flush=True )
                    print(time.time() - start, flush=True)

                barcode = 0
                umi = 0
                gene = 0
                included = 0
                included_tissue = 0

                mapq = read.mapping_quality
                if mapq == 255 and read.has_tag("CB") and read.has_tag("UB") and read.has_tag("GN"):

                    barcode = read.get_tag("CB")
                    umi = read.get_tag("UB")
                    gene = read.get_tag("GN")

                if barcode and umi and gene:
                    included = 1
                    if barcode in tissue_barcodes:
                        included_tissue=1

                f_out.write(f"{i},{included},{included_tissue},{barcode},{umi},{gene},{r}\n")
                
    print( f"Finished, took {time.time()-start}", flush=True )
                

if __name__ == "__main__":

    # User inputs
    parser = argparse.ArgumentParser(description="BAM preprocess")
    parser.add_argument("--bam_path", required=True, help="Path to BAM file")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--output_path", required=True, help="Output path. Preprocessed csv file will be stored there.")

    args = parser.parse_args()

    main(
        bam_path=args.bam_path,
        tissue_barcodes_path=args.tissue_barcodes_path,
        output_path = args.output_path,
        )


