import os
import time
import argparse
from pathlib import Path
import random
import pysam

random.seed(42)

def main(
        bam_path: str | Path,
        tissue_barcodes_path: str | Path,
        output_path: str | Path,
        fastq_read_count: int,
):
    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_path} exists.")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Load tissue barcodes
    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)

    # Create output CSV
    start = time.time()
    reads_included = 0
    
    with open(output_path, "w") as f_out:
        f_out.write("index,included,included_tissue,barcode,umi,gene,random\n")
        
        # Get included reads from BAM
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                r = random.random()
                barcode = 0
                umi = 0
                gene = 0
                included = 0
                included_tissue = 0

                if read.has_tag("CB") and read.has_tag("UB") and read.has_tag("GN"):
                    if read.has_tag("sM") and read.get_tag("sM") >= 0 and read.has_tag("NH") and read.get_tag("NH") == 1:
                        barcode = read.get_tag("CB")
                        umi = read.get_tag("UB")
                        gene = read.get_tag("GN")
                        
                if barcode and umi and gene:
                    included = 1
                    reads_included += 1
                    if barcode in tissue_barcodes:
                        included_tissue = 1

                    f_out.write(f"{reads_included},{included},{included_tissue},{barcode},{umi},{gene},{r}\n")
                    if reads_included % 5_000_000 == 0 and reads_included > 0:
                        print(f"{reads_included:,} rows written: {round((time.time()-start)/60,2)} min", flush=True)

        # Pad missing reads
        n_missing = fastq_read_count - reads_included
        if n_missing > 0:
            print(f"\n {reads_included}/{fastq_read_count} already included in output CSV.")
            print(f"Adding {n_missing:,} placeholder rows for missing reads...", flush=True)
            for i in range(n_missing):
                reads_included += 1
                f_out.write(f"{reads_included},0,0,None,None,None,{random.random()}\n")
                if reads_included % 5_000_000 == 0 and reads_included > 0:
                    print(f"{reads_included:,} rows written: {round((time.time()-start)/60,2)} min", flush=True)

    print(f"\nDone! {reads_included:,} rows written: {round((time.time()-start)/60,2)} min", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BAM preprocess")
    parser.add_argument("--bam_path", required=True, help="Path to BAM file")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--output_path", required=True, help="Output CSV path")
    parser.add_argument("--fastq_read_count", type=int, required=True, help="Total number of reads from FASTQ")

    args = parser.parse_args()

    main(
        bam_path=args.bam_path,
        tissue_barcodes_path=args.tissue_barcodes_path,
        output_path=args.output_path,
        fastq_read_count=args.fastq_read_count,
    )


