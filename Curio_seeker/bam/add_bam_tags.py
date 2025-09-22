#!/usr/bin/env python3
import os
import argparse
import pyarrow.parquet as pq
import pandas as pd
import pysam
import subprocess

# User inputs
parser = argparse.ArgumentParser(description="Add tags to Seeker BAM")
parser.add_argument("--parquet_dir", required=True, help="Directory containing parquet files")
parser.add_argument("--bam_in", required=True, help="Path to Seeker BAM file")
parser.add_argument("--bam_out", required=True, help="Path to output BAM file with tags")

args = parser.parse_args()
parquet_dir = args.parquet_dir
bam_in = args.bam_in
bam_out = args.bam_out

# Create output folder
os.makedirs(os.path.dirname(bam_out), exist_ok=True)

# Count total reads in input BAM
print("Input BAM:", bam_in, flush=True)
print("Counting total reads in input BAM...", flush=True)
total_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_in]))
print(f"Total reads in input BAM: {total_reads:,}", flush=True)

# Load parquet files
print("\nParquet directory:", parquet_dir, flush=True)
print("Loading parquet files...", flush=True)
tables = []
for file in os.listdir(parquet_dir):
    if file.endswith(".parquet"):
        path = os.path.join(parquet_dir, file)
        print(f"    Reading: {path}...", flush=True)
        tables.append(pq.read_table(path).to_pandas())

df = pd.concat(tables, ignore_index=True)
print(f"\nTotal number of rows: {len(df):,}", flush=True)

# Deduplicate reads by read1_id
print("\nDeduplicating reads...", flush=True)
dedup = df.drop_duplicates(subset=["read1_id"], keep="first").set_index("read1_id") # Only keep the first occurrence (all duplicated reads should be multimapping). This is somewhat risky.
print(f"Number of unique reads: {len(dedup):,}", flush=True)

# Process BAM and add tags
print("\nAnnotating BAM...", flush=True)

not_found = 0

with pysam.AlignmentFile(bam_in, "rb") as bam_infile, \
     pysam.AlignmentFile(bam_out, "wb", template=bam_infile) as bam_outfile:

    for i, read in enumerate(bam_infile):
        
        if i % 5_000_000 == 0 and i > 0:
            print(f"    Processed {i:,} reads...", flush=True)
            
        rid = read.query_name
        
        try:
            row = dedup.loc[rid]
        except KeyError:
            row = None
            not_found += 1

        if row is not None:
            if pd.notna(row["BM"]):
                read.set_tag("BM", row["BM"], value_type="Z")
            if pd.notna(row["cellular_barcode"]):
                read.set_tag("CB", row["cellular_barcode"], value_type="Z")
            if pd.notna(row["molecular_barcode"]):
                read.set_tag("UB", row["molecular_barcode"], value_type="Z")
            if pd.notna(row["read1_seq"]):
                read.set_tag("R1", row["read1_seq"], value_type="Z")
            if pd.notna(row["read2_seq"]):
                read.set_tag("R2", row["read2_seq"], value_type="Z")
            if pd.notna(row["r1_proper_structure_matched"]):
                read.set_tag("PS", int(row["r1_proper_structure_matched"]), value_type="i")
        bam_outfile.write(read)
        
print(f"\nBAM written to {bam_out}", flush=True)
print(f"Number of reads from input BAM not found in parquet files: {not_found:,}", flush=True)
        
print("\nCounting total reads in output BAM...", flush=True)
total_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_out]))
print(f"Total reads in output BAM: {total_reads:,}", flush=True)
        
### Colum descriptions according to spatial_bioinformatics_support@takarabio.com
# read1_id - fastq read id
# read1_seq - fastq sequence of R1
# cellular_barcode - cell barcode sequence
# molecular_barcode - molecular barcode sequence
# illumina_primer - is illumina primer detected
# proper_structure - R1 contains proper sequence structure
# BB - child cell barcode sequence
# BM - parent cell barcode sequence
# BD - number nucleotide differences BB to BM
# molecular_barcode_count_distinct - unique count of molecular barcodes
# r1_proper_structure_matched - R1 proper structure + cell bc + molecular bc
# read2_seq - fastq sequence of R2
# reference_id - id of reference
# reference_start - start of reference
# XS - assigned reference tag
# XT - assigned reference gene symbol