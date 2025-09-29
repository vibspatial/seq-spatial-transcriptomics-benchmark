import os
import time
import argparse
from pathlib import Path
import pandas as pd
import random

random.seed(42)

def main(
        umi_gene_path: str | Path,
        features_path: str | Path,
        tissue_barcodes_path: str | Path,
        output_path: str | Path,
        fastq_read_count: int,
):
    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_path} exists.")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print("Reading umi_gene.tsv.gz...", flush=True)
    df = pd.read_csv(
        umi_gene_path,
        sep="\t",
        names=["barcode", "umi", "gene_id", "read_count"],
        compression="gzip"
    )

    print("Reading features.tsv...", flush=True)
    features = pd.read_csv(
        features_path,
        sep="\t",
        names=["gene_id", "gene", "type"],
    )
    gene_lookup = dict(zip(features["gene_id"], features["gene"]))

    print(f"Collapsed rows: {len(df):,}", flush=True)

    # Replace gene_id with gene name
    df["gene"] = df["gene_id"].map(gene_lookup)
    df.drop(columns=["gene_id"], inplace=True)

    # Expand rows by read_count
    print("Exploding rows by read_count...", flush=True)
    df_exploded = df.loc[df.index.repeat(df["read_count"])].reset_index(drop=True)

    print(f"Rows after explosion: {len(df_exploded):,}", flush=True)

    # Load tissue barcodes
    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)

    # Flags
    df_exploded["included"] = (
        df_exploded["barcode"].notna() &
        df_exploded["umi"].notna() &
        df_exploded["gene"].notna()
    ).astype(int)

    df_exploded["included_tissue"] = (
        (df_exploded["included"] == 1) &
        df_exploded["barcode"].isin(tissue_barcodes)
    ).astype(int)

    df_exploded["random"] = [random.random() for _ in range(len(df_exploded))]

    # Pad missing reads
    n_missing = fastq_read_count - len(df_exploded)
    if n_missing > 0:
        print(f"Adding {n_missing:,} placeholder rows for missing reads...", flush=True)
        filler = pd.DataFrame({
            "barcode": [None] * n_missing,
            "umi": [None] * n_missing,
            "gene": [None] * n_missing,
            "included": [0] * n_missing,
            "included_tissue": [0] * n_missing,
            "random": [random.random() for _ in range(n_missing)]
        })
        df_exploded = pd.concat([df_exploded, filler], ignore_index=True)

    print(f"Final row count: {len(df_exploded):,} (should match FASTQ)", flush=True)

    # Write output
    print(f"\nWriting CSV to {output_path}", flush=True)
    start = time.time()
    with open(output_path, "w") as f_out:
        f_out.write("index,included,included_tissue,barcode,umi,gene,random\n")
        for i, row in df_exploded.iterrows():
            f_out.write(f"{i},{row.included},{row.included_tissue},{row.barcode},{row.umi},{row.gene},{row.random}\n")
            if i % 5_000_000 == 0 and i > 0:
                print(f"{i:,} rows written: {round((time.time()-start)/60,2)} min", flush=True)

    print("Done!", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Explode umi_gene.tsv.gz to per-read format")
    parser.add_argument("--umi_gene_path", required=True, help="Path to umi_gene.tsv.gz")
    parser.add_argument("--features_path", required=True, help="Path to features.tsv")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--output_path", required=True, help="Output CSV path")
    parser.add_argument("--fastq_read_count", type=int, required=True, help="Total number of reads from FASTQ")

    args = parser.parse_args()

    main(
        umi_gene_path=args.umi_gene_path,
        features_path=args.features_path,
        tissue_barcodes_path=args.tissue_barcodes_path,
        output_path=args.output_path,
        fastq_read_count=args.fastq_read_count,
    )