import os
import argparse
import dask.dataframe as dd
from pathlib import Path
import json
import uuid
import shutil

def main(
    input_path: str | Path,
    output_dir: str | Path,
):
    os.makedirs(output_dir, exist_ok=True)

    print(f"Reading preprocessed CSV: {input_path}", flush=True)
    df_csv = dd.read_csv(
        input_path,
        usecols=[
            "included",
            "included_tissue",
            "barcode",
            "umi",
            "gene",
            "random",
        ],
        dtype={
            "included": "bool",
            "included_tissue": "bool",
            "barcode": "string",
            "umi": "string",
            "gene": "string",
            "random": "float32",
        },
    )
        
    print("Writing full dataset to parquet file...", flush=True) 
    _parquet_path = os.path.join(output_dir, f"full_{uuid.uuid4()}.parquet") 
    df_csv.to_parquet(_parquet_path) 
    print("Done writing parquet file", flush=True)
    
    df = dd.read_parquet(_parquet_path)
    
    total_reads = df.shape[0].compute()
    print("Number of reads:", total_reads, flush=True)
    
    # Calculate fractions of reads that contribute to UMI counts
    def _compute_fraction(df: dd.DataFrame, col: str) -> float:
        return (df[col] == 1).sum().compute() / total_reads

    fractions = {
        "included": _compute_fraction(df, "included"),
        "included_tissue": _compute_fraction(df, "included_tissue"),
    }
    
    fractions_path = os.path.join(output_dir, "fractions.json")
    with open(fractions_path, "w") as f:
        json.dump(fractions, f, indent=2)
    print(f"UMI read fractions saved to {fractions_path}", flush=True)
    
    # Create histograms
    def _create_histogram(
        df: dd.DataFrame,
        included_col: str,
        output_file: str | Path,
    ):
        # Keep only reads that contribute to UMI counts
        df = df[df[included_col] == 1]

        # Group by barcode/UMI/gene to count number of reads per UMI
        print("Grouping by barcode/UMI/gene to count reads per UMI...", flush=True)
        umi_counts = df.groupby(["barcode", "umi", "gene"]).size()

        # Build histogram
        print("Building histogram for preseq...", flush=True)
        histogram = (
            umi_counts.map_partitions(lambda s: s.value_counts())
            .compute()
            .groupby(level=0) 
            .sum()
            .sort_index()
        )
        histogram.to_csv(output_file, sep="\t", header=False)

        print(f"Histogram saved to {output_file}", flush=True)
    
    _create_histogram(df, "included", os.path.join(output_dir, "histogram_all.txt"))
    _create_histogram(df, "included_tissue", os.path.join(output_dir, "histogram_tissue.txt"))    
    
    # Clean up parquet 
    if os.path.exists(_parquet_path): 
        shutil.rmtree(_parquet_path)

    print("Done!", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert preprocessed CSV to preseq histogram")
    parser.add_argument("--input_path", required=True, help="Path to preprocessed CSV file")
    parser.add_argument("--output_dir", required=True, help="Directory where output will be written")
    args = parser.parse_args()
    main(args.input_path, args.output_dir)
