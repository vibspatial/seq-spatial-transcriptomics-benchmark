import os
import uuid
import time
import shutil

import dask.dataframe as dd
from dask.dataframe import DataFrame
import pandas as pd
from pathlib import Path
import argparse

def calculate_metrics_all(
        input_path: str | Path,
        output_dir: str | Path,
        all_barcodes_path: str | Path,
        tissue_area: float = 1,
        step_size: int = 10000000,
) -> DataFrame:
    os.makedirs(output_dir, exist_ok=True)

    df = dd.read_csv(
        input_path,
        dtype={
            "included": "bool",
            "included_tissue": "bool",
            "barcode": "string",
            "umi": "string",
            "gene": "string",
            "random": "float32",
        },
    )

    # Back it to parquet for faster processing
    _parquet_path = os.path.join(output_dir, f"full_{uuid.uuid4()}.parquet")
    print("Writing to parquet file...", flush=True)
    df.to_parquet(_parquet_path)
    df = dd.read_parquet(_parquet_path)
    print("Done writing parquet file", flush=True)

    total_size = len(df)
    print("Total size:", total_size, flush=True)

    sample_sizes = []
    for i in range(1, total_size // step_size + 1):
        sample_sizes.append(i * step_size)

    if total_size not in sample_sizes:
        sample_sizes.append(total_size)
    print("All sample sizes:", sample_sizes, flush=True)

    for _sample_size in sample_sizes:
        print(f"Sample size {_sample_size:,}", flush=True)
        _output_dir_sample_size = os.path.join(output_dir, f"sample_size_{_sample_size}")
        os.makedirs(_output_dir_sample_size)
        if _sample_size == total_size:
            _sample_size = 0  # do not sample if equal to total size
        calculate_metrics(
            df=df,
            output_dir=_output_dir_sample_size,
            all_barcodes_path=all_barcodes_path,
            tissue_area=tissue_area,
            sample_size=_sample_size,
        )

    # Clean up parquet
    if os.path.exists(_parquet_path):
        shutil.rmtree(_parquet_path)

    return df


def calculate_metrics(
        df: DataFrame,
        output_dir: str | Path,
        all_barcodes_path: str | Path,
        tissue_area: float = 1,
        sample_size: int = 0,  # with sample_size 0, we do not sample
):
    t0 = time.time()

    with open(all_barcodes_path) as f:
        all_barcodes = set(line.strip() for line in f)

    # sample
    if sample_size:
        start = time.time()
        frac = sample_size / len(df)
        print(f"\nSampling {sample_size:,} reads...", flush=True)
        df = df[df["random"] <= frac]
        _parquet_path = os.path.join(output_dir, f"sample_{uuid.uuid4()}.parquet")
        df.to_parquet(_parquet_path)
        df = dd.read_parquet(_parquet_path)
        print(f"Sampling done - {round(time.time() - start)} seconds.", flush=True)

    # compute metrics
    df_included = df[df["included"]]
    print("Grouping UMIs...", flush=True)
    start = time.time()
    df_included_grouped = df_included.groupby(["barcode", "umi", "gene"]).size().compute()
    print(f"Grouping UMIs done - {round(time.time() - start)} seconds", flush=True)

    df_included_tissue = df[df["included_tissue"]]
    print("Grouping UMIs in tissue...", flush=True)
    start = time.time()
    df_included_tissue_grouped = df_included_tissue.groupby(["barcode", "umi", "gene"]).size().compute()
    print(f"Grouping UMIs in tissue done - {round(time.time() - start)} seconds", flush=True)

    print("Calculating capture efficiency metrics...", flush=True)
    start = time.time()
    total_umis = df_included_grouped.sum()
    unique_umis = len(df_included_grouped)

    tissue_unique_umis = len(df_included_tissue_grouped)
    gene_counts = df_included_grouped.groupby("gene").sum()
    unique_gene_counts = df_included_grouped.groupby("gene").count()

    # Same, but in tissue
    gene_counts_tissue = df_included_tissue_grouped.groupby("gene").sum()
    unique_gene_counts_tissue = df_included_tissue_grouped.groupby("gene").count()

    gene_counts = gene_counts.to_frame(name="gene_counts").reset_index()
    unique_gene_counts = unique_gene_counts.to_frame(name="unique_gene_counts").reset_index()
    gene_counts_tissue = gene_counts_tissue.to_frame(name="gene_counts_tissue").reset_index()
    unique_gene_counts_tissue = unique_gene_counts_tissue.to_frame(name="unique_gene_counts_tissue").reset_index()

    gene_counts_merged = (
        gene_counts
        .merge(unique_gene_counts, on="gene", how="outer")
        .merge(gene_counts_tissue, on="gene", how="outer")
        .merge(unique_gene_counts_tissue, on="gene", how="outer")
    )

    assert tissue_area > 0, "Please provide tissue_area > 0."

    for _name in ["gene_counts", "unique_gene_counts", "gene_counts_tissue", "unique_gene_counts_tissue"]:
        gene_counts_merged[f"{_name}_corrected"] = gene_counts_merged[_name] / tissue_area

    gene_counts_merged.to_csv(os.path.join(output_dir, "gene_counts.csv"))

    median_gene_counts = gene_counts_merged["gene_counts"].median()
    median_unique_gene_counts = gene_counts_merged["unique_gene_counts"].median()
    median_gene_counts_tissue = gene_counts_merged["gene_counts_tissue"].median()
    median_unique_gene_counts_tissue = gene_counts_merged["unique_gene_counts_tissue"].median()

    median_reads_per_unique_umi = df_included_grouped.median()
    median_reads_per_unique_umi_tissue = df_included_tissue_grouped.median()
    print(f"Calculating capture efficiency metrics done - {round(time.time() - start)} seconds", flush=True)
    def safe_div(a, b): return a / b if b > 0 else 0

    total_reads = len(df)
    assert total_reads > 0, "Total number of reads is zero."

    # put it in dictionary, then write to csv
    row = {
        "number_of_spots_with_at_least_1_UMI_in_tissue": df_included_tissue_grouped.index.get_level_values("barcode").nunique(),
        "percentage_spots_recovered":df_included_tissue_grouped.index.get_level_values("barcode").nunique()/len(all_barcodes)*100,
        "total_number_of_reads": total_reads,
        "total_number_of_genes": len(gene_counts),
        "total_number_of_UMIs": total_umis,
        "total_number_of_unique_UMIs": unique_umis,
        "number_of_genes_in_tissue": len(gene_counts_tissue),
        "number_of_unique_UMIs_in_tissue": tissue_unique_umis,
        "percentage_of_genes_in_tissue": safe_div(len(gene_counts_tissue) ,len(gene_counts))*100,
        "percentage_of_unique_UMIs_in_tissue": safe_div(tissue_unique_umis, unique_umis)*100,
        "sequencing_saturation_reads": (1-safe_div(unique_umis, total_reads))*100 if total_reads>0 else 0,
        "sequencing_saturation_UMIs": (1-safe_div(unique_umis, total_umis))*100,
        "unique_UMIs_per_mm2_in_tissue": safe_div(tissue_unique_umis, tissue_area),
        "median_UMIs_per_gene": median_gene_counts,
        "median_UMIs_per_gene_per_mm2": safe_div(median_gene_counts,tissue_area),
        "median_unique_UMIs_per_gene": median_unique_gene_counts,
        "median_unique_UMIs_per_gene_per_mm2": safe_div(median_unique_gene_counts, tissue_area),
        "median_unique_UMIs_per_gene_in_tissue": median_unique_gene_counts_tissue,
        "median_unique_UMIs_per_gene_in_tissue_per_mm2": safe_div(median_unique_gene_counts_tissue, tissue_area),
        "percentage_of_reads_with_usable_UMI": safe_div(total_umis, total_reads)*100,
        "ratio_reads_per_unique_UMI": safe_div(total_reads, unique_umis),
        "median_reads_per_unique_UMI": median_reads_per_unique_umi,
        "median_reads_per_unique_UMI_tissue": median_reads_per_unique_umi_tissue,
    }
    
    for key, value in row.items():
        print(f"{key}: {value:,}", flush=True)

    # clean up parquet
    if sample_size>0 and os.path.exists(_parquet_path):
        shutil.rmtree(_parquet_path)
    pd.DataFrame([row]).to_csv(
        os.path.join(output_dir, "metrics.csv"),
        index=False,
    )

    t1 = time.time()
    print(f"Sample size {sample_size:,} finished! - {round((t1 - t0) / 60, 2)} minutes", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="UMI/gene metrics calculation")
    parser.add_argument("--input_path", required=True, help="Path to per-read CSV created by preprocess.py")
    parser.add_argument("--all_barcodes_path", required=True, help="Path to all_barcodes.txt")
    parser.add_argument("--output_dir", required=True, help="Output directory where metrics and intermediate files will be stored.")
    parser.add_argument("--tissue_area", type=float, required=False, default=1, help="Tissue area in mm2.")
    parser.add_argument("--step_size", type=int, required=False, default=10000000, help="Step size.")
    args = parser.parse_args()

    calculate_metrics_all(
        input_path=args.input_path,
        all_barcodes_path=args.all_barcodes_path,
        output_dir=args.output_dir,
        tissue_area=args.tissue_area,
        step_size=args.step_size,
    )
