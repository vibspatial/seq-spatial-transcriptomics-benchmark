import pysam
import pandas as pd
import numpy as np
import random
import os
from collections import Counter
import argparse
import subprocess
import gc
import glob
import pyarrow.parquet as pq
import pyarrow as pa

# User inputs
parser = argparse.ArgumentParser(description="BAM metrics calculation")
parser.add_argument("--bam_path", required=True, help="Path to BAM file")
parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
parser.add_argument("--output_path", required=True, help="Output folder for metrics")
parser.add_argument("--tissue_area_mm2", type=float, required=True, help="Tissue area in mm^2")
parser.add_argument("--step_size", type=int, default=10_000_000, help="Step size for subsampling reads")
parser.add_argument("--testing", action='store_true', help="Enable testing mode with smaller BAM")

args = parser.parse_args()

bam_path = args.bam_path
tissue_barcodes_path = args.tissue_barcodes_path
output_path = args.output_path
tissue_area_mm2 = args.tissue_area_mm2
step_size = args.step_size
testing = args.testing

random.seed(42)

# Create output folders
os.makedirs(output_path, exist_ok=True)
os.makedirs(os.path.join(output_path, "gene_counts_per_subsample"), exist_ok=True)
os.makedirs(os.path.join(output_path, "batches_tmp"), exist_ok=True)

# Subsample BAM to a smaller size for testing
if testing:
    bam_out = f"{output_path}/bam_small.bam"
    num_reads = 1_000_000
    step_size = 200_000
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pysam.AlignmentFile(bam_out, "wb", template=bam) as out:
            for i, read in enumerate(bam):
                if i >= num_reads:
                    break
                out.write(read)
    print(f"Written {num_reads} reads to {bam_out}", flush=True)
    
    bam_path = bam_out

# Load barcodes in tissue
with open(tissue_barcodes_path) as f:
    tissue_barcodes = set(line.strip() for line in f)
    
# Count total reads
print("Counting total reads...", flush=True)
total_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_path]))
print(f"Total reads: {total_reads:,}", flush=True)

# Define subsample levels
subsample_levels = list(range(step_size, total_reads, step_size))
if subsample_levels[-1] != total_reads:
    subsample_levels.append(total_reads)
fractions = {n: n / total_reads for n in subsample_levels}

# Initialize metrics
all_umis = []
batch_index = 0

metrics = {
    n: {
        "total_reads": 0,
        "q30_read_2_bases": 0,
        "q30_read_2_total": 0,
        "reads_Assigned": 0,
        "reads_Unassigned_NoFeatures": 0,
        "reads_Unassigned_MultiMapping": 0,
        "reads_Unassigned_Ambiguity": 0,
    }
    for n in subsample_levels
}

# Get metrics from BAM
print("Processing BAM...", flush=True)
with pysam.AlignmentFile(bam_path, "rb") as bam:
    for i, read in enumerate(bam):
        # Assign random number to this read
        r = random.random()  
                    
        # Get UMI (if applicable)
        if read.has_tag("PS") and read.has_tag("XS") and read.has_tag("BM") and read.has_tag("UB") and read.has_tag("XT"):
            if read.get_tag("PS") == 1 and read.get_tag("XS") == "Assigned":
                barcode = read.get_tag("BM")
                umi = read.get_tag("UB")
                gene = read.get_tag("XT")
                min_subsample = next((n for n in subsample_levels if r <= fractions[n]), subsample_levels[-1])

                all_umis.append((barcode, umi, gene, min_subsample))
        
        # Save to disk every 5 million reads
        if i % 5_000_000 == 0 and i > 0:
            print(f"Processed {i:,} reads...", flush=True)
            gc_return = gc.collect()
            
            df = pd.DataFrame(all_umis, columns=["barcode", "umi", "gene", "min_subsample"])
            df.to_parquet(os.path.join(output_path, "batches_tmp", f"batch_{batch_index}.parquet"), index=False, compression="snappy")
            print(f"Saved batch {batch_index} ({len(all_umis)} UMIs)", flush=True)
            batch_index += 1
            all_umis = []
                
        # Get per-read metrics
        for n in subsample_levels:
            if r > fractions.get(n, 1):
                continue

            m = metrics[n]
            
            # Total reads
            m["total_reads"] += 1

            # Q30 metrics for full read
            if read.query_qualities:
                m["q30_read_2_bases"] += sum(q >= 30 for q in read.query_qualities if q > 2)
                m["q30_read_2_total"] += sum(q > 2 for q in read.query_qualities)

            # Mapping quality metrics    
            if read.has_tag("XS") and read.get_tag("XS") == "Assigned":
                m["reads_Assigned"] += 1
            if read.has_tag("XS") and read.get_tag("XS") == "Unassigned_NoFeatures":
                m["reads_Unassigned_NoFeatures"] += 1
            if read.has_tag("XS") and read.get_tag("XS") == "Unassigned_MultiMapping":
                m["reads_Unassigned_MultiMapping"] += 1
            if read.has_tag("XS") and read.get_tag("XS") == "Unassigned_Ambiguity":
                m["reads_Unassigned_Ambiguity"] += 1
                
# Save remaining UMIs
if all_umis:
    print(f"Processed {total_reads} reads...", flush=True)
    df = pd.DataFrame(all_umis, columns=["barcode", "umi", "gene", "min_subsample"])
    df.to_parquet(os.path.join(output_path, "batches_tmp", f"batch_{batch_index}.parquet"), index=False, compression="snappy")
    print(f"Saved batch {batch_index} ({len(all_umis)} UMIs)", flush=True)

# Get all batch files
batch_files = glob.glob(os.path.join(output_path, "batches_tmp", "*.parquet"))

# Compute final metrics
rows = []

def safe_div(a, b): return a / b if b > 0 else 0

for i, n in enumerate(subsample_levels):
    gc_return = gc.collect()
    
    print(f"Computing metrics for subsample {i+1}/{len(subsample_levels)}.", flush=True)
    m = metrics[n]
    
    # Subsample UMIs
    unique_umis_n = Counter()
    
    for f in batch_files:
        table = pq.read_table(f, columns=["barcode", "umi", "gene", "min_subsample"])
        filtered_table = table.filter(pa.compute.less_equal(table["min_subsample"], n))
        for b, u, g in zip(filtered_table["barcode"].to_pylist(),
                            filtered_table["umi"].to_pylist(),
                            filtered_table["gene"].to_pylist()):
            unique_umis_n[(b,u,g)] += 1
        
    # Count UMIs
    total_umis = sum(unique_umis_n.values())
    unique_umis = len(unique_umis_n)
    
    tissue_unique_umis_n = Counter()
    for (barcode, umi, gene), count in unique_umis_n.items():
        if barcode in tissue_barcodes:
            tissue_unique_umis_n[(barcode, umi, gene)] = count
            
    tissue_unique_umis = len(tissue_unique_umis_n)
    
    # Per-gene UMIs in tissue
    gene_counts = Counter()
    tissue_gene_counts = Counter()
    unique_gene_counts = Counter()
    unique_gene_counts_tissue = Counter()

    for (barcode, umi, gene), count in unique_umis_n.items():
        gene_counts[gene] += count
        unique_gene_counts[gene] += 1
        if barcode in tissue_barcodes:
            tissue_gene_counts[gene] += count
            unique_gene_counts_tissue[gene] += 1
    
    all_genes = sorted(
        set(gene_counts.keys()) |
        set(tissue_gene_counts.keys()) |
        set(unique_gene_counts.keys()) |
        set(unique_gene_counts_tissue.keys()))

    gene_df = pd.DataFrame([
        {
            "gene": g,
            "total_counts": gene_counts.get(g, 0),
            "total_counts_in_tissue": tissue_gene_counts.get(g, 0),
            "unique_counts": unique_gene_counts.get(g, 0),
            "unique_counts_in_tissue": unique_gene_counts_tissue.get(g, 0),
            "total_counts_corrected": safe_div(gene_counts.get(g, 0), tissue_area_mm2),
            "total_counts_in_tissue_corrected": safe_div(tissue_gene_counts.get(g, 0), tissue_area_mm2),
            "unique_counts_corrected": safe_div(unique_gene_counts.get(g, 0), tissue_area_mm2),
            "unique_counts_in_tissue_corrected": safe_div(unique_gene_counts_tissue.get(g, 0), tissue_area_mm2)
        }
        for g in all_genes
    ])
    gene_df.to_csv(os.path.join(output_path, f"gene_counts_per_subsample/gene_counts_{n}.csv"), index=False)

    # Median UMIs per gene
    median_umis_per_gene = np.median(list(gene_counts.values())) if gene_counts else 0
    median_tissue_umis_per_gene = np.median(list(tissue_gene_counts.values())) if tissue_gene_counts else 0
    median_unique_umis_per_gene = np.median(list(unique_gene_counts.values())) if unique_gene_counts else 0
    median_unique_umis_per_gene_tissue = np.median(list(unique_gene_counts_tissue.values())) if unique_gene_counts_tissue else 0

    # Median reads per unique UMI
    median_reads_per_unique_umi = np.median(list(unique_umis_n.values())) if unique_umis else 0
    median_reads_per_unique_umi_tissue = np.median(list(tissue_unique_umis_n.values())) if tissue_unique_umis_n else 0

    # Create dataframe
    row = {
        "Q30_bases_read_2": (safe_div(m["q30_read_2_bases"], m["q30_read_2_total"])) * 100,
        "reads_Assigned": (safe_div(m["reads_Assigned"], m["total_reads"])) * 100,
        "reads_mapped_Unassigned_NoFeatures": (safe_div(m["reads_Unassigned_NoFeatures"], m["total_reads"])) * 100,
        "reads_mapped_Unassigned_MultiMapping": (safe_div(m["reads_Unassigned_MultiMapping"], m["total_reads"])) * 100,
        "reads_Unassigned_Ambiguity": (safe_div(m["reads_Unassigned_Ambiguity"], m["total_reads"])) * 100,
        "number_of_spots_with_at_least_1_UMI_in_tissue": len({b for b, _, _ in tissue_unique_umis_n}),
        "percentage_spots_recovered": safe_div(len({b for b, _, _ in tissue_unique_umis_n}), len(tissue_barcodes)) * 100,
        "total_number_of_reads": m["total_reads"],
        "total_number_of_reads_target": n,
        "total_number_of_genes": len(gene_counts),
        "total_number_of_UMIs": total_umis,
        "total_number_of_unique_UMIs": unique_umis,
        "number_of_genes_in_tissue": len(tissue_gene_counts),
        "number_of_unique_UMIs_in_tissue": tissue_unique_umis,
        "percentage_of_genes_in_tissue": safe_div(len(tissue_gene_counts), len(gene_counts)) * 100 if gene_counts else 0,
        "percentage_of_unique_UMIs_in_tissue": safe_div(tissue_unique_umis, unique_umis) * 100 if unique_umis else 0,
        "sequencing_saturation_reads": (1 - safe_div(unique_umis, m["total_reads"])) * 100 if m["total_reads"] > 0 else 0,
        "sequencing_saturation_umis": (1 - safe_div(unique_umis, total_umis)) * 100 if total_umis > 0 else 0,
        "unique_UMIs_per_mm2_in_tissue": safe_div(tissue_unique_umis, tissue_area_mm2),
        "median_UMIs_per_gene": median_umis_per_gene,
        "median_UMIs_per_gene_per_mm2": safe_div(median_umis_per_gene, tissue_area_mm2),
        "median_unique_UMIs_per_gene": median_unique_umis_per_gene,
        "median_unique_UMIs_per_gene_per_mm2": safe_div(median_unique_umis_per_gene, tissue_area_mm2),
        "median_unique_UMIs_per_gene_in_tissue": median_unique_umis_per_gene_tissue,
        "median_unique_UMIs_per_gene_in_tissue_per_mm2": safe_div(median_unique_umis_per_gene_tissue, tissue_area_mm2),
        "percentage_of_reads_with_usable_UMI": (safe_div(total_umis, m["total_reads"])) * 100,
        "ratio_reads_per_unique_UMI": safe_div(m["total_reads"], unique_umis),
        "median_reads_per_unique_UMI": median_reads_per_unique_umi,
        "median_reads_per_unique_UMI_in_tissue": median_reads_per_unique_umi_tissue
    }
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(output_path, "bam_metrics.csv"), index=False)
print("Metrics saved to ", os.path.join(output_path, "bam_metrics.csv"), flush=True)