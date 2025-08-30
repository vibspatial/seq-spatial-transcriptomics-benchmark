import pysam
import pandas as pd
import numpy as np
import random
import os
from collections import Counter

# User inputs
bam_path = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/SpaceRanger/outs/possorted_genome_bam.bam"
all_barcodes_path = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/harpy/all_barcodes.txt"
tissue_barcodes_path = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/harpy/barcodes_in_tissue.txt"
probe_set_path = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/SpaceRanger/outs/probe_set.csv"
tissue_area_mm2 = 15.46
step_size = 10_000_000  # 10M read increments per subsample
output_folder = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output"
random.seed(42)

# Create output folders
os.makedirs(output_folder, exist_ok=True)
os.makedirs(os.path.join(output_folder, "gene_counts_per_subsample"), exist_ok=True)

# Subsample BAM to a smaller size for testing
testing = True
if testing:
    bam_out = f"{output_folder}/bam_small.bam"
    num_reads = 1_000_000
    step_size = 200_000
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        with pysam.AlignmentFile(bam_out, "wb", template=bam) as out:
            for i, read in enumerate(bam):
                if i >= num_reads:
                    break
                out.write(read)
    print(f"Written {num_reads} reads to {bam_out}")
    
    bam_path = bam_out

# Load all barcodes
with open(all_barcodes_path) as f:
    all_barcodes = set(line.strip() for line in f)

# Load barcodes in tissue
with open(tissue_barcodes_path) as f:
    tissue_barcodes = set(line.strip() for line in f)
    
# Get included genes from probe set
probe_set_df = pd.read_csv(probe_set_path, comment='#')
included_genes = set(probe_set_df[probe_set_df['included'] == True]['gene_id'])

# Count total reads
print("Counting total reads...")
with pysam.AlignmentFile(bam_path, "rb") as bam:
    total_reads = sum(1 for _ in bam)
print(f"Total reads: {total_reads:,}")

# Define subsample levels
subsample_levels = list(range(0, total_reads, step_size))
if subsample_levels[-1] != total_reads:
    subsample_levels.append(total_reads)
fractions = {n: n / total_reads for n in subsample_levels if n > 0}

# Initialize metrics
metrics = {
    n: {
        "total_reads": 0,
        "total_umis": 0,
        "unique_umis": set(),
        "genes": set(),
        "barcodes": set(),
        "q30_read_bases": 0,
        "q30_read_total": 0,
        "q30_barcode_bases": 0,
        "q30_barcode_total": 0,
        "q30_umi_bases": 0,
        "q30_umi_total": 0,
        "reads_mapped": 0,
        "reads_confident": 0,
        "reads_confident_filtered": 0,
        "reads_half_mapped": 0,
        "reads_split_mapped": 0,
        "tissue": {
            "umis": 0,
            "unique_umis": set(),
            "genes": set(),
            "barcodes": set(),
            "gene_counts": {}
        },
        "umi_counts": [],
        "umi_counts_tissue": []
    }
    for n in subsample_levels
} # TODO: Add mapping metrics

# Get metrics from BAM
print("Processing BAM...")
with pysam.AlignmentFile(bam_path, "rb") as bam:
    for i, read in enumerate(bam):
        if i % 5_000_000 == 0 and i > 0:
            print(f"Processed {i:,} reads...")

        r = random.random()  # assign random number to this read

        for n in subsample_levels:
            if n == 0 or r > fractions.get(n, 1):
                continue

            m = metrics[n]
            
            # Total reads
            m["total_reads"] += 1

            # Q30 metrics for read, barcode and UMI
            if read.query_qualities:
                m["q30_read_bases"] += sum(q >= 30 for q in read.query_qualities if q > 2)
                m["q30_read_total"] += sum(q > 2 for q in read.query_qualities)

            if read.has_tag("CY"):
                q_barcode = read.get_tag("CY")
                q_barcode_scores = [ord(c) - 33 for c in q_barcode]
                m["q30_barcode_bases"] += sum(q >= 30 for q in q_barcode_scores if q > 2)
                m["q30_barcode_total"] += sum(q > 2 for q in q_barcode_scores)

            if read.has_tag("UY"):
                q_umi = read.get_tag("UY")
                q_umi_scores = [ord(c) - 33 for c in q_umi]
                m["q30_umi_bases"] += sum(q >= 30 for q in q_umi_scores if q > 2)
                m["q30_umi_total"] += sum(q > 2 for q in q_umi_scores)

                
            # Mapping quality metrics   
            UMI_count = False 
            mapq = read.mapping_quality 
            if mapq > 0:
                m["reads_mapped"] += 1
            if mapq == 255:
                m["reads_confident"] += 1
                if read.has_tag("GX") and read.get_tag("GX") in included_genes:
                    m["reads_confident_filtered"] += 1
                    UMI_count = True
            if mapq == 1:
                m["reads_half_mapped"] += 1
            if mapq == 3:
                m["reads_split_mapped"] += 1
                
            # Barcodes, UMIs, Genes
            if UMI_count and read.has_tag("CB") and read.has_tag("UB") and read.has_tag("GN"):
                barcode = read.get_tag("CB")
                umi = read.get_tag("UB")
                gene = read.get_tag("GN")

                m["total_umis"] += 1
                m["unique_umis"].add((barcode, umi, gene))
                m["genes"].add(gene)
                m["barcodes"].add(barcode)

                # Count per-UMI reads
                m["umi_counts"].append(umi)

                # Tissue-specific metrics
                if barcode in tissue_barcodes:
                    m["tissue"]["umis"] += 1
                    m["tissue"]["unique_umis"].add((barcode, umi, gene))
                    m["tissue"]["genes"].add(gene)
                    m["tissue"]["barcodes"].add(barcode)
                    m["umi_counts_tissue"].append(umi)
                    m["tissue"]["gene_counts"][gene] = m["tissue"]["gene_counts"].get(gene, 0) + 1

# Compute final metrics
rows = []
for n in subsample_levels:
    m = metrics[n]

    total_umis = m["total_umis"]
    unique_umis = len(m["unique_umis"])
    tissue_umis = m["tissue"]["umis"]
    tissue_unique_umis = len(m["tissue"]["unique_umis"])

    def safe_div(a, b): return a / b if b > 0 else 0
    
    # Per-gene UMIs in tissue
    gene_counts = m["tissue"]["gene_counts"]
    gene_df = pd.DataFrame([
        {"gene": g, "total_counts": c, "total_counts_corrected": c / tissue_area_mm2}
        for g, c in gene_counts.items()
    ])
    gene_df.to_csv(os.path.join(output_folder, f"gene_counts_per_subsample/gene_counts_{n}.csv"), index=False)
    
    # Median & mean UMIs per gene
    median_umis_per_gene = np.median(list(gene_counts.values())) if gene_counts else 0
    mean_umis_per_gene = np.mean(list(gene_counts.values())) if gene_counts else 0

    row = {
        "Q30_bases_read": (safe_div(m["q30_read_bases"], m["q30_read_total"])) * 100,
        "Q30_bases_spatial_barcode": (safe_div(m["q30_barcode_bases"], m["q30_barcode_total"])) * 100,
        "Q30_bases_UMI": (safe_div(m["q30_umi_bases"], m["q30_umi_total"])) * 100,
        "reads_mapped": (safe_div(m["reads_mapped"], m["total_reads"])) * 100,
        "reads_mapped_confidently": (safe_div(m["reads_confident"], m["total_reads"])) * 100,
        "reads_mapped_confidently_to_filtered": (safe_div(m["reads_confident_filtered"], m["total_reads"])) * 100,
        "reads_half_mapped": (safe_div(m["reads_half_mapped"], m["total_reads"])) * 100,
        "reads_split_mapped": (safe_div(m["reads_split_mapped"], m["total_reads"])) * 100,
        "number_of_spots_with_at_least_1_UMI_in_tissue": len({b for b,_,_ in m["tissue"]["unique_umis"]}),
        "percentage_spots_recovered": (safe_div(len(m["tissue"]["barcodes"]), len(tissue_barcodes))) * 100,
        "total_number_of_reads": m["total_reads"],
        "total_number_of_reads_target": n,
        "total_number_of_genes": len(m["genes"]),
        "total_number_of_UMIs": total_umis,
        "total_number_of_unique_UMIs": unique_umis,
        "number_of_genes_in_tissue": len(m["tissue"]["genes"]),
        "number_of_unique_UMIs_in_tissue": tissue_unique_umis,
        "percentage_of_genes_in_tissue": (safe_div(len(m["tissue"]["genes"]), len(m["genes"]))) * 100,
        "percentage_of_unique_UMIs_in_tissue": (safe_div(tissue_unique_umis, unique_umis)) * 100,
        "sequencing_saturation_reads": (1 - safe_div(unique_umis, total_umis)) * 100,
        "sequencing_saturation_umis": (1 - safe_div(unique_umis, m["total_reads"])) * 100,
        "unique_UMIs_per_mm2_in_tissue": (safe_div(tissue_unique_umis, tissue_area_mm2)) * 100,
        "median_UMIs_per_gene": median_umis_per_gene,
        "median_UMIs_per_gene_per_mm2": safe_div(median_umis_per_gene, tissue_area_mm2),
        "mean_UMIs_per_gene": mean_umis_per_gene,
        "mean_UMIs_per_gene_per_mm2": safe_div(mean_umis_per_gene, tissue_area_mm2),
        "percentage_of_reads_with_usable_UMI": (safe_div(total_umis, m["total_reads"])) * 100,
        "ratio_reads_per_unique_UMI": safe_div(m["total_reads"], unique_umis),
        "median_reads_per_UMI": np.median(list(Counter(m["umi_counts"]).values())) if Counter(m["umi_counts"]) else 0,
        "median_reads_per_UMI_in_tissue": np.median(list(Counter(m["umi_counts_tissue"]).values())) if Counter(m["umi_counts_tissue"]) else 0
    }
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(output_folder, "subsample_metrics.csv"), index=False)
print("Metrics saved to ", os.path.join(output_folder, "subsample_metrics.csv"))