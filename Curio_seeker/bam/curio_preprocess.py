import os
import time
import argparse
from pathlib import Path

import pandas as pd
import Levenshtein
import pyarrow.parquet as pq

import random
random.seed(42)

def main(
        parquet_dir : str | Path,
        tissue_barcodes_path: str | Path,
        output_path: str| Path,
        size: int = 0,
):
    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_path} exists.")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    print("Reading parquet files...", flush=True)
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
    dedup = df.drop_duplicates(subset=["read1_id"], keep="first").reset_index(drop=True) # Only keep the first occurrence (all duplicated reads should be multimapping). This is somewhat risky.
    print(f"Number of unique reads: {len(dedup):,}", flush=True)

    # Load barcodes in tissue
    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)
        
    # UMI Deduplication
    print("\nDeduplicating UMIs")
    def umi_dedup_directional(df: pd.DataFrame, edit_distance: int = 1) -> pd.DataFrame:
        corrected = {}

        for (spatial_barcode, gene), subdf in df.groupby(["BM", "XT"]):
            umi_counts = subdf["molecular_barcode"].value_counts().to_dict()
            sorted_umis = sorted(umi_counts.items(), key=lambda x: x[1], reverse=True)

            assigned = {}

            for umi, count in sorted_umis:
                if umi in assigned: # The most abundant UMI is considered the correct UMI. Directional.
                    continue
                assigned[umi] = umi

                for other_umi, _ in sorted_umis:
                    if other_umi in assigned:
                        continue
                    if Levenshtein.distance(umi, other_umi) <= edit_distance: # Equal to Hamming for strings of the same length
                        assigned[other_umi] = umi
            
            for read_id, umi in zip(subdf["read1_id"], subdf["molecular_barcode"]):
                corrected[read_id] = assigned[umi]

        df["corrected_umi"] = df["read1_id"].map(corrected)
        return df

    mask = (dedup['XS'] == 'Assigned') & (dedup['r1_proper_structure_matched'] == True)
    df_for_umi_assigned = dedup[mask].copy()

    deduped = umi_dedup_directional(df_for_umi_assigned, edit_distance=1) # Hamming distance < 2
    
    read_to_corrected_umi = pd.Series(deduped["corrected_umi"].values, index=deduped['read1_id']).to_dict()
    dedup['corrected_umi'] = dedup['read1_id'].map(read_to_corrected_umi)

    # Create preprocessed CSV
    print(f"\nCreating preprocessed CSV", flush=True)
    start = time.time()
    with open(output_path, "w") as f_out:
        f_out.write("index,included,included_tissue,barcode,umi,gene,"
                    "reads_Assigned,reads_Unassigned_NoFeatures,reads_Unassigned_MultiMapping,reads_Unassigned_Ambiguity,reads_NA,random\n")
        
        for i, row in dedup.iterrows():
            r = random.random()
            step = 5000000
            if i%step == 0:
                print(f"{i:,} reads processed: {round((time.time() - start) / 60, 2)} minutes", flush=True)
            if size:
                if i > size:
                    break

            barcode = 0
            umi = 0
            gene = 0
            included = 0
            included_tissue = 0
            
            if row['r1_proper_structure_matched'] == True and row['XS'] == 'Assigned' and row['BM'] in tissue_barcodes:
                # For Curio Seeker, all barcodes are in tissue. 
                # We know there is some filtering step that removes a subset of barcodes even though they match the barcode whitelist, but we don't know why. This mimics this step.
                barcode = row['BM']
                umi = row['corrected_umi']
                gene = row['XT']
                    
            if barcode and umi and gene:
                included = 1 # for easier subsampling
                if barcode in tissue_barcodes:
                    included_tissue = 1

            # featureCounts metrics
            reads_Assigned = 0
            reads_Unassigned_NoFeatures = 0
            reads_Unassigned_MultiMapping = 0
            reads_Unassigned_Ambiguity = 0
            reads_NA = 0
            
            if row['XS'] == 'Assigned': # Read successfully maps to a feature
                reads_Assigned = 1
            if row['XS'] == 'Unassigned_NoFeatures': # Read maps to the genome but does not overlap any annotated feature (could be in intergenic region)
                reads_Unassigned_NoFeatures = 1
            if row['XS'] == 'Unassigned_MultiMapping': # Read maps to multiple locations in the genome
                reads_Unassigned_MultiMapping = 1
            if row['XS'] == 'Unassigned_Ambiguity': # Read maps to multiple features and cannot be assigned
                reads_Unassigned_Ambiguity = 1
            if pd.isna(row['XS']): # Read is unmapped
                reads_NA = 1
                
            f_out.write(f"{i},{included},{included_tissue},{barcode},{umi},{gene},"
                        f"{reads_Assigned},{reads_Unassigned_NoFeatures},{reads_Unassigned_MultiMapping},{reads_Unassigned_Ambiguity},{reads_NA},{r}\n")
    
    print(f"{i:,} reads processed: {round((time.time() - start) / 60, 2)} minutes", flush=True)        
    print(f"Done!", flush=True)
                
if __name__ == "__main__":

    # User inputs
    parser = argparse.ArgumentParser(description="Curio Seeker preprocessing")
    parser.add_argument("--parquet_dir", required=True, help="Directory containing parquet files")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--output_path", required=True, help="Output path. Preprocessed csv file will be stored there.")
    parser.add_argument("--size", type=int, required=False, default=0, help="If specified, will only preprocess the first 'size' reads.")

    args = parser.parse_args()

    main(
        parquet_dir=args.parquet_dir,
        tissue_barcodes_path=args.tissue_barcodes_path,
        output_path=args.output_path,
        size=args.size,
    )



