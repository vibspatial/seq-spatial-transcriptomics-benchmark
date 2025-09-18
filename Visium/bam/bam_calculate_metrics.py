import os
import uuid
import time
import shutil


import dask.dataframe as dd
import pandas as pd
from pathlib import Path
import argparse


def calculate_metrics(
        input_path : str | Path,
        output_dir : str | Path,
        tissue_barcodes_path: str | Path,
        tissue_area: int = 1,
        sample_size: int = 0,

):
    os.makedirs(output_dir, exist_ok=True )
    df=dd.read_csv(
    input_path,
    dtype={
        "random_number" : "float32",
        "included": "bool",
        "included_tissue": "bool",
        "barcode": "string",
        "umi": "string",
        "gene": "string",
        "q30_read_bases": "int32",
        "q30_read_total": "int32",
        "q30_barcode_bases": "int32",
        "q30_barcode_total": "int32",
        "q30_umi_bases": "int32",
        "q30_umi_total": "int32",
        "reads_mapped": "bool",
        "reads_confident": "bool",
        "reads_confident_filtered": "bool",
        "reads_half_mapped": "bool",
        "reads_split_mapped": "bool",
    },
    #blocksize="1GB
    )

    with open(tissue_barcodes_path) as f:
        tissue_barcodes = set(line.strip() for line in f)

    # sample
    if sample_size:
        start = time.time()
        frac = sample_size / len( df )
        print( f"Start sampling {sample_size} reads.")
        df=df.sample( frac=frac )
        _parquet_path = os.path.join( output_dir, f"sample_{uuid.uuid4()}_{sample_size}.parquet"  )
        df.to_parquet( _parquet_path )
        # read back lazy
        df=dd.read_parquet( _parquet_path )
        print("End sampling.")
        print( f"took: {time.time() - start} seconds" )
        

    # compute metrics
    df_included=df[ df[ "included" ] ]
    print( "Start grouping UMIs:" )
    start=time.time()
    df_included_grouped=df_included.groupby( [ "barcode", "umi", "gene" ] ).size().compute()
    print( "End grouping UMIs." )
    print( f"took: {time.time() - start} seconds" )
    print( "Calculate total UMIs." )
    total_umis=df_included_grouped.sum()
    print( "Calculate unique UMIs" )
    unique_umis=len(df_included_grouped)

    print( "total UMIs:", total_umis )
    print( "total unique UMIs:", unique_umis )


    df_included_tissue = df[df[ "included_tissue" ]]
    print( "Start grouping UMIs in tissue." )
    start = time.time()
    df_included_tissue_grouped=df_included_tissue.groupby( [ "barcode", "umi", "gene" ] ).size().compute()
    print("End grouping UMIs in tissue." )
    print( f"took: {time.time() - start} seconds" )
    print( "Calculate total unique UMIs in tissue." )
    tissue_unique_umis = len( df_included_tissue_grouped )

    print( "total unique UMIs in tissue:", tissue_unique_umis  )


    print( "Start calculating gene counts." )
    gene_counts=df_included_grouped.groupby( "gene" ).sum()
    print( "End calculating gene counts." )
    print( "Start calculating unique gene counts." )
    unique_gene_counts = df_included_grouped.groupby( "gene" ).count()
    print( "End calculating unique gene counts." )

    # Same, but in tissue
    print( "Start calculating gene counts in tissue." )
    gene_counts_tissue=df_included_tissue_grouped.groupby( "gene" ).sum()
    print( "End calculating gene counts in tissue." )
    print( "Start calculating unique gene counts in tissue." )
    unique_gene_counts_tissue = df_included_tissue_grouped.groupby( "gene" ).count()  # 3s for 2M 18s for 5M
    print( "End calculating unique gene counts in tissue." )

    gene_counts=gene_counts.to_frame( name="gene_counts" ).reset_index()
    unique_gene_counts=unique_gene_counts.to_frame( name="unique_gene_counts" ).reset_index()

    gene_counts_tissue=gene_counts_tissue.to_frame( name="gene_counts_tissue" ).reset_index()
    unique_gene_counts_tissue=unique_gene_counts_tissue.to_frame( name="unique_gene_counts_tissue" ).reset_index()

    gene_counts_merged = (
        gene_counts
        .merge(unique_gene_counts, on='gene', how='outer')
        .merge(gene_counts_tissue, on='gene', how='outer')
        .merge(unique_gene_counts_tissue, on='gene', how='outer')
    )

    assert tissue_area>0, "Please provide tissue_area>0."

    for _name in [ "gene_counts", "unique_gene_counts", "gene_counts_tissue", "unique_gene_counts_tissue" ]:
        gene_counts_merged[ f"{_name}_corrected" ] = gene_counts_merged[ _name ] / tissue_area


    gene_counts_merged.to_csv( os.path.join( output_dir, f"gene_counts_sample_size_{sample_size}.csv" ) )

    median_gene_counts=gene_counts_merged[ "gene_counts" ].median()
    median_unique_gene_counts=gene_counts_merged["unique_gene_counts" ].median()
    median_gene_counts_tissue=gene_counts_merged[ "gene_counts_tissue" ].median()
    median_unique_gene_counts_tissue=gene_counts_merged[ "unique_gene_counts_tissue" ].median()
    print( "Median gene counts (total; unique; total in tissue; unique in tissue):", median_gene_counts, median_unique_gene_counts, median_gene_counts_tissue, median_unique_gene_counts_tissue )

    median_reads_per_unique_umi=df_included_grouped.median()  # Am I sure this is correct. Check with Julien what is definition of unique UMI
    median_reads_per_unique_umi_tissue=df_included_tissue_grouped.median()

    print("Median reads per unique umi:", median_reads_per_unique_umi, median_reads_per_unique_umi_tissue )


    print( "Start calculating metrics on all reads." )
    start = time.time()
    # metrics to calculate for all reads. Check with Julien, this is for all reads, right?
    columns_to_sum = [
        "q30_read_bases",
        "q30_read_total",
        "q30_barcode_bases",
        "q30_barcode_total",
        "q30_umi_bases",
        "q30_umi_total",
        "reads_mapped",
        "reads_confident",
        "reads_confident_filtered",
        "reads_half_mapped",
        "reads_split_mapped",
        ]

    reads_sum=df[  columns_to_sum  ].sum().compute()
    print( "End Calculating metrics on all reads." )
    print( "Took:", time.time() - start )

    def safe_div(a, b): return a / b if b > 0 else 0

    total_reads = len( df )

    assert total_reads > 0, "Total number of reads is zero."

    # put it in dictionary, then write to csv
    row = {
        "Q30_bases_read":safe_div(reads_sum["q30_read_bases"], reads_sum["q30_read_total"])*100,
        "Q30_bases_spatial_barcode": safe_div(reads_sum["q30_barcode_bases"], reads_sum["q30_barcode_total"]) * 100,
        "Q30_bases_UMI" : safe_div(reads_sum["q30_umi_bases"], reads_sum["q30_umi_total"]) * 100,
        "reads_mapped" : safe_div(reads_sum["reads_mapped"], total_reads) * 100,
        "reads_mapped_confident" : safe_div(reads_sum["reads_confident"], total_reads ) * 100,
        "reads_mapped_confident_to_filtered" : safe_div(reads_sum["reads_confident_filtered"], total_reads ) * 100,
        "reads_half_mapped" : safe_div(reads_sum["reads_half_mapped"], total_reads) * 100,
        "reads_split_mapped" : safe_div(reads_sum["reads_split_mapped"], total_reads) * 100,
        "number_of_spots_with_at_least_1_UMI_in_tissue": df_included_tissue_grouped.index.get_level_values("barcode").nunique(),
        "percentage_spots_recovered":df_included_tissue_grouped.index.get_level_values("barcode").nunique()/len( tissue_barcodes )*100,
        "total_number_of_reads" : total_reads,
        "total_number_of_genes" : len( gene_counts ),
        "total_number_of_UMIs": total_umis,
        "total_number_of_unique_UMIs": unique_umis,
        "number_of_genes_in_tissue": len(gene_counts_tissue),
        "number_of_unique_UMIs_in_tissue": tissue_unique_umis,
        "percentage_of_genes_in_tissue": safe_div( len( gene_counts_tissue ) ,len( gene_counts ) )*100,
        "percentage_of_unique_UMIs_in_tissue": safe_div(  tissue_unique_umis, unique_umis )*100,
        "sequencing_saturation_reads" : (1-safe_div( unique_umis, total_reads ))*100 if total_reads>0 else 0,
        "sequencing_saturation_UMIs": (1-safe_div( unique_umis, total_umis ))*100,
        "unique_UMIs_per_mm2_in_tissue": safe_div( tissue_unique_umis, tissue_area ),
        "median_UMIs_per_gene": median_gene_counts,
        "median_UMIs_per_gene_per_mm2": safe_div( median_gene_counts,tissue_area ),
        "median_unique_UMIs_per_gene": median_unique_gene_counts,
        "median_unique_UMIs_per_gene_per_mm2": safe_div( median_unique_gene_counts, tissue_area ),
        "median_unique_UMIs_per_gene_in_tissue": median_unique_gene_counts_tissue,
        "median_unique_UMIs_per_gene_in_tissue_per_mm2": safe_div(median_unique_gene_counts_tissue, tissue_area),
        "percentage_of_reads_with_usable_UMI": safe_div( total_umis, total_reads)*100,
        "ratio_reads_per_unique_UMI": safe_div( total_reads, unique_umis ),
        "median_reads_per_unique_UMI": median_reads_per_unique_umi,
        "median_reads_per_unique_UMI_tissue": median_reads_per_unique_umi_tissue,
    }

    # clean up parquet
    if sample_size>0 and os.path.exists(_parquet_path):
        shutil.rmtree(_parquet_path)

    pd.DataFrame( [ row ] ).to_csv(
        os.path.join( output_dir, f"metrics_{sample_size}.csv" ),
        index=False,
        )

if __name__ == "__main__":

    # User inputs
    parser = argparse.ArgumentParser(description="BAM metrics calculation")
    parser.add_argument("--input_path", required=True, help="Path to csv file created by bam_preprocess.py")
    parser.add_argument("--tissue_barcodes_path", required=True, help="Path to tissue_barcodes.txt")
    parser.add_argument("--output_dir", required=True, help="Output directory where metrics and intermediate files will be stored.")
    parser.add_argument("--tissue_area", type=int, required=False, default=1, help="Tissue area in mm2.")
    parser.add_argument("--sample_size", type=int, required=False, default=0, help="If >0 specified will calculate metrics for a random subset of 'sample size' reads.")
    args = parser.parse_args()

    calculate_metrics(
        input_path=args.input_path,
        tissue_barcodes_path=args.tissue_barcodes_path,
        output_dir = args.output_dir,
        tissue_area=args.tissue_area,
        sample_size=args.sample_size,
        )