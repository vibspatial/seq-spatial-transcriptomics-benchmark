### Import packages

import sparrow as sp
import spatialdata as sd
import pandas as pd
import numpy as np
import dask as da
import os
import gc
import geopandas as gpd
from matplotlib import pyplot as plt
import seaborn as sns
import h5py
import json

from sparrow.io._visium_hd import visium_hd

from spatialdata.transformations import get_transformation
from sparrow.image._image import _get_spatial_element

import scanpy as sc
from sparrow.table import add_table_layer
from sparrow.points import add_points_layer
from sparrow.shape import add_shapes_layer
from sparrow.utils._keys import _INSTANCE_KEY, _REGION_KEY

import dask.dataframe as dd

from skimage import io

### Create sdata

sdata = visium_hd(
    path='/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs', 
    bin_size=[2,8,16], 
    dataset_id="", # Had to adapt spatialdata io visium reader according to https://github.com/scverse/spatialdata-io/issues/218
    bins_as_squares=True, 
    fullres_image_file='/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/IMAGES/SCA002_HE_hires.tif', # Adapted harpy visium hd wrapper to take this argument, I also had to set scale_factors=None in line 361
    output= None) # I run into buffer issues if I try to write the sdata to disk
# I set filtered_counts_file=False in harpy visium hd wrapper

se=_get_spatial_element(sdata, layer="_full_image")

arr=se.data.rechunk("auto")
arr_reshaped = arr.transpose(1, 2, 0)

sdata=sp.im.add_image_layer(sdata, arr=arr_reshaped, transformations=get_transformation(sdata[ "_full_image" ], get_all=True), output_layer="_full_image", overwrite=True, scale_factors=(2, 2, 2, 2))
sdata.write('/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/sdata.zarr')

### Calculate QC metrics

for table_layer in ['square_002um', 'square_008um', 'square_016um']:
    adata = sdata.tables[table_layer].copy()
    sc.pp.calculate_qc_metrics(adata, layer=None, inplace=True)
    
    sdata = add_table_layer(
        sdata,
        adata=adata,
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY][0],
        overwrite=True,
    )
    
### Convert transcripts to image for Bigwarp
    
def convert_transcripts_to_img(adata, x, y, intensity_column, output_path):
    # Create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Extract the x, y, and intensity columns
    x = adata.obs[x].values
    y = adata.obs[y].values
    intensity = adata.obs[intensity_column].values

    # Scale the intensity values to the 8-bit range [0, 255]
    intensity_scaled = 255 * (intensity - intensity.min()) / (intensity.max() - intensity.min())
    intensity_scaled = intensity_scaled.astype(np.uint8)

    # Create an empty image
    width = x.max()+1
    height = y.max()+1
    img = np.zeros((height, width), dtype=np.uint8)

    # Assign the scaled intensity values to the image xy-positions
    img[y, x] = intensity_scaled
    
    # Save image
    io.imsave(output_path, img)
    
for col in ['n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'log1p_n_genes_by_counts']:
    convert_transcripts_to_img(
        adata = sdata.tables['square_002um'],
        x = 'array_col',
        y = 'array_row',
        intensity_column = col, 
        output_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/transcript_images/square_002um_{col}.tiff')
    
### Add Bigwarp landmarks

def add_landmarks(sdata, landmarks_csv, x_column, y_column, fullres_microns_per_pixel, points_layer):
    df_landmarks = pd.read_csv(landmarks_csv, header=None)
    df_pxl = pd.DataFrame({
            'name': df_landmarks[0],
            'bigwarp_x': df_landmarks[x_column] * ((1 / fullres_microns_per_pixel) * 1000),
            'bigwarp_y': df_landmarks[y_column] * ((1 / fullres_microns_per_pixel) * 1000),
        })
    ddf_landmarks = dd.from_pandas(df_pxl, npartitions=2)

    sdata = add_points_layer(
        sdata, 
        ddf = ddf_landmarks, 
        output_layer = points_layer, 
        coordinates = {"x": "bigwarp_x", "y": "bigwarp_y"}, 
        transformations=None, 
        overwrite=True)
    
    return sdata

sdata = add_landmarks(
    sdata, 
    landmarks_csv = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/2_cytassist_to_slidescanner/landmarks.csv',
    x_column = 4,
    y_column = 5,
    fullres_microns_per_pixel = 0.22009446422471257,
    points_layer = 'landmarks_JM'
)

sdata = add_landmarks(
    sdata, 
    landmarks_csv = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/3_counts_to_slidescanner/landmarks.csv',
    x_column = 4,
    y_column = 5,
    fullres_microns_per_pixel = 0.22009446422471257,
    points_layer = 'landmarks_counts_JM'
)

sdata = add_landmarks(
    sdata, 
    landmarks_csv = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/4_cytassist_to_slidescanner_KH/landmarks_KH.csv',
    x_column = 4,
    y_column = 5,
    fullres_microns_per_pixel = 0.121971056, # NOTE: landmarks_KH.csv was created by using the tissue_hires_image.png as fixed image (instead of the SCA002_HE_hires.tif) so we need to take the scaling factor into account.
    points_layer = 'landmarks_KH'
)

### Add Cytassist coordinates

def add_cytassist_coordinates(sdata, table_layer, tissue_positions_path, feature_slice_path, bin_size, spatial_key):
    # Copy adata
    adata = sdata.tables[table_layer].copy()
    
    # Load the tissue positions
    df = pd.read_parquet(tissue_positions_path).set_index('barcode')
    
    # Load the metadata and extract affine transformation matrix
    with h5py.File(feature_slice_path, "r") as f5:
        metadata = json.loads(dict(f5.attrs)["metadata_json"])
    affine_matrix = np.array(metadata['transform_matrices']['spot_colrow_to_cytassist_colrow'])

    # Apply the affine transformation matrix
    if bin_size == 2:
        coords = np.hstack([df[['array_col', 'array_row']].values, np.ones((df.shape[0], 1))])
    elif bin_size > 2:
        coords = np.hstack([df[['array_col', 'array_row']].values * (bin_size/2), np.ones((df.shape[0], 1))])
    if bin_size < 2:
        raise ValueError("bin_size must be 2 or greater.")
    
    transformed_coords = coords @ affine_matrix.T

    # Save the transformed coordinates
    df_transformed = pd.DataFrame({
        'x': transformed_coords[:, 0],
        'y': transformed_coords[:, 1],
    }, index=df.index)
    
    # Add bigwarp coordinates temporarily to adata.obs (to make sure the barcodes match)
    adata.obs = adata.obs.join(df_transformed, how='left')
    
    # Add coordinates to obsm
    coordinates = adata.obs[['x', 'y']].to_numpy()
    adata.obsm[spatial_key] = coordinates
    adata.obs = adata.obs.drop(['x', 'y'], axis=1)
    
    # Overwrite table layer
    sdata = add_table_layer(
        sdata,
        adata=adata,
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY][0],
        overwrite=True,
    )
    
    return sdata

for bin_size in [2, 8, 16]:
    table_layer = f'square_{bin_size:03}um'

    sdata = add_cytassist_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs/binned_outputs/{table_layer}/spatial/tissue_positions.parquet', 
        feature_slice_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs/feature_slice.h5', 
        bin_size = bin_size, 
        spatial_key = 'spatial_cytassist')
    
    sdata = add_cytassist_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311/outs/binned_outputs/{table_layer}/spatial/tissue_positions.parquet', 
        feature_slice_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311/outs/feature_slice.h5', 
        bin_size = bin_size, 
        spatial_key = 'spatial_cytassist_311_json')
    
    sdata = add_cytassist_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311_no_json/outs/binned_outputs/{table_layer}/spatial/tissue_positions.parquet', 
        feature_slice_path = f'/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311_no_json/outs/feature_slice.h5', 
        bin_size = bin_size, 
        spatial_key = 'spatial_cytassist_311_no_json')
    
### Add Bigwarp registration results

def add_bigwarp_coordinates(sdata, table_layer, tissue_positions_path, bigwarp_csv_path, fullres_microns_per_pixel, spatial_key):
    gc.collect()
    
    # Copy adata
    adata = sdata.tables[table_layer].copy()
    
    # Load the tissue positions
    df_tissue = pd.read_parquet(tissue_positions_path)
    
    # Read the CSV file containing transformed coordinates from bigwarp
    df_bigwarp = pd.read_csv(bigwarp_csv_path)

    # Convert coordinates to pixels
    scale_factor = (1 / fullres_microns_per_pixel) * 1000
    df_bigwarp[['x', 'y']] = df_bigwarp[['x', 'y']].to_numpy() * scale_factor

    # Add the new coordinates to the tissue positions data
    df_tissue = (df_tissue
                .join(df_bigwarp)
                .set_index('barcode')
                .filter(items=['x', 'y']))

    # Add bigwarp coordinates temporarily to adata.obs (to make sure the barcodes match)
    adata.obs = adata.obs.join(df_tissue, how='left')
    
    # Add coordinates to obsm
    adata.obsm[spatial_key] = adata.obs[['x', 'y']].to_numpy()
    adata.obs = adata.obs.drop(['x', 'y'], axis=1)
    
    # Reduce memory usage before overwriting table layer
    gc.collect()
    del df_tissue
    del df_bigwarp
    
    # Overwrite table layer
    sdata = add_table_layer(
        sdata,
        adata=adata,
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY][0],
        overwrite=True,
    )
    
    return sdata

paths_dict = {
    'square_002um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/2_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/2_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/3_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/3_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/4_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_002um/4_bigwarp_tps_transformed_coordinates.csv'],
    'square_008um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/2_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/2_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/3_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/3_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/4_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_008um/4_bigwarp_tps_transformed_coordinates.csv'],
    'square_016um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/Visium_HD/1/SCA002_reseq/COUNTS/outs/binned_outputs/square_016um/spatial/tissue_positions.parquet',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/2_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/2_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/3_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/3_bigwarp_tps_transformed_coordinates.csv',
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/4_bigwarp_affine_transformed_coordinates.csv', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Bigwarp/transformed_barcode_positions/square_016um/4_bigwarp_tps_transformed_coordinates.csv'],
}

for table_layer, paths in paths_dict.items():

    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[1], 
        fullres_microns_per_pixel = 0.22009446422471257, 
        spatial_key = 'spatial_bigwarp_affine_JM'
    )
    
    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[2], 
        fullres_microns_per_pixel = 0.22009446422471257,
        spatial_key = 'spatial_bigwarp_tps_JM'
    )
    
    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[3], 
        fullres_microns_per_pixel = 0.22009446422471257, 
        spatial_key = 'spatial_bigwarp_affine_counts_JM'
    )
    
    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[4], 
        fullres_microns_per_pixel = 0.22009446422471257,
        spatial_key = 'spatial_bigwarp_tps_counts_JM'
    )
    
    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[5], 
        fullres_microns_per_pixel = 0.121971056, # NOTE: landmarks_KH.csv was created by using the tissue_hires_image.png as fixed image (instead of the SCA002_HE_hires.tif) so we need to take the scaling factor into account.
        spatial_key = 'spatial_bigwarp_affine_KH'
    )
    
    sdata = add_bigwarp_coordinates(
        sdata, 
        table_layer = table_layer, 
        tissue_positions_path = paths[0],
        bigwarp_csv_path = paths[6], 
        fullres_microns_per_pixel = 0.121971056, # NOTE: landmarks_KH.csv was created by using the tissue_hires_image.png as fixed image (instead of the SCA002_HE_hires.tif) so we need to take the scaling factor into account.
        spatial_key = 'spatial_bigwarp_tps_KH'
    )
    
### Add Spaceranger registration results    
    
def add_spaceranger_coordinates(sdata, table_layer, new_tissue_positions_path, spatial_key):
    # Copy adata
    adata = sdata.tables[table_layer].copy()
    
    # Read the CSV file containing transformed coordinates from bigwarp
    df_spaceranger = pd.read_parquet(new_tissue_positions_path)
    
    # Add the new coordinates to the tissue positions data
    df_spaceranger = (df_spaceranger
                .set_index('barcode')
                .filter(items=['pxl_col_in_fullres', 'pxl_row_in_fullres']))

    # Add spaceranger coordinates temporarily to adata.obs (to make sure the barcodes match)
    adata.obs = adata.obs.join(df_spaceranger, how='left')
    
    # Add coordinates to obsm
    coordinates = adata.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy()
    adata.obsm[spatial_key] = coordinates
    adata.obs = adata.obs.drop(['pxl_col_in_fullres', 'pxl_row_in_fullres'], axis=1)
    
    # Overwrite table layer
    sdata = add_table_layer(
        sdata,
        adata=adata,
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY][0],
        overwrite=True,
    )
    
    return sdata

paths_dict = {
    'square_002um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311_no_json/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet'],
    'square_008um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311_no_json/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet'],
    'square_016um': [
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311/outs/binned_outputs/square_016um/spatial/tissue_positions.parquet', 
        '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/Loupe_Spaceranger/Spaceranger_311_no_json/outs/binned_outputs/square_016um/spatial/tissue_positions.parquet'],
}

for table_layer, paths in paths_dict.items():

    sdata = add_spaceranger_coordinates(
        sdata, 
        table_layer = table_layer, 
        new_tissue_positions_path = paths[0],
        spatial_key = 'spatial_spaceranger_311_json'
    )
    
    sdata = add_spaceranger_coordinates(
        sdata, 
        table_layer = table_layer, 
        new_tissue_positions_path = paths[1],
        spatial_key = 'spatial_spaceranger_311_no_json'
    )

### Calculate distances between different registration attempts
    
def calculate_distances(sdata, table_layer, spatial_key_1, spatial_key_2, column_tag, fullres_microns_per_pixel):
    gc.collect()
    
    # Copy adata
    adata = sdata.tables[table_layer].copy()
    
    # Get coordinates in micron
    coords_1 = adata.obsm[spatial_key_1] * fullres_microns_per_pixel
    coords_2 = adata.obsm[spatial_key_2] * fullres_microns_per_pixel
    
    # Calculate distance between coordinates
    dx = coords_1[:, 0] - coords_2[:, 0]
    dy = coords_1[:, 1] - coords_2[:, 1]
    distances = np.sqrt(dx**2 + dy**2)
    
    adata.obs[f'{column_tag}_distance'] = distances
    adata.obs[f'{column_tag}_dx'] = dx
    adata.obs[f'{column_tag}_dy'] = dy
    
    # Overwrite table layer
    sdata = add_table_layer(
        sdata,
        adata=adata,
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY][0],
        overwrite=True,
    )
    
    return sdata

for table_layer in ['square_002um', 'square_008um', 'square_016um']:
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_bigwarp_affine_JM',
        column_tag = '1',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_bigwarp_tps_JM',
        column_tag = '2',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_JM', 
        spatial_key_2 = 'spatial_bigwarp_tps_JM',
        column_tag = '3',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_spaceranger_311_json',
        column_tag = '4',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_spaceranger_311_no_json',
        column_tag = '5',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_spaceranger_311_json', 
        spatial_key_2 = 'spatial_spaceranger_311_no_json',
        column_tag = '6',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_JM', 
        spatial_key_2 = 'spatial_bigwarp_affine_KH',
        column_tag = '7',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_tps_JM', 
        spatial_key_2 = 'spatial_bigwarp_tps_KH',
        column_tag = '8',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_KH', 
        spatial_key_2 = 'spatial_bigwarp_tps_KH',
        column_tag = '9',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_bigwarp_affine_counts_JM',
        column_tag = '10',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial', 
        spatial_key_2 = 'spatial_bigwarp_tps_counts_JM',
        column_tag = '11',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_JM', 
        spatial_key_2 = 'spatial_bigwarp_affine_counts_JM',
        column_tag = '12',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_JM', 
        spatial_key_2 = 'spatial_bigwarp_tps_counts_JM',
        column_tag = '13',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_tps_JM', 
        spatial_key_2 = 'spatial_bigwarp_tps_counts_JM',
        column_tag = '14',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )
    
    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_tps_KH', 
        spatial_key_2 = 'spatial_bigwarp_tps_counts_JM',
        column_tag = '15',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_bigwarp_affine_counts_JM', 
        spatial_key_2 = 'spatial_bigwarp_tps_counts_JM',
        column_tag = '16',
        fullres_microns_per_pixel = 0.22009446422471257, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_cytassist', 
        spatial_key_2 = 'spatial_cytassist_311_json',
        column_tag = '17',
        fullres_microns_per_pixel = 4.625, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_cytassist', 
        spatial_key_2 = 'spatial_cytassist_311_no_json',
        column_tag = '18',
        fullres_microns_per_pixel = 4.625, 
    )

    sdata = calculate_distances(
        sdata, 
        table_layer = table_layer, 
        spatial_key_1 = 'spatial_cytassist_311_json', 
        spatial_key_2 = 'spatial_cytassist_311_no_json',
        column_tag = '19',
        fullres_microns_per_pixel = 4.625, 
    )
    
### Plot distances and displacement field
    
def violin_plot_distances(df, column, output_path):
    # Create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create violin plot
    df = df[df['in_tissue'] == 1]
    
    plt.figure(figsize=(4, 8))
    sns.violinplot(data=df, y=column)
    plt.ylim(0, 35) 
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_displacement_field(adata, spatial_key, column_tag, grid_size, figsize, landmark_points_list, output_path):
    # Create output directory
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Get coordinates
    df = adata.obs
    df['x'] = adata.obsm[spatial_key][:, 0]
    df['y'] = adata.obsm[spatial_key][:, 1]

    # subset the grid
    df = df[df['in_tissue'] == 1]
    df = df[(df['array_col'] % grid_size == 0) & (df['array_row'] % grid_size == 0)]
    
    # Create the plot
    plt.figure(figsize=figsize)
    
    # Plot the landmark points
    colors = ['r', 'k', 'b']
    for idx, landmark_points in enumerate(landmark_points_list):
        plt.scatter(landmark_points['x'], landmark_points['y'], s=50, c=colors[idx], marker='+')
        for i, row in landmark_points.iterrows():
            label = row['name'].replace('Pt-', '')
            plt.text(row['x'], row['y'] - 200, label, fontsize=13, ha='center', c=colors[idx])
            
    # Plot the displacement vectors
    max_magnitude = df[f'{column_tag}_distance'].max() # To normalize the arrows
    quiver = plt.quiver(df['x'], df['y'], df[f'{column_tag}_dx']/max_magnitude, df[f'{column_tag}_dy']/max_magnitude, df[f'{column_tag}_distance'], angles='xy', scale_units='xy', scale=0.001, width=0.001, headwidth=3, headlength=5, headaxislength=4.5, minshaft=1, minlength=1, alpha=1)
    
    # Add colorbar
    cbar = plt.colorbar(quiver, orientation='vertical')
    cbar.set_label('Distance (μm)', rotation=270, labelpad=15)
    
    # Set plot limits and labels
    plt.ylim(figsize[1], 0)
    plt.axis('equal')

    # Save the figure
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()

# Create new shapes layer for plotting in tissue
for table_layer in ['square_016um', 'square_008um', 'square_002um']:
    
    gdf = sdata.shapes[f"_{table_layer}"]
    cell_IDs = sdata.tables[table_layer].obs[sdata.tables[table_layer].obs['in_tissue'] == 1][_INSTANCE_KEY]
    filtered_gdf = gdf.loc[cell_IDs]
    
    transformations = get_transformation(sdata.shapes[f"_{table_layer}"])
    
    sdata = add_shapes_layer(
        sdata, 
        input = filtered_gdf, 
        output_layer = f"_{table_layer}_in_tissue",
        transformations = None,
        overwrite = True
    )
    
# Define plotting parameters
se=_get_spatial_element(sdata, layer="_full_image")

shape_x = se.shape[2]
shape_y = se.shape[1]

conversion_factor = 2000
figsize = (shape_x/conversion_factor,shape_y/conversion_factor)

grid_sizes = {'square_002um': 64, 'square_008um': 16, 'square_016um': 8}

# Create plots
for table_layer in ['square_002um', 'square_008um', 'square_016um']:
    grid_size = grid_sizes[table_layer]
    
    for column_tag in range(1,20):
        violin_plot_distances(
            df = sdata.tables[table_layer].obs, 
            column = f'{column_tag}_distance',
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/violin_plots/violin_plot_{column_tag}'
        )
        
        sp.pl.plot_shapes(
            sdata,
            shapes_layer = f"_{table_layer}_in_tissue",
            img_layer = "_full_image",
            channel = 0,
            linewidth = 0,
            alpha = 0.95,
            column = f'{column_tag}_distance',
            cmap = 'viridis',
            figsize = figsize,
            table_layer = table_layer,
            output = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/distance_plots/2_distance_plot_{column_tag}'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/0_displacement_field_plot_{column_tag}_no_landmarks'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [sdata.points['landmarks_JM']], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/1_displacement_field_plot_{column_tag}_landmarks_JM'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [sdata.points['landmarks_counts_JM']], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/2_displacement_field_plot_{column_tag}_landmarks_counts_JM'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [sdata.points['landmarks_KH']], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/3_displacement_field_plot_{column_tag}_landmarks_KH'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [sdata.points['landmarks_JM'], sdata.points['landmarks_KH']], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/4_displacement_field_plot_{column_tag}_landmarks_JM_and_KH'
        )
        
        plot_displacement_field(
            adata = sdata.tables[table_layer], 
            spatial_key = 'spatial',
            column_tag = column_tag,
            grid_size = grid_size, 
            figsize = figsize,
            landmark_points_list = [sdata.points['landmarks_JM'], sdata.points['landmarks_counts_JM']], 
            output_path = f'E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/plots/{table_layer}/displacement_field_plots/5_displacement_field_plot_{column_tag}_landmarks_JM_and_counts_JM'
        )
        
# Create adata
output_path = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/1/SCA002_reseq/registration_tests/harpy/adata'
os.makedirs(output_path, exist_ok=True)

for table_layer in sdata.tables:
    adata = sdata.tables[table_layer].copy()
    adata.write(f'{output_path}/{table_layer}.h5ad')