from __future__ import annotations

from pathlib import Path

import os
import geopandas as gpd
import pandas as pd
import numpy as np
import scanpy as sc
from dask_image.imread import imread
from spatialdata import SpatialData, read_zarr
from spatialdata.models import ShapesModel
from spatialdata.models import TableModel
from spatialdata.transformations import get_transformation, remove_transformation, set_transformation
# from spatialdata_io._constants._constants import VisiumKeys
# from spatialdata_io.readers.visium import visium as sdata_visium

import harpy as hp
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

# TODO: Add support for reading in fluorescence image, segmentation, ...
# TODO: Spot sizes are not the actual size they should be in micrometer
# TODO: Match barcode names for all levels

def bmkmanu(
    path: str | Path,
    dataset_id: str | None = None,
    filtered_counts: bool = False,
    add_in_tissue_column: bool = True,
    levels: list[int] | None = None,
    technology: str | None = None,
    bc_path: str | Path | None = None,
    output: str | Path | None = None,
) -> SpatialData:
    """
    Read *BMKGENE* BMKMANU S1000/S3000 formatted dataset.

    Parameters
    ----------
    path
        Path to directory containing the BMKMANU S1000/S3000 output.
    dataset_id
        Unique identifier of the dataset.
    filtered_counts
        If `True`, only the filtered counts (`{path}/BSTMAtrix/05.AllheStat/BSTViewer_project/subdata` directories) will be read. 
        If `False`, the raw counts (`{path}/BSTMAtrix/05.AllheStat/BSTViewer_project/level_matrix` directories) will be read.
    add_in_tissue_column
        If `True`, an `in_tissue` column will be added to the `obs` DataFrame of the table layers.
    levels
        List of levels to include in the dataset. If `None`, all levels will be included.
    image_path
        Path to the registered H&E image from BMChiper.
    technology
        Specifies whether BMKGENE S1000 or S3000 was used. Should be set to either 'S1000' or 'S3000'.
    bc_path
        Path to .bc file for level 1 (L1) outputted by BSTMAtrix level_bc2bc.pl script. This maps the level names from the data stored in 05.AllheStat/BSTViewer_project/subdata 
        and 05.AllheStat/BSTViewer_project/level_matrix to the barcode names in the config_prefix folder. This allows us to remove hypothetical spots for level 1 that are not actually 
        present on the capture slide. If None, hypothetical spots will be kept in.
    output
        The path where the resulting `SpatialData` object will be backed. If None, it will not be backed to a zarr store.
    
    """
    
    # Create SpatialData object from H&E image
    img_path = os.path.join(path, '05.AllheStat/BSTViewer_project/he.tif')

    sdata = hp.io.create_sdata(
        input=img_path, 
        img_layer=f'{dataset_id}_he', 
        dims=['c', 'y', 'x'],
        chunks=2048,
        to_coordinate_system=dataset_id,
        scale_factors=[2, 2, 2, 2]
    )
    
    # Get all available levels (= all available folders)
    if filtered_counts:
        counts_path = os.path.join(path, '05.AllheStat/BSTViewer_project/subdata')
    else:
        counts_path = os.path.join(path, '05.AllheStat/BSTViewer_project/level_matrix')
    
    if levels is not None and filtered_counts:
        level_names = [f"L{level}_heAuto" for level in levels]
    elif levels is not None and not filtered_counts:
        level_names = [f"level_{level}" for level in levels]
    else:
        level_names = [name for name in os.listdir(counts_path) if os.path.isdir(os.path.join(counts_path, name))]
    
    level_paths = {name: os.path.join(counts_path, name) for name in level_names}
    
    # Add levels to SpatialData object
    he_shape = hp.im._image._get_spatial_element(sdata, f'{dataset_id}_he').shape
    
    he_width = he_shape[2]
    he_height = he_shape[1]

    for level_name, level_path in level_paths.items():
        
        if filtered_counts:
            level = int(level_name.split('_')[0].replace('L', ''))
        else:
            level = int(level_name.split('_')[1])
        
        sdata = _add_BMK_level(
            sdata,
            dataset_id,
            level,
            level_name,
            level_path,
            he_width,
            he_height,
            technology,
            add_in_tissue_column,
            bc_path,
        )
        
    if output is not None:
        sdata.write(output)
        sdata = read_zarr(sdata.path)

    return sdata

def _calculate_scaling_factor(width, height, technology):
    std_width = 1000
    
    if technology == 'S1000':
        std_height = std_width / (46 * 31) * (46 * 36 * np.sqrt(3) / 2.0)
    elif technology == 'S3000':
        std_height = std_width / (42 * 46) * (43 * 52 * np.sqrt(3) / 2.0)
    else:
        raise ValueError('technology should be set to either `S1000` or `S3000`')
    
    if std_width / std_height > width / height:
        scale = width / std_width
    else:
        scale = height / std_height
        
    return scale
    
def _calculate_radius(coordinates, width, height, subset_size=1000):
    x_center = width/2
    y_center = height/2

    x_min = x_center - (subset_size / 2)
    x_max = x_center + (subset_size / 2)
    y_min = y_center - (subset_size / 2)
    y_max = y_center + (subset_size / 2)

    subset_coords = coordinates[
        (coordinates[:, 0] >= x_min) & (coordinates[:, 0] <= x_max) &
        (coordinates[:, 1] >= y_min) & (coordinates[:, 1] <= y_max)
    ]
    
    distances = np.linalg.norm(subset_coords[:, np.newaxis] - subset_coords[np.newaxis, :], axis=-1)
    
    np.fill_diagonal(distances, np.inf) 
    min_distance = np.min(distances)

    radius = (min_distance / 2) * 0.98
    return radius

def _add_BMK_level(sdata, dataset_id, level, level_name, level_path, width, height, technology, add_in_tissue_column, bc_path):
    # Calculate scaling factor
    scaling_factor = _calculate_scaling_factor(width, height, technology)

    # Read in expression data and coordinates as adata
    adata = sc.read_10x_mtx(level_path, var_names='gene_symbols', cache=True)
    barcode_pos = pd.read_csv(f'{level_path}/barcodes_pos.tsv.gz', compression='gzip', sep='\t', names=["barcode", "array_col", "array_row"], header=None).set_index('barcode')
    adata.obs = adata.obs.join(barcode_pos, how="left")
    adata.obs['barcode'] = adata.obs.index
    
    coordinates = barcode_pos[['array_col','array_row']].to_numpy() * scaling_factor
    adata.obsm["spatial"] = coordinates

    # Add in_tissue column (if applicable)
    if add_in_tissue_column:
        tissue_barcodes_path = os.path.join(os.path.dirname(os.path.dirname(level_path)), f'subdata/L{level}_heAuto/barcodes_pos.tsv.gz')
        tissue_barcodes = pd.read_csv(tissue_barcodes_path, compression='gzip', sep='\t', names=["barcode", "array_col", "array_row"], header=None)["barcode"].values.tolist()
        adata.obs['in_tissue'] = adata.obs.index.isin(tissue_barcodes).astype(int)
        
    # Add read data
    barcodes_read_path = os.path.join(os.path.dirname(os.path.dirname(level_path)), f'level_matrix/level_{level}/barcodes_read.tsv.gz')
    barcodes_read = pd.read_csv(barcodes_read_path, compression='gzip', sep='\t', names=["barcode", "read_counts"], header=None).set_index("barcode")
    adata.obs = adata.obs.join(barcodes_read, how="left")
    adata.obs["read_counts"] = adata.obs["read_counts"].fillna(0).astype(int)
    
    # Remove hypothetical spots for level 1, if bc_path is given
    if level == 1 and bc_path is not None:
        bc = pd.read_csv(bc_path, sep='\t', names=['barcode_name', 'barcode'], header=None)

        adata.obs.index = adata.obs.index.astype(str)
        bc['barcode'] = bc['barcode'].astype(str)

        valid_barcodes = adata.obs.index.intersection(bc['barcode'])
        adata = adata[valid_barcodes].copy()

        bc = bc.set_index('barcode')
        adata.obs['barcode_name'] = adata.obs.index.map(bc['barcode_name'])

        
    # Add region key and change indexes
    adata.obs[_REGION_KEY] = f"{dataset_id}_{level_name}_labels"
    adata.obs[_INSTANCE_KEY] = [int(x.split('_')[1]) + 1 for x in adata.obs.index] # The + 1 is to ensure that the index starts from 1
    adata.obs.index = adata.obs[_INSTANCE_KEY]
    adata.obs.index.name = 'cells'
    
    # Add shapes and labels layer based on coordinates
    coordinates = adata.obsm["spatial"]
    radius = _calculate_radius(coordinates, width, height)
    radii = np.full(len(coordinates), radius)
    spots = ShapesModel.parse(coordinates, geometry=0, radius=radii)
    spots.index = adata.obs.index
    polygons = spots.buffer(radius, cap_style="round")
    
    sdata = hp.sh.add_shapes_layer(
        sdata, 
        gpd.GeoDataFrame(geometry=polygons), 
        output_layer=f"{dataset_id}_{level_name}", 
        transformations=get_transformation(sdata.images[f"{dataset_id}_he"], get_all=True),
        overwrite=True
    )
    
    sdata = hp.im.rasterize(
        sdata,
        shapes_layer=f"{dataset_id}_{level_name}", 
        output_layer=f"{dataset_id}_{level_name}_labels",
        chunks=2048,
        overwrite=True,
    )
    
    # Add adata as table layer
    sdata = hp.tb.add_table_layer(
        sdata, 
        adata=adata.copy(), 
        output_layer=level_name, 
        region=[F"{dataset_id}_{level_name}_labels"], 
        overwrite=True
    )

    return sdata