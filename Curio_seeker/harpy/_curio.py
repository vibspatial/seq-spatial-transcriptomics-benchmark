from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
from spatialdata import SpatialData, read_zarr
from spatialdata.models import TableModel
from spatialdata.transformations import get_transformation, set_transformation, remove_transformation
from spatialdata_curio import curio as sdata_curio

import harpy as hp
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

from shapely.affinity import translate as shapely_translate
import numpy as np

def curio(
    path: str | Path,
    dataset_id: str | None = None,
    output: str | Path | None = None,
) -> SpatialData:
    """
    Read *Curio Seeker* formatted dataset.

    Wrapper around `spatialdata.io.readers.curio.curio`, but with the resulting table annotated by a labels layer.

    Parameters
    ----------
    path
        Path to directory containing the *10x Genomics* Visium output.
    dataset_id
        Unique identifier of the dataset.
    output
        The path where the resulting `SpatialData` object will be backed. If None, it will not be backed to a zarr store.
    """
    sdata = sdata_curio(
        path=path,
    )

    for table_layer in [*sdata.tables]:
        adata = sdata[table_layer].copy()
        adata.var_names_make_unique()
        adata.X = adata.X.tocsc()

        _old_instance_key = sdata[table_layer].uns[TableModel.ATTRS_KEY][TableModel.INSTANCE_KEY]
        _old_region_key = sdata[table_layer].uns[TableModel.ATTRS_KEY][TableModel.REGION_KEY_KEY]
        adata.obs[_REGION_KEY] = pd.Categorical([f"{dataset_id}_labels"] * len(adata.obs), categories=[f"{dataset_id}_labels"])
        adata.obs[_INSTANCE_KEY] = range(1, len(adata.obs) + 1)
        adata.obs['barcode'] = adata.obs[_old_instance_key]

        adata.uns.pop(TableModel.ATTRS_KEY)
        adata = TableModel.parse(
            adata,
            region_key=_REGION_KEY,
            region=adata.obs[_REGION_KEY].cat.categories.to_list(),
            instance_key=_INSTANCE_KEY,
        )
        
        sdata.shapes[dataset_id] = sdata.shapes['cells']
        del sdata.shapes['cells']

        assert len(sdata.shapes[dataset_id]) == len(
            adata
        ), f"Shapes layer '{dataset_id}' and corresponding table '{table_layer}' should have same length."

        sdata.shapes[dataset_id].index = (
            adata.obs.set_index(_old_instance_key).loc[sdata.shapes[dataset_id].index, _INSTANCE_KEY].values
        )
        sdata.shapes[dataset_id].index.name = _INSTANCE_KEY

        if _old_region_key in adata.obs.columns:
            adata.obs.drop(columns=_old_region_key, inplace=True)

        if _old_instance_key in adata.obs.columns:
            adata.obs.drop(columns=_old_instance_key, inplace=True)

        # Convert Points to Polygons and translate so all values are positive
        radius = _calculate_radius(coordinates=adata.obsm["spatial"])
        polygons = sdata.shapes[dataset_id].buffer(radius, cap_style="round")

        min_x, min_y, max_x, max_y, = polygons.total_bounds
        offset_x, offset_y = -min_x, -min_y
        
        polygons = polygons.apply(
            lambda geom: shapely_translate(geom, xoff=offset_x, yoff=offset_y)
        )

        adata.obsm["spatial"] = adata.obsm["spatial"] + np.array([offset_x, offset_y])

        sdata = hp.sh.add_shapes_layer(
            sdata, 
            gpd.GeoDataFrame(geometry=polygons), 
            output_layer=dataset_id, 
            overwrite=True
        )
        
        transformations = get_transformation(sdata.shapes[dataset_id], get_all=True)
        set_transformation(
            sdata.shapes[dataset_id], transformations["global"], to_coordinate_system=dataset_id
        )
        remove_transformation(sdata.shapes[dataset_id], to_coordinate_system="global")

        # Create labels layer
        sdata = hp.im.rasterize(
            sdata,
            shapes_layer=dataset_id,
            output_layer=f"{dataset_id}_labels",
            chunks=5000,
            overwrite=True,
        )

        del sdata[table_layer]

        sdata[table_layer] = adata

    if output is not None:
        sdata.write(output)
        sdata = read_zarr(sdata.path)

    return sdata

def _calculate_radius(coordinates, subset_size=1000):
    x_center, y_center = np.median(coordinates, axis=0)

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