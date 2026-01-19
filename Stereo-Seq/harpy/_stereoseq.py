from __future__ import annotations

from pathlib import Path

import os
import geopandas as gpd
import pandas as pd
import numpy as np
from spatialdata import SpatialData, read_zarr
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation, remove_transformation, set_transformation

import harpy as hp
import gzip
from scipy import sparse
import anndata as ad
import dask.array as da
from dask_image import imread
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

from spatialdata import rasterize_bins, rasterize_bins_link_table_to_labels
from spatialdata.models import TableModel


def stereoseq(
    regist_img_path: str | Path,
    dataset_id: str | None,
    raw_h5ad_path: str | Path | None = None,
    tissue_h5ad_path: str | Path | list[str | Path] | None = None,
    tissue_cut_path: str | Path | None = None,
    output: str | Path | None = None,
) -> SpatialData:
    """
    Read *Stereo=Seq* formatted dataset.
    
    Parameters
    ----------
    regist_img_path
        Path to registered image created by the SAW count pipeline (`<SN>_<stainType>_regist.tif`).
    dataset_id
        Unique identifier of the dataset.
    raw_h5ad_path
        Path to .h5ad output created from `<SN>.raw.gef` by `saw convert gef2h5ad`.
    tissue_h5ad_path
        Path or list of paths to .h5ad output created from `<SN>.tissue.gef` by `saw convert gef2h5ad`. 
        If multiple paths are provided, these should correspond to .h5ad files with different bin sizes.
    tissue_cut_path
        Path to the tissue segmentation mask created by the SAW count pipeline (`<SN>_<stainType>_tissue_cut.tif`).
    output
        The path where the resulting `SpatialData` object will be backed. If None, it will not be backed to a zarr store.
    """

    # Create sdata
    sdata = hp.io.create_sdata(
        input=regist_img_path,
        img_layer=f"{dataset_id}_image",
        dims=["c", "y", "x"],
        chunks=1024,
        to_coordinate_system=dataset_id,
        scale_factors=(2, 2, 2, 2),
    )

    # Add SAW tissue segmentation mask
    if tissue_cut_path is not None:
        dask_array = imread.imread(tissue_cut_path).squeeze()
        sdata = hp.im.add_labels_layer(
            sdata,
            arr=dask_array,
            output_layer=f"{dataset_id}_tissue_cut",
            chunks=1024,
            transformations=get_transformation(sdata.images[f"{dataset_id}_image"], get_all=True),
            scale_factors=(2, 2, 2, 2),
        )

    # Process h5ad files
    if raw_h5ad_path is not None:
        adata = ad.read_h5ad(raw_h5ad_path)
        if "bin_size" not in adata.uns:
            raise ValueError(f"{path} does not contain adata.uns['bin_size']")
        if "spatial" not in adata.obsm:
            raise ValueError(f"{path} does not contain adata.obsm['spatial']")

        prefix = f"{dataset_id}_raw"
        sdata = _process_SAW_h5ad(sdata, dataset_id, adata, prefix, adata.uns["bin_size"])

    if tissue_h5ad_path is not None:
        if not isinstance(tissue_h5ad_path, list):
            tissue_h5ad_path = [tissue_h5ad_path]

        h5ad_files = []
        for path in tissue_h5ad_path:
            adata = ad.read_h5ad(path)
            if "bin_size" not in adata.uns:
                raise ValueError(f"{path} does not contain adata.uns['bin_size']")
            if "spatial" not in adata.obsm:
                raise ValueError(f"{path} does not contain adata.obsm['spatial']")
            h5ad_files.append((adata, adata.uns["bin_size"]))

        bin_sizes = [b for _, b in h5ad_files]
        if len(bin_sizes) != len(set(bin_sizes)):
            raise ValueError(
                "Multiple tissue h5ad files share the same bin_size. "
                "Please check your input files."
            )

        for adata, bin_size in h5ad_files:
            prefix = f"{dataset_id}_bin_{bin_size}"
            sdata = _process_SAW_h5ad(sdata, dataset_id, adata, prefix, bin_size)

    # Save output
    if output is not None:
        sdata.write(output)
        sdata = read_zarr(sdata.path)

    return sdata


def _process_SAW_h5ad(sdata, dataset_id, adata, prefix, bin_size):
    # Add Anndata as table
    adata.X = adata.X.tocsc()

    adata.obs["temporary_region_key"] = f"{prefix}_shapes" # Temporarily set to shapes layer for rasterize_bins
    adata.obs["temporary_instance_key"] = range(len(adata.obs)) # Temporarily needs to be 0-based to get correctly handled by rasterize_bins
    adata.obs["array_col"] = (adata.obs["x"] / bin_size).astype(int)
    adata.obs["array_row"] = (adata.obs["y"] / bin_size).astype(int)

    sdata[f"{prefix}_table"] = TableModel.parse(
        adata,
        region=[f"{prefix}_shapes"],
        region_key="temporary_region_key",
        instance_key="temporary_instance_key",
    )

    # Create square polygons based on coordinates
    coordinates = adata.obsm["spatial"]
    radius = bin_size / 2
    radii = np.full(len(coordinates), radius)
    spots = ShapesModel.parse(coordinates, geometry=0, radius=radii)
    spots.index = adata.obs["temporary_instance_key"]
    polygons = spots.buffer(radius, cap_style=3)

    # Add shapes layer
    sdata = hp.sh.add_shapes_layer(
        sdata,
        gpd.GeoDataFrame(geometry=polygons),
        output_layer=f"{prefix}_shapes",
        transformations=get_transformation(sdata.images[f"{dataset_id}_image"], get_all=True),
        overwrite=True,
    )

    # add labels layer (rasterized bins).
    labels_element = rasterize_bins(
        sdata,
        bins=f"{prefix}_shapes",
        table_name=f"{prefix}_table",
        row_key="array_row",
        col_key="array_col",
        value_key=None,
        return_region_as_labels=True,
    )

    sdata[f"{prefix}_labels"] = labels_element
    rasterize_bins_link_table_to_labels(
        sdata=sdata,
        table_name=f"{prefix}_table",
        rasterized_labels_name=f"{prefix}_labels",
    )

    # Correct region and instance key
    adata = sdata[f"{prefix}_table"].copy()
    _relabeled_instance_key = adata.uns[TableModel.ATTRS_KEY][TableModel.INSTANCE_KEY]

    adata.obs.rename(columns={"temporary_region_key": _REGION_KEY, _relabeled_instance_key: _INSTANCE_KEY}, inplace=True)
    adata.uns.pop(TableModel.ATTRS_KEY)

    adata = TableModel.parse(
        adata,
        region_key=_REGION_KEY,
        region=adata.obs[_REGION_KEY].cat.categories.to_list(),
        instance_key=_INSTANCE_KEY,
    )

    sdata[f"{prefix}_shapes"].index = (
        adata.obs.set_index("temporary_instance_key").loc[sdata[f"{prefix}_shapes"].index, _INSTANCE_KEY].values
    )
    sdata[f"{prefix}_shapes"].index.name = _INSTANCE_KEY

    for col in ["temporary_region_key", "temporary_instance_key"]:
        if col in adata.obs.columns:
            adata.obs.drop(columns=col, inplace=True)

    del sdata[f"{prefix}_table"]
    sdata[f"{prefix}_table"] = adata

    return sdata