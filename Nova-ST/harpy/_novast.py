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
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

def novast(
    path: str | Path,
    dataset_id: str,
    seq1st_id: str,
    chip_id: str,
    path_img: str | Path,
    output: str | Path | None = None,
) -> SpatialData:
    """
    Read *Nova-ST* formatted dataset.
    
    Parameters
    ----------
    path
        Path to directory containing NovaScope output.
    dataset_id
        Unique identifier of the dataset. This should match the `run_id` from the NovaScope config file.
    seq1st_id
        Identifier of the seq1st directories. This should match the `id` under `seq1st` in the NovaScope config file.
    chip_id
        Identifier for the Nova-ST chip. This should match the `chip` from the NovaScope config file.
    path_img
        Path to registered image.
    output
        The path where the resulting `SpatialData` object will be backed. If None, it will not be backed to a zarr store.
    """

    # Create sdata
    sdata = hp.io.create_sdata(
        input=path_img,
        img_layer=f'{dataset_id}_image',
        dims=['c', 'y', 'x'],
        chunks=1024,
        to_coordinate_system=dataset_id,
        scale_factors=(2, 2, 2, 2),
    )

    # Read barcodes
    barcodes = pd.read_csv(
        os.path.join(path, f"align/{seq1st_id}/{chip_id}/{dataset_id}/sge/barcodes.tsv.gz"),
        sep="\t", header=None,
        names=[
            "barcode", "barcode_index", "full_barcode_index",
            "lane_id", "tile_id", "x", "y", "counts"
        ]
    )
    count_split = barcodes["counts"].str.split(",", expand=True).astype(int)
    count_split.columns = ["Gene", "GeneFull", "Spliced", "Unspliced", "Ambiguous"]
    barcodes = pd.concat([barcodes.drop(columns=["counts"]), count_split], axis=1)

    barcodes = barcodes.set_index("barcode")
    barcodes.index.names = [None]

    barcodes = barcodes.drop(columns=["lane_id", "tile_id"]) # Drop lane_id and tile_id columns since they are just set to 1 and do not contain any useful information

    # Read features
    features = pd.read_csv(
        os.path.join(path, f"align/{seq1st_id}/{chip_id}/{dataset_id}/sge/features.tsv.gz"),
        sep="\t", header=None,
        names=["ensembl_id", "gene_symbol", "feature_index", "counts"]
    )

    feat_split = features["counts"].str.split(",", expand=True).astype(int)
    feat_split.columns = ["Gene", "GeneFull", "Spliced", "Unspliced", "Ambiguous"]

    features = pd.concat([features.drop(columns=["counts"]), feat_split], axis=1)
    features = features.set_index("gene_symbol")
    features.index.names = [None]

    # Read matrices
    rows = []
    cols = []
    values = {k: [] for k in ["Gene", "GeneFull", "Spliced", "Unspliced", "Ambiguous"]}

    with gzip.open(os.path.join(path, f"align/{seq1st_id}/{chip_id}/{dataset_id}/sge/matrix.mtx.gz"), "rt") as f:
        for line in f:
            # Lines starting with '%' are comments
            if line.startswith("%"):
                continue
            
            parts = line.strip().split()

            # First non-comment line is the size line
            if len(parts) == 3:
                n_rows, n_cols, nonzero = map(int, parts)
                continue
            
            r = int(parts[0]) - 1
            c = int(parts[1]) - 1
            vals = list(map(int, parts[2:]))

            rows.append(r)
            cols.append(c)
            for k, v in zip(values.keys(), vals):
                values[k].append(v)

    layers = {}
    for k, vals in values.items():
        layers[k] = sparse.coo_matrix(
            (vals, (cols, rows)),
            shape=(n_cols, n_rows),
        ).tocsr()

    # Create AnnData with Gene as main matrix
    adata = ad.AnnData(
        X=layers["Gene"],
        obs=barcodes,
        var=features
    )

    # Add the other matrices as layers
    for k in ["GeneFull", "Spliced", "Unspliced", "Ambiguous"]:
        adata.layers[k] = layers[k]

    # Gene: Counts reads mapping entirely to exonic regions of annotated transcripts, capturing fully spliced transcripts. This procedure represents traditional gene expression analyses and was the default in earlier CellRanger versions.
    # GeneFull: Counts reads overlapping entire gene locus, i.e. including exonic and intronic regions. This captures both unspliced (primary) and spliced transcripts and providing a more comprehensive view of gene activity.
    # Velocyto: Generates separate count matrices for spliced, unspliced, and ambiguous reads, following the rules from (La Manno et al., 2018). This enables RNA velocity analyses to infer dynamic cellular processes.
        
    # Remove indices
    adata.obs = adata.obs.drop(columns=["barcode_index", "full_barcode_index"])
    adata.var = adata.var.drop(columns=["feature_index"])

    # Add coordinates to obsm
    coordinates = barcodes[["x","y"]].to_numpy() / 1000 # Divide by 1000 to get coordinates in micrometer
    adata.obsm["spatial"] = coordinates
    adata.obs = adata.obs.drop(columns=["x", "y"])

    # Add region key and change indexes
    adata.obs[_REGION_KEY] = f"{dataset_id}_labels"
    adata.obs[_INSTANCE_KEY] = range(1, len(adata.obs) + 1)

    # Add shapes and labels layer based on coordinates
    coordinates = adata.obsm["spatial"]
    radius = 0.15
    radii = np.full(len(coordinates), radius)
    spots = ShapesModel.parse(coordinates, geometry=0, radius=radii)
    spots.index = adata.obs[_INSTANCE_KEY]
    polygons = spots.buffer(radius, cap_style="round")

    sdata = hp.sh.add_shapes_layer(
        sdata, 
        gpd.GeoDataFrame(geometry=polygons), 
        output_layer=f"{dataset_id}_shapes", 
        transformations=get_transformation(sdata.images[f"{dataset_id}_he"], get_all=True),
        overwrite=True
    )

    sdata = hp.im.rasterize(
        sdata,
        shapes_layer=f"{dataset_id}_shapes", 
        output_layer=f"{dataset_id}_labels",
        chunks=2048,
        overwrite=True,
        scale_factors=(2, 2, 2, 2),
    )

    # Add adata as table layer
    sdata = hp.tb.add_table_layer(
        sdata, 
        adata=adata.copy(), 
        output_layer=f"{dataset_id}_table", 
        region=[f"{dataset_id}_labels"], 
        overwrite=True
    )

    if output is not None:
        sdata.write(output)
        sdata = read_zarr(sdata.path)

    return sdata