import harpy as hp
import pandas as pd
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

def assign_region(centroid, mask, mask_shape):
    x, y = int(round(centroid[0])), int(round(centroid[1]))
    
    if 0 <= y < mask_shape[0] and 0 <= x < mask_shape[1]:
        return 1 if mask[y, x] == 1 else 0
    return 0  # Outside mask bounds

def add_in_tissue(sdata, table_layer, tissue_mask):
    adata = sdata.tables[table_layer].copy()
    coordinates = adata.obsm['spatial']
    coordinates_df = pd.DataFrame(coordinates, columns=["x", "y"], index=adata.obs.index)

    mask_shape = tissue_mask.shape
    adata.obs["in_tissue_annotation"] = coordinates_df.apply(
        lambda row: assign_region((row["x"], row["y"]), tissue_mask, mask_shape), axis=1
    )

    sdata = hp.tb.add_table_layer(
        sdata,
        adata=adata.copy(),
        output_layer=table_layer,
        region=adata.obs[_REGION_KEY].cat.categories.to_list(),
        overwrite=True,
    )
    
    return sdata