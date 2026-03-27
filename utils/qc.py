import scanpy as sc
import harpy as hp
from harpy.utils._keys import _INSTANCE_KEY, _REGION_KEY

def qc_metrics(sdata, table_layer, species):
    adata = sdata.tables[table_layer].copy()

    adata.var_names_make_unique()
    if species == 'human':
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        qc_vars = ["mt"]
    elif species == 'mouse':
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        qc_vars = ["mt"]
    elif species == 'rice':
        qc_vars = []
    else:
        raise ValueError("Species must be either 'human', 'mouse' or 'rice'.")

    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=[50], inplace=True)
    
    # Replace NaN with 0 for all columns starting with 'pct'
    for col in adata.obs.columns:
        if col.startswith("pct"):
            adata.obs[col] = adata.obs[col].fillna(0)
    
    # Add adata as table layer
    region = adata.obs[_REGION_KEY].cat.categories.to_list()
    
    sdata = hp.tb.add_table_layer(
        sdata, 
        adata=adata, 
        output_layer=table_layer, 
        region=region, 
        overwrite=True
    )

    return sdata