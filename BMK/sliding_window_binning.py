# Environment: harpy

# Import packages
print('*** importing packages ***', flush=True)
import sparrow as sp
import os
import numpy as np
import spatialdata as sd
from sparrow.table import add_table_layer
from skimage import io, img_as_ubyte
from scipy.sparse import csr_matrix
from sklearn.neighbors import BallTree
from joblib import Parallel, delayed

# Monitor memory usage
import psutil
def print_memory_usage():
    process = psutil.Process(os.getpid())
    print(f'Memory usage: {process.memory_info().rss / 1024 ** 3:.2f} GB', flush=True)
print_memory_usage()

# Specify paths
print('*** specifying paths ***', flush=True)
sdata_path = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/BMK/exp74/harpy/sdata.zarr'
output_path = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/BMK/exp74/harpy'
temp_folder = '/data/projects/s00/2024-07-UTBenchmark-SpC/data/processed/BMK/exp74/harpy/temp'

# Read sdata
print('*** reading sdata ***', flush=True)
sdata = sd.read_zarr(sdata_path)
print_memory_usage()

# Define hires_img
print('*** defining hires_img ***', flush=True)
hires_img = sdata['raw_image'].compute().data

# Define sliding window function
def sliding_window_binning(adata, radius):

    # Get coordinates and expression data
    coordinates = adata.obsm['spatial']
    expression_data = adata.X 
    if not isinstance(expression_data, csr_matrix):
        expression_data = csr_matrix(expression_data)

    # Calculate pairwise distances
    tree = BallTree(coordinates)

    # Function to aggregate expression data for a single spot
    def aggregate_spot(i):
        # Query the ball tree to find all spots within the radius
        nearby_spots = tree.query_radius([coordinates[i]], r=radius)[0]

        # Sum the expression data of nearby spots
        aggregated_expression = expression_data[nearby_spots].sum(axis=0)
        return aggregated_expression

    # Aggregate expression
    aggregated_expression = Parallel(n_jobs=-3, temp_folder=temp_folder)(
        delayed(aggregate_spot)(i) for i in range(coordinates.shape[0])
    )

    print('aggregate expression DONE!', flush=True)
    print_memory_usage()

    # Convert the list of aggregated sparse arrays into a sparse matrix
    aggregated_expression = csr_matrix(np.vstack(aggregated_expression))
    print('vstack DONE!', flush=True)
    print_memory_usage()

    # Create a new AnnData object and assign the aggregated expression data
    new_adata = adata.copy()
    new_adata.X = aggregated_expression
    print('new_adata DONE!', flush=True)
    print_memory_usage()
    
    return new_adata

# Calculate sliding window expression values
print('*** calculating sliding window expression values ***', flush=True)
level_name = 'L1_heAuto'
new_adata = sliding_window_binning(sdata.tables[f'table_{level_name}'], radius=20)
sdata = add_table_layer(sdata, adata=new_adata.copy(), output_layer=f'table_{level_name}_sliding_window', region=[f"labels_spots_{level_name}"], overwrite=True)

# Leiden clustering
print('*** leiden clustering ***', flush=True)
import gc # garbage collection to avoid out-of-memory issues.

for table_name in [level_name, f'{level_name}_sliding_window']:
    print(f'** running {table_name} **', flush=True)
    # Preprocess data
    print('** preprocessing **', flush=True)
    sdata = sp.tb.preprocess_transcriptomics(
        sdata=sdata,
        size_norm=False, # scanpy.sc.pp.normalize_total is used for normalization
        labels_layer=f"labels_spots_{level_name}",
        table_layer=f"table_{table_name}",
        output_layer=f"table_{table_name}_processed",
        min_counts=50,
        min_cells=5,
        highly_variable_genes=False,
        max_value_scale=10,
        n_comps = 50,
        update_shapes_layers=False,
        overwrite=True
    )

    sp.pl.preprocess_transcriptomics(
        sdata, 
        table_layer=f"table_{table_name}_processed",
        output=os.path.join(output_path, f"{table_name}_preprocess")
    )
    print_memory_usage()
    gc.collect()

    # Leiden clustering
    leiden_resolutions = [0.6]

    for res in leiden_resolutions:
        res_str = str(res).replace('.', '')
        print(f'** clustering at resolution {res} **', flush=True)
        sdata = sp.tb.leiden(
            sdata, 
            labels_layer=f"labels_spots_{level_name}",
            table_layer=f"table_{table_name}_processed",
            output_layer=f"table_{table_name}_processed",
            calculate_umap=True,
            calculate_neighbors=True,
            rank_genes=True,
            n_neighbors=35,
            n_pcs=30,
            resolution=res,
            key_added=f'leiden_{res_str}',
            index_names_var=None,
            index_positions_var=None,
            random_state=100,
            overwrite=True
        )

        sp.pl.cluster(
            sdata, 
            table_layer=f"table_{table_name}_processed",
            key_added=f'leiden_{res_str}',
            output=os.path.join(output_path, f'UMAP_{table_name}_{res_str}')
        )

    print_memory_usage()
    gc.collect()

    # Save raw adata with H&E added
    print('** saving raw adata **', flush=True)
    adata = sdata.tables[f'table_{table_name}'].copy()

    library_id = level_name
    adata.uns["spatial"] = {library_id: {}}
    adata.uns["spatial"][library_id]["images"] = {}
    adata.uns['spatial'][library_id]['images']['hires'] = img_as_ubyte(hires_img.transpose(1, 2, 0))
    adata.uns['spatial'][library_id]['use_quality'] = 'hires' 
    adata.uns['spatial'][library_id]['scalefactors'] = {"tissue_hires_scalef": 1}

    adata.write(f'{output_path}/adata_{table_name}.h5ad')
    
    # Save processed adata with H&E added
    print('** saving processed adata **', flush=True)
    adata = sdata.tables[f'table_{table_name}_processed'].copy()

    library_id = level_name
    adata.uns["spatial"] = {library_id: {}}
    adata.uns["spatial"][library_id]["images"] = {}
    adata.uns['spatial'][library_id]['images']['hires'] = img_as_ubyte(hires_img.transpose(1, 2, 0))
    adata.uns['spatial'][library_id]['use_quality'] = 'hires' 
    adata.uns['spatial'][library_id]['scalefactors'] = {"tissue_hires_scalef": 1}

    adata.write(f'{output_path}/adata_{table_name}_processed.h5ad')

print("*** DONE! ***", flush=True)
