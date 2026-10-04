# Diffusion analysis

## Setup environment

```bash
uv venv --python=3.12
source .venv/bin/activate
uv pip install spatialdata spatialdata_plot matplotlib numpy geopandas shapely ipykernel
```

## Multiple boundary sections

Each `Experiment.boundary_point_pairs` entry is computed separately and saved in
`plots/<technology>/<dataset_id>/analysis/crop{i}/diffusion.npz`. Choose sections
far enough apart that their eligible analysis bins are disjoint at the configured
maximum diffusion distance. Pooling assumes this and performs no overlap detection
or deduplication.

All configured sections' signed distances and transcript counts are concatenated
per `(technology, dataset_id)` before normalization, KDE peak detection/alignment,
or summary statistics. Transcript counts remain the weights. Comparison plots in
`plots/combined/` show one distribution or summary row per experiment.

Spatial inspection plots stay in the crop folders. When enabled, individual
histograms and boxplots use pooled distances and are saved in `analysis/pooled/`.
Old distribution plots in crop folders are not updated. No additional pooled data
file is needed, and existing single-crop NPZ files remain usable.

After adding, changing, or reordering boundary pairs, enable `COMPUTE_DIFFUSION`
in the notebook to regenerate their results. Only explicitly configured crops
enter the summaries; leftover crop folders are ignored.

Run the pooling regression checks from this directory:

```bash
.venv/bin/python -m unittest discover -s tests -v
```
