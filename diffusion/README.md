# Spatial accuracy analysis

## Setup environment

```bash
uv venv --python=3.12
source .venv/bin/activate
uv pip install spatialdata spatialdata_plot matplotlib numpy geopandas shapely ipykernel
```

## Regression testing

Run the regression checks from this directory:

```bash
.venv/bin/python -m unittest discover -s tests -v
```
