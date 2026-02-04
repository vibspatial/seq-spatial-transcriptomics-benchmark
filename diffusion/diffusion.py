from dataclasses import dataclass


@dataclass(frozen=True)
class LevelInfo:
    name: str 
    table_name: str  
    labels_layer: str 


@dataclass(frozen=True)
class Experiment:
    technology: str  # spatial omics technology used to analyze sample
    species: str  # species of sample
    tissue: str  # tissue type of sample
    dataset_id: str  # shorthand ID of dataset
    sdata_path: str  # path to SpatialData zarr with analysis output
    he_image: str | None  # name of H&E image layer
    tissue_mask: str | None  # name of tissue mask shapes layer
    levels: list[LevelInfo]  # bin sizes, they MUST be ordered from finest to coarsest
    crop_rects: list[tuple[int, int, int, int]]  # list of crop rectangles (xmin, ymin, xmax, ymax) in pixel coordinates
