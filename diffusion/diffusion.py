from dataclasses import dataclass, field


@dataclass(frozen=True)
class LevelInfo:
    name: str 
    table_name: str  
    labels_layer: str 


@dataclass(frozen=True)
class Experiment:
    technology: str  # spatial omics technology used to analyze sample
    tissue: str  # tissue type of sample
    dataset_id: str  # shorthand ID of dataset
    sdata_path: str  # relative path to (subsetted) SpatialData zarr with analysis output
    he_image: str | None  # name of H&E image layer
    pixel_size_um: float # pixel size of H&E image in microns; units of coordinate system used for e.g. spot centroids
    tissue_mask: str | None  # name of tissue mask shapes layer
    levels: list[LevelInfo]  # bin sizes, they MUST be ordered from finest to coarsest
    boundary_point_pairs: list[tuple[tuple[float, float], tuple[float, float]]] = field(default_factory=list)  # list of (approx_start_xy, approx_end_xy) pairs near tissue boundary
