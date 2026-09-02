from __future__ import annotations

import json
import math
from pathlib import Path
from time import perf_counter

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection, LineCollection, PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, QhullError
from scipy.spatial import cKDTree
import shapely
from shapely import affinity
from shapely.errors import GEOSException
from shapely.geometry import Point, box
from shapely.ops import nearest_points, unary_union

import logging
logging.getLogger('ome_zarr').setLevel(logging.ERROR) # This mutes unnecessary Zarr warning (see: https://github.com/scverse/spatialdata/issues/400)


def _elapsed(start_time, label):
    """Print and return elapsed wall time for a pipeline stage."""
    seconds = perf_counter() - start_time
    print(f"  [timing] {label}: {seconds:.2f} s")
    return seconds


# size_um corresponds to the diameter for circles and the side length for squares.
SPOT_CONFIG = {
    "Nova-ST": {"shape": "circle", "size_um": 0.30},
    "Stereo-Seq": {"shape": "circle", "size_um": 0.22},
    "Curio Seeker": {"shape": "circle", "size_um": 10.0},
    "Visium (FFPE)": {"shape": "circle", "size_um": 55.0},
    "Visium HD (FFPE)": {"shape": "square", "size_um": 2.0},
    "Visium HD (FF)": {"shape": "square", "size_um": 2.0},
    "Visium HD 3'": {"shape": "square", "size_um": 2.0},
    "BMKMANU S1000": {"shape": "circle", "size_um": 2.5},
    "BMKMANU S3000": {"shape": "circle", "size_um": 2.5},
}

CURIO_SOURCE_SPOT_RADIUS_PX = 10.0

PLOT_METRICS = {
    "n_neighbors": "Delaunay neighbors",
    "mean_neighbor_distance_um": "Mean neighbor distance (µm)",
    "voronoi_area_um2": "Voronoi area (µm²)",
    "equivalent_voronoi_diameter_um": "Equivalent Voronoi diameter (µm)",
    "spot_fraction_of_voronoi": "Spot fraction of Voronoi cell",
    "total_counts": "Recovered UMI counts",
}


# 1. Read the capture area and spot positions

def transform_geometry(geometry, element, coordinate_system):
    from spatialdata.transformations import get_transformation

    transform = get_transformation(element, get_all=True).get(coordinate_system)
    if transform is None:
        raise KeyError(f"Element has no transformation to {coordinate_system!r}")
    if transform.__class__.__name__ == "Identity":
        return geometry
    if hasattr(transform, "to_affine_matrix"):
        matrix = np.asarray(
            transform.to_affine_matrix(("x", "y"), ("x", "y")), dtype=float
        )
        parameters = [
            matrix[0, 0], matrix[0, 1], matrix[1, 0],
            matrix[1, 1], matrix[0, 2], matrix[1, 2],
        ]
        return affinity.affine_transform(geometry, parameters)
    if hasattr(transform, "translation"):
        return affinity.translate(geometry, *transform.translation[:2])
    raise TypeError(f"Unsupported transformation: {type(transform).__name__}")


def get_capture_area(sdata, sample):
    dataset_id = sample["dataset_id"]
    layer_name = sample.get("capture_area_layer", f"{dataset_id}_capture_area")
    if layer_name not in sdata.shapes:
        raise KeyError(f"Capture-area layer {layer_name!r} is absent")

    layer = sdata.shapes[layer_name]
    geometry = unary_union(layer.geometry)
    geometry = transform_geometry(geometry, layer, dataset_id)
    microns_per_pixel = sample["microns_per_pixel"]
    geometry_um = affinity.scale(
        geometry,
        xfact=microns_per_pixel,
        yfact=microns_per_pixel,
        origin=(0, 0),
    )
    return geometry_um


def get_tissue_area(sdata, sample, capture_area=None):
    dataset_id = sample["dataset_id"]
    layer_name = sample.get("tissue_area_layer", f"{dataset_id}_tissue")
    if layer_name not in sdata.shapes:
        if sample["technology"] == "Curio Seeker":
            if capture_area is None:
                capture_area = get_capture_area(sdata, sample)
            print("  Curio Seeker: treating the complete capture area as tissue")
            return capture_area
        print(f"  Warning: tissue-area layer {layer_name!r} is absent")
        return None
    layer = sdata.shapes[layer_name]
    geometry = transform_geometry(unary_union(layer.geometry), layer, dataset_id)
    scale = sample["microns_per_pixel"]
    return affinity.scale(geometry, xfact=scale, yfact=scale, origin=(0, 0))


def _umi_table_name(sample):
    if "umi_table_layer" in sample:
        return sample["umi_table_layer"]
    technology = sample["technology"]
    dataset_id = sample["dataset_id"]
    if technology.startswith("Visium HD"):
        return "square_002um"
    if technology.startswith("BMKMANU"):
        return f"{dataset_id}_level_1_table"
    if technology == "Stereo-Seq":
        return f"{dataset_id}_raw_table"
    return "table"


def read_umi_observations(sdata, sample):
    table_name = _umi_table_name(sample)
    if table_name not in sdata.tables:
        raise KeyError(f"UMI table {table_name!r} is absent")
    table = sdata.tables[table_name]
    if "total_counts" not in table.obs:
        raise KeyError(f"Table {table_name!r} has no 'total_counts' column")
    if "spatial" not in table.obsm:
        raise KeyError(f"Table {table_name!r} has no obsm['spatial'] coordinates")

    # Copy before unit conversion; AnnData may expose obsm as a writable view.
    observed_xy = np.asarray(table.obsm["spatial"], dtype=float)[:, :2].copy()
    observed_xy *= sample["microns_per_pixel"]
    observed_counts = table.obs["total_counts"].to_numpy(dtype=float)
    valid = np.isfinite(observed_xy).all(axis=1) & np.isfinite(observed_counts)
    id_column = sample.get("umi_spot_id_column")
    if id_column is None:
        preferred_columns = {
            "Curio Seeker": "barcode",
            "BMKMANU S1000": "barcode_name",
            "BMKMANU S3000": "barcode_name",
        }
        candidate = preferred_columns.get(sample["technology"])
        if candidate in table.obs:
            id_column = candidate
    if id_column is not None:
        if id_column not in table.obs:
            raise KeyError(f"Table {table_name!r} has no {id_column!r} identifier column")
        all_observed_ids = table.obs[id_column].astype(str).to_numpy()
    else:
        all_observed_ids = np.asarray(table.obs_names.astype(str))
    observed_ids = all_observed_ids[valid]
    if sample["technology"].startswith("BMKMANU"):
        observed_ids = _canonical_bmk_identifiers(observed_ids)
    elif sample["technology"] == "Nova-ST":
        # NovaScope uses the first 30 bases of the 32-base in-situ whitelist
        # barcode as the expression-table identifier.
        observed_ids = pd.Index(observed_ids).astype(str).str.slice(0, 30).to_numpy()
    observed_xy, observed_counts = observed_xy[valid], observed_counts[valid]

    return pd.DataFrame({
        "spot_id": observed_ids,
        "x_um": observed_xy[:, 0],
        "y_um": observed_xy[:, 1],
        "total_counts": observed_counts,
    })


def add_umi_counts(spots, observations, sample, tissue_area):
    technology = sample["technology"]
    known_identifier_technology = (
        technology == "Curio Seeker"
        or technology == "Nova-ST"
        or technology.startswith("Visium")
        or technology.startswith("BMKMANU")
    )
    identifier_mode = sample.get("umi_identifier_matching", "auto")
    if identifier_mode not in {"auto", True, False}:
        raise ValueError("umi_identifier_matching must be 'auto', True, or False")

    spot_ids = spots["spot_id"].astype(str)
    if technology == "Nova-ST":
        spot_ids = spot_ids.str.slice(0, 30)
    observed_id_matches = observations["spot_id"].isin(set(spot_ids))
    identifier_overlap = observed_id_matches.mean() if len(observations) else 0.0
    minimum_overlap = sample.get("umi_identifier_minimum_overlap", 0.8)
    if identifier_mode is True and identifier_overlap < minimum_overlap:
        raise ValueError(
            f"Only {identifier_overlap:.1%} of expression-table spot identifiers "
            f"match the whitelist; expected at least {minimum_overlap:.1%}."
        )
    use_identifiers = identifier_mode is True or (
        identifier_mode == "auto" and identifier_overlap >= minimum_overlap
    )
    if known_identifier_technology and identifier_mode != False and not use_identifiers:
        raise ValueError(
            f"Only {identifier_overlap:.1%} of {technology} expression-table "
            "spot identifiers match the complete whitelist; expected at least "
            f"{minimum_overlap:.1%}. Check 'umi_spot_id_column' and barcode_file."
        )

    if use_identifiers:
        counts_by_id = observations.groupby("spot_id")["total_counts"].sum()
        matched = spot_ids.isin(counts_by_id.index)
        spots["total_counts"] = (
            spot_ids.map(counts_by_id).fillna(0).to_numpy(dtype=float)
        )
        spots["in_tissue"] = (
            True if technology == "Curio Seeker"
            else shapely.covers(tissue_area, spots.geometry.to_numpy())
            if tissue_area is not None else False
        )
        observed_matched = int(observed_id_matches.sum())
        print(
            f"  Mapped UMI counts by exact spot identifier for "
            f"{observed_matched:,}/{len(observations):,} observed positions; "
            f"{matched.sum():,}/{len(spots):,} whitelist spots had recovered UMIs"
        )
        return spots

    print(
        f"  Exact identifier matching unavailable for {technology} "
        f"({identifier_overlap:.1%} overlap); using coordinates"
    )

    observed_xy = observations[["x_um", "y_um"]].to_numpy(dtype=float)
    observed_counts = observations["total_counts"].to_numpy(dtype=float)
    spot_xy = spots[["x_um", "y_um"]].to_numpy(dtype=float)
    distances, indices = cKDTree(spot_xy).query(observed_xy, k=1)
    tolerance_um = sample.get(
        "umi_mapping_tolerance_um",
        max(0.05, sample["microns_per_pixel"] * 0.75),
    )
    matched = distances <= tolerance_um
    total_counts = np.zeros(len(spots), dtype=float)
    np.add.at(total_counts, indices[matched], observed_counts[matched])
    spots["total_counts"] = total_counts
    spots["in_tissue"] = (
        shapely.covers(tissue_area, spots.geometry.to_numpy())
        if tissue_area is not None else False
    )
    print(
        f"  Mapped UMI counts for {matched.sum():,}/{len(observed_counts):,} "
        f"observed positions (tolerance {tolerance_um:g} µm)"
    )
    return spots


def read_spot_centres(sample, sdata, capture_area):
    technology = sample["technology"]
    dataset_id = sample["dataset_id"]
    barcode_file = Path(sample["barcode_file"])
    microns_per_pixel = sample["microns_per_pixel"]

    # Reconstruct the initial in situ whitelist independently of tissue RNA-seq.
    if technology.startswith("BMKMANU"):
        return _read_bmk_whitelist_centres(
            sample,
            sdata,
            barcode_file,
        )

    if not barcode_file.exists():
        raise FileNotFoundError(barcode_file)

    if technology == "Visium (FFPE)":
        df = pd.read_csv(barcode_file)
        return _spaceranger_centres(df)

    if technology.startswith("Visium HD"):
        df = pd.read_parquet(barcode_file)
        return _spaceranger_centres(df)

    if technology == "Nova-ST":
        return _read_novast_centres(barcode_file, capture_area, microns_per_pixel)

    if technology == "Stereo-Seq":
        return _read_stereoseq_centres(barcode_file, capture_area, microns_per_pixel)

    if technology == "Curio Seeker":
        return _read_curioseeker_centres(sample, barcode_file)

    raise ValueError(f"Unsupported technology: {technology}")


def _spaceranger_centres(df):
    spot_ids = df["barcode"].astype(str) if "barcode" in df else df.index.astype(str)
    return gpd.GeoDataFrame(
        {"spot_id": spot_ids},
        geometry=gpd.points_from_xy(
            df["pxl_col_in_fullres"], df["pxl_row_in_fullres"]
        ),
    )


def _capture_bounds_px(capture_area, microns_per_pixel):
    return tuple(value / microns_per_pixel for value in capture_area.bounds)


def _read_novast_centres(barcode_file, capture_area, microns_per_pixel):
    columns = ["barcode", "lane_id", "tile_id", "x", "y", "count"]
    min_x, min_y, max_x, max_y = _capture_bounds_px(capture_area, microns_per_pixel)
    selected = []
    for chunk in pd.read_csv(
        barcode_file, sep="\t", header=None, names=columns, chunksize=1_000_000
    ):
        # NovaScope coordinates are nanometres; image coordinates are µm.
        inside = (
            chunk["x"].between(min_x * 1000, max_x * 1000)
            & chunk["y"].between(min_y * 1000, max_y * 1000)
        )
        if inside.any():
            selected.append(chunk.loc[inside])
    frame = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame(columns=columns)
    return gpd.GeoDataFrame(
        {"spot_id": frame["barcode"].astype(str)},
        geometry=gpd.points_from_xy(frame["x"] / 1000, frame["y"] / 1000),
    )


def _read_stereoseq_centres(barcode_file, capture_area, microns_per_pixel):
    import h5py

    min_x, min_y, max_x, max_y = _capture_bounds_px(capture_area, microns_per_pixel)
    with h5py.File(barcode_file, "r") as handle:
        matrix = handle["bpMatrix_1"]
        x0, y0 = max(0, math.floor(min_x)), max(0, math.floor(min_y))
        x1 = min(matrix.shape[1], math.ceil(max_x))
        y1 = min(matrix.shape[0], math.ceil(max_y))
        rows, columns = np.nonzero(matrix[y0:y1, x0:x1, 0])
    x, y = columns + x0, rows + y0
    return gpd.GeoDataFrame(
        {"spot_id": np.arange(len(x)).astype(str)},
        geometry=gpd.points_from_xy(x, y),
    )


def _curio_matched_locations_file(sample):
    if "matched_bead_location_file" in sample:
        path = Path(sample["matched_bead_location_file"])
    else:
        # .../<run>/harpy/sdata.zarr -> .../<run>/curioseeker/results/OUTPUT/<id>/
        run_dir = Path(sample["sdata_path"]).parent.parent
        dataset_id = sample["dataset_id"]
        path = (
            run_dir
            / "curioseeker"
            / "results"
            / "OUTPUT"
            / dataset_id
            / f"{dataset_id}_MatchedBeadLocation.csv"
        )
    if not path.exists():
        raise FileNotFoundError(
            "The Curio matched-bead location file is required to register the "
            f"complete whitelist to the SpatialData capture area: {path}"
        )
    return path


def _read_curioseeker_centres(sample, barcode_file):
    df = pd.read_csv(
        barcode_file, sep="\t", header=None, names=["barcode", "y", "x"]
    )
    matched_file = _curio_matched_locations_file(sample)
    matched = pd.read_csv(matched_file, index_col=0)
    required_columns = {"SPATIAL_1", "SPATIAL_2"}
    if not required_columns.issubset(matched.columns):
        raise ValueError(
            f"{matched_file} must contain columns {sorted(required_columns)}"
        )

    x_origin = matched["SPATIAL_2"].max() + CURIO_SOURCE_SPOT_RADIUS_PX
    y_origin = matched["SPATIAL_1"].min() - CURIO_SOURCE_SPOT_RADIUS_PX
    x = x_origin - df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float) - y_origin
    print(
        "  Registered Curio whitelist using recovered-bead origin from "
        f"{matched_file.name}"
    )
    return gpd.GeoDataFrame(
        {"spot_id": df["barcode"].astype(str)},
        geometry=gpd.points_from_xy(x, y),
    )


def _read_bmk_whitelist_centres(sample, sdata, signal_file):
    if not signal_file.exists():
        raise FileNotFoundError(signal_file)
    if "barcode_index_file" not in sample:
        raise KeyError(
            "BMKMANU samples require 'barcode_index_file' in sample_metadata"
        )
    barcode_index_file = Path(sample["barcode_index_file"])
    if not barcode_index_file.exists():
        raise FileNotFoundError(barcode_index_file)

    technology = sample["technology"]
    if technology == "BMKMANU S1000":
        block_rows, block_columns = 46, 46
        spot_rows, spot_columns = 35, 30
        block_row_pitch, block_column_pitch = 36, 31
    elif technology == "BMKMANU S3000":
        block_rows, block_columns = 43, 42
        spot_rows, spot_columns = 51, 45
        block_row_pitch, block_column_pitch = 52, 46
    else:
        raise ValueError(f"Unsupported BMKMANU technology: {technology}")

    expected_positions = block_rows * block_columns * spot_rows * spot_columns
    barcode_lookup = _read_bmk_barcode_index(barcode_index_file)
    barcode_base = max(barcode_lookup.values()) + 1
    barcode_codes = _decode_bmk_signal(signal_file, barcode_lookup, expected_positions)
    whitelist_indices = _unique_decoded_bmk_positions(barcode_codes)
    whitelist_codes = barcode_codes[whitelist_indices].astype(np.uint64)
    barcode_3 = whitelist_codes % barcode_base
    remainder = whitelist_codes // barcode_base
    barcode_2 = remainder % barcode_base
    barcode_1 = remainder // barcode_base
    spot_ids = np.asarray([
        f"bc1_{first}-bc2_{second}-bc3_{third}"
        for first, second, third in zip(barcode_1, barcode_2, barcode_3)
    ])

    within_block_count = spot_rows * spot_columns
    block = whitelist_indices // within_block_count
    block_row = block // block_columns
    block_column = block % block_columns
    within_block = whitelist_indices % within_block_count
    spot_row = within_block // spot_columns
    spot_column = within_block % spot_columns

    x_step = 1000.0 / (block_columns * block_column_pitch)
    y_step = x_step * math.sqrt(3) / 2
    x_standard = x_step * (
        block_column_pitch * block_column
        + spot_column
        + 1
        - 0.5 * (spot_row % 2)
    )
    y_standard = y_step * (
        block_row_pitch * block_row
        + spot_row
        + 1
    )

    image_height, image_width = _bmk_image_size(sdata, sample["dataset_id"])
    standard_height = y_step * block_rows * block_row_pitch
    if 1000.0 / standard_height > image_width / image_height:
        scaling_factor = image_width / 1000.0
    else:
        scaling_factor = image_height / standard_height

    decoded_count = int(np.count_nonzero(barcode_codes != 0))
    undecoded_count = expected_positions - decoded_count
    ambiguous_count = decoded_count - len(whitelist_indices)
    print(f"  Theoretical BMK positions: {expected_positions:,}")
    print(f"  Exactly decoded positions: {decoded_count:,}")
    print(f"  Removed positions without an exact decode: {undecoded_count:,}")
    print(f"  Removed non-unique decoded positions: {ambiguous_count:,}")
    print(f"  Strict in situ whitelist positions: {len(whitelist_indices):,}")

    return gpd.GeoDataFrame(
        {
            "spot_id": spot_ids,
            "signal_index": whitelist_indices,
            "barcode_code": barcode_codes[whitelist_indices],
        },
        geometry=gpd.points_from_xy(
            x_standard * scaling_factor,
            y_standard * scaling_factor,
        ),
    )


def _canonical_bmk_identifiers(values):
    identifiers = pd.Series(values, dtype="string")
    return identifiers.str.replace(
        r"V2[AB](?=bc[123]_)", "", regex=True
    ).to_numpy()


def _read_bmk_barcode_index(barcode_index_file):
    lookup = {}
    with barcode_index_file.open() as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0] == "Index":
                continue
            fluorescence = [int(value) for value in fields[2:14]]
            code = "".join(
                "0" if fluorescence[cycle] == 1
                else "1" if fluorescence[cycle + 6] == 1
                else "2"
                for cycle in range(6)
            )
            lookup[code] = int(fields[0])
    if not lookup:
        raise ValueError(f"No barcode codes found in {barcode_index_file}")
    return lookup


def _decode_bmk_signal(signal_file, barcode_lookup, expected_positions):
    base = max(barcode_lookup.values()) + 1
    barcode_codes = np.zeros(expected_positions, dtype=np.uint32)
    position_count = 0
    with signal_file.open() as handle:
        for position_count, line in enumerate(handle, start=1):
            if position_count > expected_positions:
                raise ValueError(
                    f"{signal_file} has more than {expected_positions:,} positions"
                )
            signal = line.strip()
            if len(signal) != 18:
                raise ValueError(
                    f"Expected an 18-character signal at line {position_count:,}"
                )
            components = (
                barcode_lookup.get(signal[:6]),
                barcode_lookup.get(signal[6:12]),
                barcode_lookup.get(signal[12:]),
            )
            if all(component is not None for component in components):
                first, second, third = components
                barcode_codes[position_count - 1] = (
                    (first * base + second) * base + third
                )
    if position_count != expected_positions:
        raise ValueError(
            f"Expected {expected_positions:,} positions in {signal_file}, "
            f"found {position_count:,}"
        )
    return barcode_codes


def _unique_decoded_bmk_positions(barcode_codes):
    decoded_positions = np.flatnonzero(barcode_codes)
    order = np.argsort(barcode_codes[decoded_positions], kind="stable")
    sorted_positions = decoded_positions[order]
    sorted_codes = barcode_codes[sorted_positions]
    run_start = np.r_[0, np.flatnonzero(np.diff(sorted_codes)) + 1]
    run_stop = np.r_[run_start[1:], len(sorted_codes)]
    unique_runs = (run_stop - run_start) == 1
    return np.sort(sorted_positions[run_start[unique_runs]])


def _bmk_image_size(sdata, dataset_id):
    image = sdata.images[f"{dataset_id}_he"]
    if hasattr(image, "shape") and image.shape is not None:
        return int(image.shape[-2]), int(image.shape[-1])
    try:
        scale_zero = image["scale0"]
        array = next(iter(scale_zero.ds.data_vars.values()))
        return int(array.shape[-2]), int(array.shape[-1])
    except (KeyError, AttributeError, StopIteration) as error:
        raise TypeError(
            f"Could not determine full-resolution size of {dataset_id}_he"
        ) from error


def prepare_spots(sample, sdata):
    if sample["technology"] == "Curio Seeker":
        capture_area = None
    else:
        capture_area = get_capture_area(sdata, sample)
    centres_px = read_spot_centres(sample, sdata, capture_area)
    microns_per_pixel = sample["microns_per_pixel"]
    input_count = len(centres_px)

    valid_centres = centres_px.geometry.notna() & ~centres_px.geometry.is_empty
    invalid_count = int((~valid_centres).sum())
    centres_px = centres_px.loc[valid_centres].copy()
    if centres_px.empty:
        raise ValueError("No valid spot centres were found")

    x_um = centres_px.geometry.x.to_numpy() * microns_per_pixel
    y_um = centres_px.geometry.y.to_numpy() * microns_per_pixel
    spots = gpd.GeoDataFrame(
        centres_px.drop(columns="geometry"),
        geometry=gpd.points_from_xy(x_um, y_um),
    )
    spots["x_um"] = x_um
    spots["y_um"] = y_um
    finite_coordinates = np.isfinite(spots[["x_um", "y_um"]]).all(axis=1)
    non_finite_count = int((~finite_coordinates).sum())
    spots = spots.loc[finite_coordinates].copy()
    if spots.empty:
        raise ValueError("No spot centres with finite coordinates were found")

    if sample["technology"] == "Curio Seeker":
        radius_um = SPOT_CONFIG["Curio Seeker"]["size_um"] / 2
        centre_hull = shapely.convex_hull(
            shapely.multipoints(spots.geometry.to_numpy())
        )
        # The hull of equal circular spot footprints is the hull of their
        # centres expanded by one spot radius.
        capture_area = shapely.buffer(centre_hull, radius_um)
        print("  Created Curio capture area from the complete-whitelist hull")

    inside_capture_area = shapely.covers(
        capture_area,
        spots.geometry.to_numpy(),
    )
    outside_count = int((~inside_capture_area).sum())
    spots = spots.loc[inside_capture_area]
    duplicate_coordinates = spots.duplicated(["x_um", "y_um"])
    duplicate_count = int(duplicate_coordinates.sum())
    spots = spots.loc[~duplicate_coordinates].reset_index(drop=True)

    preparation_counts = {
        "input_spots": input_count,
        "missing_or_empty_centres": invalid_count,
        "non_finite_centres": non_finite_count,
        "outside_capture_area": outside_count,
        "duplicate_coordinates": duplicate_count,
        "retained_spots": len(spots),
    }
    print(f"  Input spot centres: {input_count:,}")
    print(f"  Removed missing or empty centres: {invalid_count:,}")
    print(f"  Removed non-finite centres: {non_finite_count:,}")
    print(f"  Removed outside capture area: {outside_count:,}")
    print(f"  Removed duplicate coordinates: {duplicate_count:,}")
    print(f"  Retained spot centres: {len(spots):,}")

    config = SPOT_CONFIG[sample["technology"]]
    size_um = config["size_um"]
    if config["shape"] == "square":
        half = size_um / 2
        spots["spot_geometry"] = [
            box(x - half, y - half, x + half, y + half)
            for x, y in zip(spots["x_um"], spots["y_um"])
        ]
        spots["spot_area_um2"] = size_um**2
    else:
        radius = size_um / 2
        spots["spot_geometry"] = [
            Point(x, y).buffer(radius)
            for x, y in zip(spots["x_um"], spots["y_um"])
        ]
        spots["spot_area_um2"] = math.pi * radius**2
    spots.attrs["spot_shape"] = config["shape"]
    spots.attrs["spot_size_um"] = size_um
    spots.attrs["technology"] = sample["technology"]
    spots.attrs["preparation_counts"] = preparation_counts
    return spots, capture_area


# 2. Delaunay neighbors and bounded Voronoi cells

def calculate_spatial_metrics(
    spots,
    capture_area,
    max_neighbors=6,
    max_distance_um=None,
):
    total_start = perf_counter()
    timings = {}
    if len(spots) < 4:
        raise ValueError("At least four non-duplicate spots are required")

    metrics = spots.copy()
    coordinates = metrics[["x_um", "y_um"]].to_numpy()

    stage_start = perf_counter()
    try:
        delaunay = Delaunay(coordinates)
    except QhullError as error:
        raise ValueError("Spot centres do not form a valid 2-D tessellation") from error
    timings["delaunay"] = _elapsed(stage_start, "Delaunay triangulation")

    stage_start = perf_counter()
    edges, neighbors = _select_delaunay_edges(
        coordinates, delaunay, max_neighbors, max_distance_um
    )
    timings["neighbor_filtering"] = _elapsed(stage_start, "Neighbor filtering")
    metrics["n_neighbors"] = [len(values) for values in neighbors]
    metrics["mean_neighbor_distance_um"] = [
        np.mean([distance for _, distance in values]) if values else np.nan
        for values in neighbors
    ]
    metrics["neighbor_indices"] = [
        ",".join(str(index) for index, _ in sorted(values)) for values in neighbors
    ]

    stage_start = perf_counter()
    cells = _build_voronoi_cells(coordinates, capture_area)
    timings["voronoi_construction_and_clipping"] = _elapsed(
        stage_start, "Vectorized Voronoi construction and clipping"
    )
    metrics["voronoi_geometry"] = cells
    metrics["voronoi_area_um2"] = shapely.area(cells)
    metrics["equivalent_voronoi_diameter_um"] = 2 * np.sqrt(
        metrics["voronoi_area_um2"] / math.pi
    )
    stage_start = perf_counter()
    spot_geometries = metrics["spot_geometry"].to_numpy()
    overlap_areas = shapely.area(shapely.intersection(spot_geometries, cells))
    metrics["spot_fraction_of_voronoi"] = (
        overlap_areas / metrics["voronoi_area_um2"].to_numpy()
    )
    timings["spot_voronoi_intersections"] = _elapsed(
        stage_start, "Vectorized spot/Voronoi intersections"
    )
    timings["metrics_total"] = perf_counter() - total_start
    print(f"  [timing] Metrics total: {timings['metrics_total']:.2f} s")
    summary = {
        "n_spots": len(metrics),
        "max_distance_um": max_distance_um,
        "max_neighbors": max_neighbors,
        "capture_area_um2": capture_area.area,
        "timings_seconds": timings,
    }
    return metrics, edges, summary


def _select_delaunay_edges(coordinates, delaunay, max_neighbors, max_distance_um):
    triangles = delaunay.simplices
    candidate_pairs = np.vstack(
        (triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]])
    )
    candidate_pairs.sort(axis=1)
    candidate_pairs = np.unique(candidate_pairs, axis=0)
    distances = np.linalg.norm(
        coordinates[candidate_pairs[:, 0]] - coordinates[candidate_pairs[:, 1]],
        axis=1,
    )

    if max_distance_um is not None:
        keep = distances <= max_distance_um
        candidate_pairs = candidate_pairs[keep]
        distances = distances[keep]

    # Sorting is only needed because a degree cap preferentially keeps short edges.
    if max_neighbors is not None:
        order = np.argsort(distances)
        candidate_pairs = candidate_pairs[order]
        distances = distances[order]

    neighbors = [[] for _ in coordinates]
    records = []
    for (first, second), distance in zip(candidate_pairs, distances):
        if max_neighbors is not None and (
            len(neighbors[first]) >= max_neighbors
            or len(neighbors[second]) >= max_neighbors
        ):
            continue
        neighbors[first].append((second, distance))
        neighbors[second].append((first, distance))
        records.append(
            {"source_index": first, "target_index": second, "distance_um": distance}
        )
    edge_table = pd.DataFrame(
        records,
        columns=["source_index", "target_index", "distance_um"],
    )
    return edge_table, neighbors


def _build_voronoi_cells(coordinates, capture_area):
    sites = shapely.multipoints(coordinates)
    cell_collection = shapely.voronoi_polygons(
        sites,
        extend_to=shapely.envelope(capture_area),
        ordered=True,
    )
    cells = shapely.get_parts(cell_collection)
    if len(cells) != len(coordinates):
        raise ValueError(
            f"Expected {len(coordinates)} Voronoi cells, received {len(cells)}"
        )

    cells = _clip_voronoi_cells(cells, capture_area)
    if np.any(shapely.is_empty(cells)):
        raise ValueError("A Voronoi cell does not overlap the capture area")
    return cells


def _clip_voronoi_cells(cells, capture_area, chunk_size=100_000):
    clipped_chunks = []
    repaired_capture_area = None
    for start in range(0, len(cells), chunk_size):
        chunk = cells[start : start + chunk_size]
        try:
            clipped = shapely.intersection(chunk, capture_area)
        except GEOSException:
            if repaired_capture_area is None:
                repaired_capture_area = shapely.make_valid(capture_area)
            clipped = shapely.intersection(
                shapely.make_valid(chunk),
                repaired_capture_area,
            )
        clipped_chunks.append(clipped)
    return np.concatenate(clipped_chunks)


# 3. Coverage, plots, and output files

def calculate_coverage_curve(
    spots,
    capture_area,
    distances_um,
    tile_size_um=1000.0,
    saturation_percent=99.9,
):
    distances_um = np.asarray(distances_um, dtype=float)
    if distances_um.size == 0:
        raise ValueError("At least one coverage distance is required")
    if np.any(distances_um < 0):
        raise ValueError("Coverage distances must be non-negative")
    if np.any(np.diff(distances_um) < 0):
        raise ValueError("Coverage distances must be in ascending order")
    if not 0 < saturation_percent <= 100:
        raise ValueError("saturation_percent must be greater than 0 and at most 100")

    initial_coverage_percent = _initial_spot_coverage_percent(
        spots,
        capture_area,
    )
    print(f"  Initial spot coverage: {initial_coverage_percent:.3f}%")
    if initial_coverage_percent >= saturation_percent:
        print(
            f"  Coverage is at least {saturation_percent:g}%; "
            "setting the complete curve to 100%"
        )
        return pd.DataFrame(
            {
                "distance_from_spot_edge_um": distances_um,
                "coverage_percent": np.full(len(distances_um), 100.0),
            }
        )

    if spots.attrs.get("spot_shape") == "circle":
        print("  Using vectorized Voronoi coverage for circular spots")
        return _calculate_circular_coverage(
            spots,
            capture_area,
            distances_um,
            saturation_percent=saturation_percent,
        )

    print("  Using tiled geometry coverage for square spots")
    coverage = _calculate_tiled_coverage(
        spots,
        capture_area,
        distances_um,
        tile_size_um,
    )
    return _round_saturated_coverage(coverage, saturation_percent)


def _initial_spot_coverage_percent(spots, capture_area):
    required_columns = {"voronoi_area_um2", "spot_fraction_of_voronoi"}
    missing_columns = required_columns.difference(spots.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(
            f"Initial coverage requires {missing}; run "
            "calculate_spatial_metrics() first"
        )
    covered_area_um2 = np.sum(
        spots["voronoi_area_um2"].to_numpy()
        * spots["spot_fraction_of_voronoi"].to_numpy()
    )
    return 100 * covered_area_um2 / capture_area.area


def _round_saturated_coverage(coverage, saturation_percent):
    saturated = coverage["coverage_percent"].to_numpy() >= saturation_percent
    if saturated.any():
        first_saturated = int(np.flatnonzero(saturated)[0])
        coverage.loc[first_saturated:, "coverage_percent"] = 100.0
    return coverage


def _calculate_circular_coverage(
    spots,
    capture_area,
    distances_um,
    chunk_size=100_000,
    saturation_percent=99.9,
):
    if "voronoi_geometry" not in spots:
        raise ValueError(
            "Circular coverage requires Voronoi cells; run "
            "calculate_spatial_metrics() first"
        )
    if "spot_size_um" not in spots.attrs:
        raise ValueError("Circular coverage requires spots.attrs['spot_size_um']")

    centres = spots.geometry.to_numpy()
    cells = spots["voronoi_geometry"].to_numpy()
    radius_um = spots.attrs["spot_size_um"] / 2
    covered_area_um2 = np.zeros(len(distances_um), dtype=float)

    for distance_index, distance_um in enumerate(distances_um):
        expanded_radius_um = radius_um + float(distance_um)
        for start in range(0, len(spots), chunk_size):
            stop = min(start + chunk_size, len(spots))
            expanded_spots = shapely.buffer(
                centres[start:stop],
                expanded_radius_um,
                quad_segs=16,
            )
            covered_parts = shapely.intersection(
                expanded_spots,
                cells[start:stop],
            )
            covered_area_um2[distance_index] += shapely.area(
                covered_parts
            ).sum()

        coverage_percent = (
            100 * covered_area_um2[distance_index] / capture_area.area
        )
        if coverage_percent >= saturation_percent:
            covered_area_um2[distance_index:] = capture_area.area
            print(
                f"  Coverage reached {coverage_percent:.3f}% at "
                f"{distance_um:g} µm; setting this and remaining values to 100%"
            )
            break

    return pd.DataFrame(
        {
            "distance_from_spot_edge_um": distances_um,
            "coverage_percent": 100 * covered_area_um2 / capture_area.area,
        }
    )


def _calculate_tiled_coverage(
    spots,
    capture_area,
    distances_um,
    tile_size_um,
):
    if tile_size_um <= 0:
        raise ValueError("tile_size_um must be greater than 0")

    maximum_distance = float(distances_um.max())
    covered_area_um2 = np.zeros(len(distances_um), dtype=float)
    spot_geometries = gpd.GeoSeries(
        spots["spot_geometry"].to_numpy(),
        index=spots.index,
    )
    spatial_index = spot_geometries.sindex
    min_x, min_y, max_x, max_y = capture_area.bounds

    for x_start in np.arange(min_x, max_x, tile_size_um):
        for y_start in np.arange(min_y, max_y, tile_size_um):
            tile = box(
                x_start,
                y_start,
                min(x_start + tile_size_um, max_x),
                min(y_start + tile_size_um, max_y),
            )
            capture_tile = capture_area.intersection(tile)
            if capture_tile.is_empty or capture_tile.area == 0:
                continue

            halo = tile.buffer(maximum_distance, cap_style="square", join_style="mitre")
            candidate_positions = spatial_index.query(halo, predicate="intersects")
            if len(candidate_positions) == 0:
                continue

            local_spot_union = unary_union(
                spot_geometries.iloc[candidate_positions]
            )
            for index, distance in enumerate(distances_um):
                expanded_spots = (
                    local_spot_union
                    if distance == 0
                    else local_spot_union.buffer(float(distance))
                )
                covered_area_um2[index] += expanded_spots.intersection(
                    capture_tile
                ).area

    coverage_percent = 100 * covered_area_um2 / capture_area.area
    return pd.DataFrame(
        {
            "distance_from_spot_edge_um": distances_um,
            "coverage_percent": coverage_percent,
        }
    )


def save_plots(
    metrics,
    edges,
    coverage,
    capture_area,
    tissue_area,
    output_dir,
    crop_size_um,
    spatial_bin_size_um,
    dpi,
):
    crop_boxes = _make_crop_boxes(capture_area, crop_size_um)
    for metric, label in PLOT_METRICS.items():
        _save_histogram(
            metrics[metric],
            label,
            output_dir / f"{metric}_histogram.png",
            dpi,
        )
        for crop_name, crop_box in crop_boxes.items():
            _save_spatial_metric_plot(
                metrics,
                edges,
                metric,
                label,
                crop_box,
                crop_name,
                tissue_area,
                output_dir / f"{metric}_spatial_{crop_name}.png",
                dpi,
            )
    _save_capture_area_binned_maps(
        metrics,
        capture_area,
        tissue_area,
        output_dir,
        dpi,
        bin_size_um=spatial_bin_size_um,
    )
    _save_cdf_plot(metrics, output_dir / "distance_cdfs.png", dpi)
    _save_umi_cdf_plot(
        metrics, output_dir / "umi_count_cdf_capture_area.png", dpi,
        title="Complete capture area",
    )
    if "in_tissue" in metrics and metrics["in_tissue"].any():
        _save_umi_cdf_plot(
            metrics.loc[metrics["in_tissue"]],
            output_dir / "umi_count_cdf_tissue.png",
            dpi,
            title="Tissue area",
        )
    _save_coverage_plot(coverage, output_dir / "spatial_coverage_curve.png", dpi)


def _make_crop_boxes(capture_area, crop_size_um):
    half = crop_size_um / 2
    center = capture_area.centroid
    min_x, min_y, _, _ = capture_area.bounds
    top_left = nearest_points(capture_area.boundary, Point(min_x, min_y))[0]

    return {
        "center": box(
            center.x - half, center.y - half,
            center.x + half, center.y + half,
        ),
        "top_left": box(
            top_left.x - (crop_size_um * 0.25), top_left.y - (crop_size_um * 0.25),
            top_left.x + (crop_size_um * 0.75), top_left.y + (crop_size_um * 0.75),
        ),
    }


def _save_histogram(values, label, output_file, dpi):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    mean = values.mean()
    median = np.median(values)
    if np.all(values == np.floor(values)) and values.max() - values.min() <= 100:
        bins = np.arange(values.min() - 0.5, values.max() + 1.5)
    elif np.isclose(values.min(), values.max(), rtol=1e-10, atol=1e-12):
        bins = 1
    else:
        bins = np.linspace(values.min(), values.max(), 51)

    figure, axis = plt.subplots(figsize=(6, 4.5))
    axis.hist(
        values,
        bins=bins,
        color="#4c91c6",
        edgecolor="#17324d",
        linewidth=0.5,
        alpha=0.9,
    )
    axis.axvline(mean, color="#d1495b", linestyle="--", linewidth=1.8,
                 label=f"Mean: {mean:.3g}")
    axis.axvline(median, color="#2a9d8f", linestyle=":", linewidth=2,
                 label=f"Median: {median:.3g}")
    axis.set(xlabel=label, ylabel="Spots")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", alpha=0.2)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _save_spatial_metric_plot(
    metrics,
    edges,
    metric,
    label,
    crop_box,
    crop_name,
    tissue_area,
    output_file,
    dpi,
):
    min_x, min_y, max_x, max_y = crop_box.bounds
    inside = (
        metrics["x_um"].between(min_x, max_x)
        & metrics["y_um"].between(min_y, max_y)
    )
    indices = np.flatnonzero(inside.to_numpy())
    if len(indices) == 0:
        print(f"Skipping empty {crop_name} crop for {metric}")
        return

    values = metrics[metric].to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    if metric == "total_counts":
        vmin, vmax = 0.0, np.nanpercentile(finite_values, 99)
    else:
        vmin, vmax = np.nanpercentile(finite_values, [1, 99])
    if vmin == vmax:
        vmax = vmin + 1
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("viridis").copy()
    display_values = values.copy()
    if metric == "total_counts":
        cmap.set_under("black")
        display_values[display_values == 0] = -max(1.0, float(vmax)) * 1e-6

    figure, axis = plt.subplots(figsize=(7, 7))
    _add_voronoi_lines(axis, metrics.iloc[indices])
    spot_collection = _add_true_size_spots(
        axis, metrics.iloc[indices], display_values[indices], norm, cmap
    )
    _add_neighbor_lines(axis, metrics, edges, crop_box)
    if tissue_area is not None:
        _add_geometry_outline(axis, tissue_area.intersection(crop_box), "red", 1.0, 6)

    colorbar = figure.colorbar(spot_collection, ax=axis, shrink=0.8)
    colorbar.set_label(label)
    _show_full_colorbar_values(colorbar)
    axis.set(
        xlim=(min_x, max_x),
        ylim=(max_y, min_y),
        xlabel="x (µm)",
        ylabel="y (µm)",
        title=f"{label} — {crop_name.replace('_', ' ')} crop",
        aspect="equal",
    )
    figure.tight_layout()
    figure.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _add_neighbor_lines(axis, metrics, edges, crop_box):
    coordinates = metrics[["x_um", "y_um"]].to_numpy()
    pairs = edges[["source_index", "target_index"]].to_numpy(dtype=int)
    segments = coordinates[pairs]
    min_x, min_y, max_x, max_y = crop_box.bounds
    visible = (
        (segments[:, :, 0].max(axis=1) >= min_x)
        & (segments[:, :, 0].min(axis=1) <= max_x)
        & (segments[:, :, 1].max(axis=1) >= min_y)
        & (segments[:, :, 1].min(axis=1) <= max_y)
    )
    axis.add_collection(LineCollection(
        segments[visible], colors="#d1495b", linewidths=0.50,
        alpha=0.8, zorder=5,
    ))


def _show_full_colorbar_values(colorbar):
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    colorbar.formatter = formatter
    colorbar.update_ticks()


def _add_voronoi_lines(axis, cropped_metrics):
    segments = []
    for geometry in cropped_metrics["voronoi_geometry"]:
        for polygon in _polygon_parts(geometry):
            segments.append(np.asarray(polygon.exterior.coords))
    axis.add_collection(LineCollection(
        segments, colors="0.35", linewidths=0.45, alpha=0.75, zorder=2
    ))


def _polygon_parts(geometry):
    if geometry is None or geometry.is_empty:
        return
    if geometry.geom_type == "Polygon":
        yield geometry
        return
    if geometry.geom_type in {"MultiPolygon", "GeometryCollection"}:
        for part in geometry.geoms:
            yield from _polygon_parts(part)


def _add_geometry_outline(axis, geometry, color, linewidth=1.0, zorder=6):
    for polygon in _polygon_parts(geometry):
        exterior = np.asarray(polygon.exterior.coords)
        axis.plot(exterior[:, 0], exterior[:, 1], color=color,
                  linewidth=linewidth, zorder=zorder)
        for interior in polygon.interiors:
            ring = np.asarray(interior.coords)
            axis.plot(ring[:, 0], ring[:, 1], color=color,
                      linewidth=linewidth, zorder=zorder)


def _save_capture_area_binned_maps(
    metrics,
    capture_area,
    tissue_area,
    output_dir,
    dpi,
    bin_size_um=200.0,
):
    min_x, min_y, max_x, max_y = capture_area.bounds
    x_edges = _bin_edges(min_x, max_x, bin_size_um)
    y_edges = _bin_edges(min_y, max_y, bin_size_um)
    n_x = len(x_edges) - 1
    n_y = len(y_edges) - 1

    x_bin = np.searchsorted(x_edges, metrics["x_um"].to_numpy(), side="right") - 1
    y_bin = np.searchsorted(y_edges, metrics["y_um"].to_numpy(), side="right") - 1
    # Include points exactly on the maximum bounds in the final bins.
    x_bin[x_bin == n_x] = n_x - 1
    y_bin[y_bin == n_y] = n_y - 1
    valid_bins = (x_bin >= 0) & (x_bin < n_x) & (y_bin >= 0) & (y_bin < n_y)
    flat_bin = y_bin[valid_bins] * n_x + x_bin[valid_bins]
    n_bins = n_x * n_y
    counts = np.bincount(flat_bin, minlength=n_bins).astype(float)

    bin_polygons = np.asarray(
        [
            box(x_edges[x], y_edges[y], x_edges[x + 1], y_edges[y + 1])
            for y in range(n_y)
            for x in range(n_x)
        ],
        dtype=object,
    )
    capture_bin_area_um2 = shapely.area(
        shapely.intersection(bin_polygons, capture_area)
    )
    density = np.divide(
        counts * 1_000_000,
        capture_bin_area_um2,
        out=np.full(n_bins, np.nan),
        where=capture_bin_area_um2 > 0,
    )

    panels = [(density, "spot_density", "Spot density", "spots/mm²")]
    for column, file_label, title, colorbar_label in (
        (
            "n_neighbors",
            "average_number_of_neighbors",
            "Average number of neighbors",
            "neighbors",
        ),
        (
            "equivalent_voronoi_diameter_um",
            "average_equivalent_voronoi_diameter",
            "Average equivalent Voronoi diameter",
            "µm",
        ),
        (
            "mean_neighbor_distance_um",
            "average_neighbor_distance",
            "Average neighbor distance",
            "µm",
        ),
        (
            "voronoi_area_um2",
            "average_voronoi_area",
            "Average Voronoi area",
            "µm²",
        ),
        (
            "spot_fraction_of_voronoi",
            "average_spot_fraction_of_voronoi",
            "Average spot fraction of Voronoi cell",
            "fraction",
        ),
        (
            "total_counts",
            "average_umi_counts",
            "Average recovered UMI counts",
            "UMIs",
        ),
    ):
        values = metrics[column].to_numpy(dtype=float)[valid_bins]
        finite = np.isfinite(values)
        sums = np.bincount(
            flat_bin[finite],
            weights=values[finite],
            minlength=n_bins,
        )
        observations = np.bincount(flat_bin[finite], minlength=n_bins)
        means = np.divide(
            sums,
            observations,
            out=np.full(n_bins, np.nan),
            where=observations > 0,
        )
        panels.append((means, file_label, title, colorbar_label))

    for values, file_label, title, colorbar_label in panels:
        figure, axis = plt.subplots(figsize=(7, 7))
        image = axis.pcolormesh(
            x_edges,
            y_edges,
            values.reshape(n_y, n_x),
            shading="flat",
            cmap="viridis",
        )
        for polygon in _polygon_parts(capture_area):
            boundary = np.asarray(polygon.exterior.coords)
            axis.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=0.5)
        if tissue_area is not None:
            _add_geometry_outline(axis, tissue_area, "red", 1.0, 4)
        colorbar = figure.colorbar(image, ax=axis, shrink=0.82)
        colorbar.set_label(colorbar_label)
        _show_full_colorbar_values(colorbar)
        axis.set(
            title=title,
            xlabel="x (µm)",
            ylabel="y (µm)",
            xlim=(min_x, max_x),
            ylim=(max_y, min_y),
            aspect="equal",
        )
        figure.tight_layout()
        figure.savefig(
            output_dir / f"{file_label}_{bin_size_um:g}um_bins.png",
            dpi=dpi,
            bbox_inches="tight",
        )
        plt.close(figure)


def _bin_edges(minimum, maximum, bin_size):
    edges = np.arange(minimum, maximum, bin_size, dtype=float)
    if edges.size == 0 or edges[0] != minimum:
        edges = np.insert(edges, 0, minimum)
    if edges[-1] < maximum:
        edges = np.append(edges, maximum)
    return edges


def _add_true_size_spots(axis, cropped_metrics, values, norm, cmap="viridis"):
    coordinates = cropped_metrics[["x_um", "y_um"]].to_numpy()
    size_um = cropped_metrics.attrs["spot_size_um"]
    if cropped_metrics.attrs.get("technology") == "Curio Seeker":
        size_um *= 0.9
    shape = cropped_metrics.attrs["spot_shape"]

    if shape == "circle":
        collection = EllipseCollection(
            widths=size_um,
            heights=size_um,
            angles=0,
            units="xy",
            offsets=coordinates,
            transOffset=axis.transData,
            cmap=cmap,
            norm=norm,
            edgecolors="black",
            linewidths=0.25,
            zorder=4,
        )
    else:
        half = size_um / 2
        vertices = [
            [(x - half, y - half), (x + half, y - half),
             (x + half, y + half), (x - half, y + half)]
            for x, y in coordinates
        ]
        collection = PolyCollection(
            vertices,
            cmap=cmap,
            norm=norm,
            edgecolors="black",
            linewidths=0.25,
            zorder=4,
        )
    collection.set_array(values)
    axis.add_collection(collection)
    return collection


def _empirical_cdf(values):
    x = np.sort(np.asarray(values, dtype=float))
    x = x[np.isfinite(x)]
    return x, np.arange(1, len(x) + 1) / len(x)


def _save_cdf_plot(metrics, output_file, dpi):
    figure, axis = plt.subplots(figsize=(6, 4.5))
    for metric, label in (
        ("equivalent_voronoi_diameter_um", "Equivalent Voronoi diameter"),
        ("mean_neighbor_distance_um", "Mean neighbor distance"),
    ):
        x, y = _empirical_cdf(metrics[metric])
        axis.plot(x, y, label=label)
    axis.set(xlabel="Distance (µm)", ylabel="Cumulative fraction", ylim=(0, 1))
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _save_umi_cdf_plot(metrics, output_file, dpi, title):
    x, y = _empirical_cdf(metrics["total_counts"])
    figure, axis = plt.subplots(figsize=(6, 4.5))
    axis.step(x, y, where="post")
    axis.set(
        xlabel="Recovered UMI counts per spot",
        ylabel="Cumulative fraction",
        ylim=(0, 1),
        title=title,
    )
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _save_coverage_plot(coverage, output_file, dpi):
    figure, axis = plt.subplots(figsize=(6, 4.5))
    axis.plot(coverage["distance_from_spot_edge_um"], coverage["coverage_percent"])
    axis.set(
        xlabel="Distance from spot edge (µm)",
        ylabel="Capture area covered (%)",
        ylim=(0, 100),
    )
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_file, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _save_tables(metrics, edges, coverage, output_dir):
    geometry_columns = {"geometry", "spot_geometry", "voronoi_geometry"}
    value_columns = [column for column in metrics if column not in geometry_columns]
    metrics[value_columns].to_parquet(output_dir / "spot_metrics.parquet", index=False)
    edges.to_parquet(output_dir / "delaunay_edges.parquet", index=False)
    voronoi_cells = gpd.GeoDataFrame(
        metrics[value_columns], geometry=metrics["voronoi_geometry"]
    )
    voronoi_cells.to_parquet(output_dir / "voronoi_cells.parquet", index=False)
    coverage.to_csv(output_dir / "spatial_coverage_curve.csv", index=False)


# 4. Run sample

def run_sample(
    sample,
    output_dir,
    max_neighbors=6,
    max_distance_um=None,
    crop_size_um=120.0,
    coverage_max_um=50.0,
    coverage_steps=51,
    coverage_tile_size_um=1000.0,
    spatial_bin_size_um=200.0,
    dpi=300,
    scalable=None,
    processing_tile_size_um=500.0,
    processing_halo_um=25.0,
):
    import spatialdata as sd

    pipeline_start = perf_counter()
    pipeline_timings = {}
    if max_neighbors is not None and max_neighbors < 1:
        raise ValueError("max_neighbors must be None or at least 1")
    if max_distance_um is not None and max_distance_um <= 0:
        raise ValueError("max_distance_um must be None or greater than 0")
    if spatial_bin_size_um <= 0:
        raise ValueError("spatial_bin_size_um must be greater than 0")

    dataset_id = sample["dataset_id"]
    technology = sample["technology"]
    sample_type = sample["sample_type"]
    replicate = sample["replicate"]
    coverage_saturation_percent = (
        98.0 if technology.startswith("Visium HD") else 99.9
    )

    if replicate:
        sample_output = Path(output_dir) / sample_type / f"{technology}_{replicate}"
    else:
        sample_output = Path(output_dir) / sample_type / technology
    sample_output.mkdir(parents=True, exist_ok=True)

    print(f"Reading {technology} - {dataset_id}")
    stage_start = perf_counter()
    sdata = sd.read_zarr(sample["sdata_path"])
    pipeline_timings["read_spatialdata"] = _elapsed(
        stage_start, "Read SpatialData"
    )

    use_scalable = (
        technology in {"Nova-ST", "Stereo-Seq"}
        if scalable is None else scalable
    )
    if use_scalable:
        if technology not in {"Nova-ST", "Stereo-Seq"}:
            raise ValueError("scalable=True is currently supported only for Nova-ST and Stereo-Seq")
        try:
            from .tiled_spatial_sampling import run_tiled_sample
        except ImportError:
            from tiled_spatial_sampling import run_tiled_sample

        capture_area = get_capture_area(sdata, sample)
        tissue_area = get_tissue_area(sdata, sample, capture_area)
        umi_observations = read_umi_observations(sdata, sample)
        print("Using disk-backed tiled processing")
        return run_tiled_sample(
            sample,
            capture_area,
            sample_output,
            max_neighbors=max_neighbors,
            max_distance_um=max_distance_um,
            processing_tile_size_um=processing_tile_size_um,
            processing_halo_um=processing_halo_um,
            coverage_max_um=coverage_max_um,
            coverage_steps=coverage_steps,
            spatial_bin_size_um=spatial_bin_size_um,
            crop_size_um=crop_size_um,
            dpi=dpi,
            saturation_percent=coverage_saturation_percent,
            tissue_area=tissue_area,
            umi_observations=umi_observations,
            umi_mapping_tolerance_um=sample.get(
                "umi_mapping_tolerance_um",
                max(0.05, sample["microns_per_pixel"] * 0.75),
            ),
        )

    print("Reading spot positions and capture area")
    stage_start = perf_counter()
    spots, capture_area = prepare_spots(sample, sdata)
    tissue_area = get_tissue_area(sdata, sample, capture_area)
    observations = read_umi_observations(sdata, sample)
    spots = add_umi_counts(spots, observations, sample, tissue_area)
    pipeline_timings["prepare_spots"] = _elapsed(stage_start, "Prepare spots")

    print(f"Calculating metrics for {len(spots):,} spots")
    metrics, edges, summary = calculate_spatial_metrics(
        spots,
        capture_area,
        max_neighbors=max_neighbors,
        max_distance_um=max_distance_um,
    )

    print("Calculating spatial coverage")
    stage_start = perf_counter()
    distances = np.linspace(0, coverage_max_um, coverage_steps)
    coverage = calculate_coverage_curve(
        metrics,
        capture_area,
        distances,
        tile_size_um=coverage_tile_size_um,
        saturation_percent=coverage_saturation_percent,
    )
    pipeline_timings["spatial_coverage"] = _elapsed(
        stage_start, "Spatial coverage"
    )

    print(f"Saving results to {sample_output}")
    stage_start = perf_counter()
    _save_tables(metrics, edges, coverage, sample_output)
    pipeline_timings["save_tables"] = _elapsed(stage_start, "Save tables")

    stage_start = perf_counter()
    save_plots(
        metrics,
        edges,
        coverage,
        capture_area,
        None if technology == "Curio Seeker" else tissue_area,
        sample_output,
        crop_size_um,
        spatial_bin_size_um,
        dpi,
    )
    pipeline_timings["save_plots"] = _elapsed(stage_start, "Save plots")

    pipeline_timings.update(summary["timings_seconds"])
    pipeline_timings["pipeline_total"] = perf_counter() - pipeline_start
    summary.update(
        dataset_id=dataset_id,
        technology=technology,
        status="complete",
        preparation_counts=spots.attrs["preparation_counts"],
        timings_seconds=pipeline_timings,
    )
    (sample_output / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(f"Done in {pipeline_timings['pipeline_total']:.2f} s")
    return summary
