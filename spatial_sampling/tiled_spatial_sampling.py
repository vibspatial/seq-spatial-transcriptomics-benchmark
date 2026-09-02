from __future__ import annotations

import json
import math
from pathlib import Path
from time import perf_counter

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from scipy.spatial import Delaunay, cKDTree
import shapely
from shapely.errors import GEOSException
from shapely.geometry import box


METRIC_COLUMNS = [
    "n_neighbors", "mean_neighbor_distance_um", "nearest_neighbor_distance_um",
    "voronoi_area_um2", "equivalent_voronoi_diameter_um",
    "spot_fraction_of_voronoi",
]


def _write_partition(frame, root, tile_y, tile_x, part):
    target = root / f"tile_y={tile_y}" / f"tile_x={tile_x}"
    target.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), target / f"part-{part:07d}.parquet")


def _stage_frame(frame, stage_dir, bounds, tile_size, part):
    min_x, min_y, _, _ = bounds
    frame["tile_x"] = np.floor((frame["x_um"] - min_x) / tile_size).astype(np.int32)
    frame["tile_y"] = np.floor((frame["y_um"] - min_y) / tile_size).astype(np.int32)
    for (tile_y, tile_x), group in frame.groupby(["tile_y", "tile_x"], sort=False):
        _write_partition(group.drop(columns=["tile_x", "tile_y"]), stage_dir, int(tile_y), int(tile_x), part)


def stage_positions(sample, capture_area, stage_dir, tile_size_um):
    complete = stage_dir / "_SUCCESS"
    if complete.exists():
        return json.loads((stage_dir / "stage_summary.json").read_text())
    stage_dir.mkdir(parents=True, exist_ok=True)
    technology = sample["technology"]
    scale = sample["microns_per_pixel"]
    bounds = capture_area.bounds
    total = retained = 0

    if technology == "Nova-ST":
        names = ["barcode", "lane_id", "tile_id", "x", "y", "count"]
        reader = pd.read_csv(
            sample["barcode_file"], sep="\t", header=None, names=names,
            usecols=[0, 3, 4], dtype={"barcode": "string", "x": "float64", "y": "float64"},
            chunksize=500_000,
        )
        for part, chunk in enumerate(reader):
            total += len(chunk)
            x = chunk["x"].to_numpy() * scale / 1000
            y = chunk["y"].to_numpy() * scale / 1000
            keep = np.isfinite(x) & np.isfinite(y)
            points = shapely.points(x[keep], y[keep])
            inside = shapely.covers(capture_area, points)
            indices = np.flatnonzero(keep)[inside]
            if not len(indices):
                continue
            out = pd.DataFrame({
                "spot_id": chunk["barcode"].iloc[indices].astype(str).to_numpy(),
                "x_um": x[indices], "y_um": y[indices],
            })
            retained += len(out)
            _stage_frame(out, stage_dir, bounds, tile_size_um, part)
    elif technology == "Stereo-Seq":
        min_x, min_y, max_x, max_y = (v / scale for v in bounds)
        with h5py.File(sample["barcode_file"], "r") as handle:
            matrix = handle["bpMatrix_1"]
            x0, y0 = max(0, math.floor(min_x)), max(0, math.floor(min_y))
            x1, y1 = min(matrix.shape[1], math.ceil(max_x)), min(matrix.shape[0], math.ceil(max_y))
            block = matrix.chunks[0] if matrix.chunks else 1024
            part = 0
            for row0 in range(y0, y1, block):
                for col0 in range(x0, x1, block):
                    values = matrix[row0:min(row0 + block, y1), col0:min(col0 + block, x1), 0]
                    rows, cols = np.nonzero(values)
                    total += values.size
                    if not len(rows):
                        part += 1
                        continue
                    x_px, y_px = cols + col0, rows + row0
                    x, y = x_px * scale, y_px * scale
                    inside = shapely.covers(capture_area, shapely.points(x, y))
                    x, y = x[inside], y[inside]
                    out = pd.DataFrame({
                        "spot_id": y_px[inside].astype(np.int64) * matrix.shape[1] + x_px[inside],
                        "x_um": x, "y_um": y,
                    })
                    retained += len(out)
                    if len(out):
                        _stage_frame(out, stage_dir, bounds, tile_size_um, part)
                    part += 1
    else:
        raise ValueError(f"Tiled staging is unsupported for {technology}")

    summary = {"input_positions_or_pixels": total, "retained_spots": retained, "tile_size_um": tile_size_um}
    (stage_dir / "stage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    complete.touch()
    return summary


def _tile_keys(stage_dir):
    keys = set()
    for path in stage_dir.glob("tile_y=*/tile_x=*"):
        keys.add((int(path.parent.name.split("=")[1]), int(path.name.split("=")[1])))
    return sorted(keys)


def _read_tiles(stage_dir, keys):
    frames = []
    for tile_y, tile_x in keys:
        path = stage_dir / f"tile_y={tile_y}" / f"tile_x={tile_x}"
        if path.exists():
            frame = pd.read_parquet(path)
            frame["_tile_y"] = tile_y
            frame["_tile_x"] = tile_x
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["spot_id", "x_um", "y_um", "_tile_y", "_tile_x"])
    # Duplicate coordinates necessarily have identical tile ownership. Keeping
    # the first here avoids a global in-memory duplicate table.
    return pd.concat(frames, ignore_index=True).drop_duplicates(["x_um", "y_um"], keep="first")


def _neighbor_data(coordinates, simplices, max_distance):
    pairs = np.sort(np.concatenate((simplices[:, [0, 1]], simplices[:, [1, 2]], simplices[:, [0, 2]])), axis=1)
    pairs = np.unique(pairs, axis=0)
    delta = coordinates[pairs[:, 0]] - coordinates[pairs[:, 1]]
    distances = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    if max_distance is not None:
        keep = distances <= max_distance
        pairs, distances = pairs[keep], distances[keep]
    return pairs, distances


def _voronoi(coordinates, clip_area):
    cells = shapely.get_parts(shapely.voronoi_polygons(
        shapely.multipoints(coordinates), extend_to=shapely.envelope(clip_area), ordered=True
    ))
    try:
        return shapely.intersection(cells, clip_area)
    except GEOSException:
        # Capture-area polygons originating in vendor files can contain invalid
        # rings (for example, a hole detached from its shell). Repair geometry
        # only on the exceptional path because make_valid is relatively costly.
        return shapely.intersection(
            shapely.make_valid(cells),
            shapely.make_valid(clip_area),
        )


def process_tiles(sample, capture_area, stage_dir, output_dir, tile_size_um, halo_um,
                  max_distance_um, tissue_area=None, umi_observations=None,
                  umi_mapping_tolerance_um=0.05):
    metrics_dir, voronoi_dir, edges_dir = (output_dir / x for x in ("spot_metrics", "voronoi_cells", "delaunay_edges"))
    done_dir = output_dir / "tile_markers"; done_dir.mkdir(parents=True, exist_ok=True)
    min_x, min_y, _, _ = capture_area.bounds
    radius = {"Nova-ST": .15, "Stereo-Seq": .11}[sample["technology"]]
    umi_tree = None
    umi_counts = None
    umi_by_id = None
    if umi_observations is not None and len(umi_observations):
        if sample["technology"] == "Nova-ST":
            umi_by_id = umi_observations.groupby("spot_id")["total_counts"].sum()
            print(
                f"  Mapping Nova-ST UMI counts by exact spot identifier "
                f"({len(umi_by_id):,} observed identifiers)",
                flush=True,
            )
        else:
            # Stereo's theoretical mask does not preserve expression-table IDs.
            print("  Mapping Stereo-Seq UMI counts by coordinate", flush=True)
            grouped = umi_observations.groupby(["x_um", "y_um"], as_index=False)["total_counts"].sum()
            umi_tree = cKDTree(grouped[["x_um", "y_um"]].to_numpy(float))
            umi_counts = grouped["total_counts"].to_numpy(float)
    keys = _tile_keys(stage_dir)
    for number, (tile_y, tile_x) in enumerate(keys, 1):
        marker = done_dir / f"tile_y{tile_y}_x{tile_x}.json"
        if marker.exists():
            metric_part = (
                metrics_dir
                / f"tile_y={tile_y}"
                / f"tile_x={tile_x}"
                / "part-0000000.parquet"
            )
            if (
                metric_part.exists()
                and {"halo_boundary_affected", "total_counts", "in_tissue"}.issubset(
                    pq.ParquetFile(metric_part).schema.names
                )
            ):
                continue
            print(
                f"  Reprocessing tile ({tile_y}, {tile_x}) to add the "
                "halo-boundary flag",
                flush=True,
            )
        core = box(min_x + tile_x * tile_size_um, min_y + tile_y * tile_size_um,
                   min_x + (tile_x + 1) * tile_size_um, min_y + (tile_y + 1) * tile_size_um)
        halo = core.buffer(halo_um, cap_style="square", join_style="mitre")
        reach = math.ceil(halo_um / tile_size_um)
        nearby = [(y, x) for y in range(tile_y - reach, tile_y + reach + 1)
                  for x in range(tile_x - reach, tile_x + reach + 1)]
        frame = _read_tiles(stage_dir, nearby)
        xy = frame[["x_um", "y_um"]].to_numpy(float)
        local = shapely.covers(halo, shapely.points(xy))
        frame, xy = frame.loc[local].reset_index(drop=True), xy[local]
        # Tile ownership comes from staging, so points exactly on a core border
        # can never be emitted by both adjacent tiles.
        in_core = (
            (frame["_tile_y"].to_numpy() == tile_y)
            & (frame["_tile_x"].to_numpy() == tile_x)
        )
        core_indices = np.flatnonzero(in_core)
        if not len(core_indices):
            marker.write_text("{}\n"); continue
        triangulation = Delaunay(xy)
        pairs, distances = _neighbor_data(xy, triangulation.simplices, max_distance_um)
        counts = np.bincount(pairs.ravel(), minlength=len(xy))
        sums = np.bincount(pairs.ravel(), weights=np.repeat(distances, 2), minlength=len(xy))
        nearest = np.full(len(xy), np.inf)
        np.minimum.at(nearest, pairs[:, 0], distances); np.minimum.at(nearest, pairs[:, 1], distances)
        cells = _voronoi(xy, shapely.intersection(capture_area, halo))
        # A core cell touching the artificial halo boundary means that points
        # beyond the fixed halo could still alter its provisional metrics.
        artificial_boundary = halo.boundary.difference(capture_area.boundary.buffer(1e-8))
        halo_boundary_affected = shapely.intersects(
            cells[core_indices], artificial_boundary
        )
        affected_count = int(np.sum(halo_boundary_affected))
        if affected_count:
            print(
                f"  Warning: tile ({tile_y}, {tile_x}) has "
                f"{affected_count:,} core spots whose Voronoi cells reach the "
                f"{halo_um:g} µm halo boundary",
                flush=True,
            )
        cell_area = shapely.area(cells)
        spot_geoms = shapely.buffer(shapely.points(xy), radius)
        spot_in_cell = shapely.area(shapely.intersection(spot_geoms, cells))
        out = frame.iloc[core_indices].drop(columns=["_tile_y", "_tile_x"]).copy()
        out["n_neighbors"] = counts[core_indices]
        out["mean_neighbor_distance_um"] = np.divide(sums[core_indices], counts[core_indices], out=np.full(len(core_indices), np.nan), where=counts[core_indices] > 0)
        out["nearest_neighbor_distance_um"] = np.where(np.isfinite(nearest[core_indices]), nearest[core_indices], np.nan)
        out["spot_area_um2"] = math.pi * radius**2
        out["voronoi_area_um2"] = cell_area[core_indices]
        out["equivalent_voronoi_diameter_um"] = 2 * np.sqrt(cell_area[core_indices] / math.pi)
        out["spot_fraction_of_voronoi"] = np.divide(spot_in_cell[core_indices], cell_area[core_indices], out=np.full(len(core_indices), np.nan), where=cell_area[core_indices] > 0)
        core_xy = xy[core_indices]
        out["total_counts"] = 0.0
        if umi_by_id is not None:
            mapping_ids = out["spot_id"].astype(str).str.slice(0, 30)
            out["total_counts"] = (
                mapping_ids.map(umi_by_id).fillna(0).to_numpy(float)
            )
        elif umi_tree is not None:
            umi_distance, umi_index = umi_tree.query(
                core_xy, distance_upper_bound=umi_mapping_tolerance_um
            )
            matched = np.isfinite(umi_distance)
            out.loc[matched, "total_counts"] = umi_counts[umi_index[matched]]
        out["in_tissue"] = (
            shapely.covers(tissue_area, shapely.points(core_xy))
            if tissue_area is not None else False
        )
        out["halo_boundary_affected"] = halo_boundary_affected
        _write_partition(out, metrics_dir, tile_y, tile_x, 0)
        vg = out.copy(); vg["geometry"] = shapely.to_wkb(cells[core_indices])
        _write_partition(vg, voronoi_dir, tile_y, tile_x, 0)
        ids = frame["spot_id"].to_numpy()
        reverse = ids[pairs[:, 1]] < ids[pairs[:, 0]]
        oriented = pairs.copy(); oriented[reverse] = oriented[reverse, ::-1]
        owner = in_core[oriented[:, 0]]
        edge_out = pd.DataFrame({"source_id": ids[oriented[owner, 0]], "target_id": ids[oriented[owner, 1]], "distance_um": distances[owner],
                                 "source_x_um": xy[oriented[owner, 0], 0], "source_y_um": xy[oriented[owner, 0], 1],
                                 "target_x_um": xy[oriented[owner, 1], 0], "target_y_um": xy[oriented[owner, 1], 1]})
        _write_partition(edge_out, edges_dir, tile_y, tile_x, 0)
        marker.write_text(json.dumps({
            "spots": len(out),
            "edges": len(edge_out),
            "halo_boundary_affected_spots": affected_count,
        }) + "\n")
        print(f"  Tile {number:,}/{len(keys):,} ({tile_y}, {tile_x}): {len(out):,} spots", flush=True)
    return metrics_dir, voronoi_dir, edges_dir


def remap_novast_umi_counts(metrics_dir, voronoi_dir, umi_observations, output_dir):
    marker = output_dir / "novast_umi_30bp_mapping.json"
    if marker.exists():
        print(f"  Reusing existing 30-bp Nova-ST UMI mapping: {marker}", flush=True)
        return

    counts_by_id = (
        umi_observations.assign(
            spot_id=umi_observations["spot_id"].astype(str).str.slice(0, 30)
        )
        .groupby("spot_id")["total_counts"]
        .sum()
    )
    matched_spots = 0
    nonzero_spots = 0
    for root in (metrics_dir, voronoi_dir):
        paths = sorted(root.glob("tile_y=*/tile_x=*/*.parquet"))
        print(f"  Remapping Nova-ST UMI counts in {root.name}", flush=True)
        for number, path in enumerate(paths, 1):
            table = pq.ParquetFile(path).read()
            frame = table.to_pandas()
            mapping_ids = frame["spot_id"].astype(str).str.slice(0, 30)
            matched = mapping_ids.isin(counts_by_id.index)
            counts = mapping_ids.map(counts_by_id).fillna(0).to_numpy(float)
            frame["total_counts"] = counts
            temporary = path.with_suffix(".parquet.tmp")
            pq.write_table(
                pa.Table.from_pandas(frame, preserve_index=False),
                temporary,
            )
            temporary.replace(path)
            if root == metrics_dir:
                matched_spots += int(matched.sum())
                nonzero_spots += int(np.count_nonzero(counts))
            print(f"    File {number:,}/{len(paths):,}", end="\r", flush=True)
        print(flush=True)

    marker.write_text(json.dumps({
        "identifier_bases": 30,
        "observed_identifiers": int(len(counts_by_id)),
        "matched_whitelist_spots": matched_spots,
        "nonzero_whitelist_spots": nonzero_spots,
    }, indent=2) + "\n")
    print(
        f"  Nova-ST UMI repair complete: {matched_spots:,} identifiers matched; "
        f"{nonzero_spots:,} spots have nonzero counts",
        flush=True,
    )


def summarize_dataset(metrics_dir, capture_area, output_dir):
    dataset = ds.dataset(metrics_dir, format="parquet", partitioning="hive")
    n = 0; affected = 0; covered = 0.0; sums = {column: 0.0 for column in METRIC_COLUMNS}; finite = {column: 0 for column in METRIC_COLUMNS}
    for batch in dataset.scanner(columns=["spot_area_um2", "halo_boundary_affected", *METRIC_COLUMNS], batch_size=250_000).to_batches():
        frame = batch.to_pandas(); n += len(frame)
        affected += int(frame["halo_boundary_affected"].sum())
        covered += float(np.nansum(frame["voronoi_area_um2"] * frame["spot_fraction_of_voronoi"]))
        for column in METRIC_COLUMNS:
            values = frame[column].to_numpy(float); good = np.isfinite(values)
            sums[column] += float(values[good].sum()); finite[column] += int(good.sum())
    summary = {"n_spots": n, "capture_area_um2": capture_area.area,
               "initial_coverage_percent": 100 * covered / capture_area.area,
               "halo_boundary_affected_spots": affected,
               "halo_boundary_affected_fraction": affected / n if n else None,
               "metric_means": {c: sums[c] / finite[c] if finite[c] else None for c in METRIC_COLUMNS}}
    (output_dir / "tiled_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def streamed_coverage(voronoi_dir, capture_area, output_dir, technology,
                      distances_um, saturation_percent=99.9):
    dataset = ds.dataset(voronoi_dir, format="parquet", partitioning="hive")
    radius = {"Nova-ST": .15, "Stereo-Seq": .11}[technology]
    results = []
    for distance in np.asarray(distances_um, float):
        covered = 0.0
        for batch in dataset.scanner(columns=["x_um", "y_um", "geometry"], batch_size=100_000).to_batches():
            frame = batch.to_pandas()
            cells = shapely.from_wkb(frame["geometry"].to_numpy())
            expanded = shapely.buffer(
                shapely.points(frame["x_um"].to_numpy(), frame["y_um"].to_numpy()),
                radius + float(distance),
            )
            covered += float(np.sum(shapely.area(shapely.intersection(expanded, cells))))
        percent = 100 * covered / capture_area.area
        results.append((distance, min(percent, 100.0)))
        print(f"  Coverage at {distance:g} µm: {percent:.3f}%", flush=True)
        if percent >= saturation_percent:
            remaining = np.asarray(distances_um, float)[len(results):]
            results[-1] = (distance, 100.0)
            results.extend((value, 100.0) for value in remaining)
            break
    coverage = pd.DataFrame(results, columns=["distance_from_spot_edge_um", "coverage_percent"])
    coverage.to_csv(output_dir / "spatial_coverage_curve.csv", index=False)
    import matplotlib.pyplot as plt
    figure, axis = plt.subplots(figsize=(6, 4.5))
    axis.plot(coverage.iloc[:, 0], coverage.iloc[:, 1])
    axis.set(xlabel="Distance from spot edge (µm)", ylabel="Capture area covered (%)", ylim=(0, 100))
    axis.grid(alpha=.25); figure.tight_layout()
    figure.savefig(output_dir / "spatial_coverage_curve.png", dpi=300, bbox_inches="tight")
    plt.close(figure)
    return coverage


def streaming_distribution_plots(metrics_dir, output_dir, sample_size=1_000_000, dpi=300):
    dataset = ds.dataset(metrics_dir, format="parquet", partitioning="hive")
    columns = ["n_neighbors", "mean_neighbor_distance_um", "voronoi_area_um2",
               "equivalent_voronoi_diameter_um", "spot_fraction_of_voronoi",
               "total_counts"]
    samples = {column: [] for column in columns}; remaining = sample_size
    # The staged order is deterministic. A uniform stride within every batch
    # avoids retaining the first spatial tiles only.
    total_rows = dataset.count_rows(); stride = max(1, math.ceil(total_rows / sample_size))
    offset = 0
    for batch in dataset.scanner(columns=columns, batch_size=250_000).to_batches():
        frame = batch.to_pandas(); take = np.arange((-offset) % stride, len(frame), stride)
        for column in columns:
            samples[column].append(frame[column].to_numpy()[take])
        offset += len(frame)
    import matplotlib.pyplot as plt
    labels = {
        "n_neighbors": "Delaunay neighbors", "mean_neighbor_distance_um": "Mean neighbor distance (µm)",
        "voronoi_area_um2": "Voronoi area (µm²)",
        "equivalent_voronoi_diameter_um": "Equivalent Voronoi diameter (µm)",
        "spot_fraction_of_voronoi": "Spot fraction of Voronoi cell",
        "total_counts": "Recovered UMI counts",
    }
    arrays = {}
    for column in columns:
        values = np.concatenate(samples[column]) if samples[column] else np.empty(0)
        values = values[np.isfinite(values)]; arrays[column] = values
        figure, axis = plt.subplots(figsize=(6, 4.5))
        if np.isclose(values.min(), values.max()):
            bins = 1
        elif np.all(values == np.floor(values)) and values.max() - values.min() <= 100:
            bins = np.arange(values.min() - .5, values.max() + 1.5)
        else:
            bins = 50
        axis.hist(values, bins=bins, color="#4c91c6", edgecolor="#17324d", linewidth=.5)
        axis.axvline(values.mean(), color="#d1495b", linestyle="--", label=f"Mean: {values.mean():.3g}")
        axis.axvline(np.median(values), color="#2a9d8f", linestyle=":", label=f"Median: {np.median(values):.3g}")
        axis.set(xlabel=labels[column], ylabel="Sampled spots"); axis.legend(frameon=False); figure.tight_layout()
        figure.savefig(output_dir / f"{column}_histogram.png", dpi=dpi, bbox_inches="tight"); plt.close(figure)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, column in zip(axes, ["equivalent_voronoi_diameter_um", "mean_neighbor_distance_um"]):
        values = np.sort(arrays[column]); axis.plot(values, np.arange(1, len(values) + 1) / len(values))
        axis.set(xlabel=labels[column], ylabel="Cumulative fraction"); axis.grid(alpha=.25)
    figure.tight_layout(); figure.savefig(output_dir / "distance_cdfs.png", dpi=dpi, bbox_inches="tight"); plt.close(figure)
    # UMI CDFs deliberately include zeros: they represent usable whitelist
    # positions from which no transcript molecule was recovered.
    for tissue_only, name, title in (
        (False, "umi_count_cdf_capture_area.png", "Complete capture area"),
        (True, "umi_count_cdf_tissue.png", "Tissue area"),
    ):
        sampled = []
        scanner_columns = ["total_counts", "in_tissue"]
        offset = 0
        for batch in dataset.scanner(columns=scanner_columns, batch_size=250_000).to_batches():
            frame = batch.to_pandas()
            if tissue_only:
                frame = frame.loc[frame["in_tissue"]]
            take = np.arange((-offset) % stride, len(frame), stride)
            sampled.append(frame["total_counts"].to_numpy(float)[take])
            offset += len(frame)
        values = np.sort(np.concatenate(sampled)) if sampled else np.empty(0)
        if not len(values):
            continue
        figure, axis = plt.subplots(figsize=(6, 4.5))
        axis.step(values, np.arange(1, len(values) + 1) / len(values), where="post")
        axis.set(xlabel="Recovered UMI counts per spot", ylabel="Cumulative fraction",
                 ylim=(0, 1), title=title)
        axis.grid(alpha=.25); figure.tight_layout()
        figure.savefig(output_dir / name, dpi=dpi, bbox_inches="tight"); plt.close(figure)


def streaming_binned_maps(metrics_dir, capture_area, tissue_area, output_dir,
                          bin_size_um, dpi=300):
    import matplotlib.pyplot as plt
    min_x, min_y, max_x, max_y = capture_area.bounds
    nx = max(1, math.ceil((max_x - min_x) / bin_size_um)); ny = max(1, math.ceil((max_y - min_y) / bin_size_um))
    x_edges = min_x + np.arange(nx + 1) * bin_size_um; y_edges = min_y + np.arange(ny + 1) * bin_size_um
    x_edges[-1], y_edges[-1] = max_x, max_y
    columns = ["n_neighbors", "equivalent_voronoi_diameter_um", "mean_neighbor_distance_um",
               "voronoi_area_um2", "spot_fraction_of_voronoi", "total_counts"]
    counts = np.zeros(nx * ny, np.int64); sums = {c: np.zeros(nx * ny) for c in columns}; observations = {c: np.zeros(nx * ny, np.int64) for c in columns}
    dataset = ds.dataset(metrics_dir, format="parquet", partitioning="hive")
    for batch in dataset.scanner(columns=["x_um", "y_um", *columns], batch_size=250_000).to_batches():
        frame = batch.to_pandas(); xb = np.clip(((frame.x_um - min_x) / bin_size_um).astype(int), 0, nx - 1); yb = np.clip(((frame.y_um - min_y) / bin_size_um).astype(int), 0, ny - 1); flat = yb * nx + xb
        counts += np.bincount(flat, minlength=nx * ny)
        for column in columns:
            values = frame[column].to_numpy(float); good = np.isfinite(values)
            sums[column] += np.bincount(flat[good], weights=values[good], minlength=nx * ny)
            observations[column] += np.bincount(flat[good], minlength=nx * ny)
    polygons = np.asarray([box(x_edges[x], y_edges[y], x_edges[x + 1], y_edges[y + 1]) for y in range(ny) for x in range(nx)], dtype=object)
    bin_area = shapely.area(shapely.intersection(polygons, capture_area))
    density = np.divide(counts * 1_000_000, bin_area, out=np.full(nx * ny, np.nan), where=bin_area > 0)
    names = {"n_neighbors": "average_number_of_neighbors", "equivalent_voronoi_diameter_um": "average_equivalent_voronoi_diameter",
             "mean_neighbor_distance_um": "average_neighbor_distance", "voronoi_area_um2": "average_voronoi_area",
             "spot_fraction_of_voronoi": "average_spot_fraction_of_voronoi",
             "total_counts": "average_umi_counts"}
    panels = [(density, "spot_density", "Spot density")]
    panels += [(np.divide(sums[c], observations[c], out=np.full(nx * ny, np.nan), where=observations[c] > 0), names[c], c.replace("_", " ")) for c in columns]
    for values, name, title in panels:
        figure, axis = plt.subplots(figsize=(7, 7)); image = axis.pcolormesh(x_edges, y_edges, values.reshape(ny, nx), shading="flat", cmap="viridis")
        if tissue_area is not None:
            for polygon in shapely.get_parts(tissue_area):
                if polygon.geom_type == "Polygon":
                    xy = np.asarray(polygon.exterior.coords)
                    axis.plot(xy[:, 0], xy[:, 1], color="red", linewidth=1.0, zorder=4)
                    for ring in polygon.interiors:
                        xy = np.asarray(ring.coords)
                        axis.plot(xy[:, 0], xy[:, 1], color="red", linewidth=1.0, zorder=4)
        figure.colorbar(image, ax=axis, shrink=.82); axis.set(title=title, xlabel="x (µm)", ylabel="y (µm)", aspect="equal", ylim=(max_y, min_y)); figure.tight_layout()
        figure.savefig(output_dir / f"{name}_{bin_size_um:g}um_bins.png", dpi=dpi, bbox_inches="tight"); plt.close(figure)


def streaming_spatial_crops(voronoi_dir, edges_dir, capture_area, tissue_area, output_dir,
                            crop_size_um, technology, dpi=300):
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection, LineCollection
    from matplotlib.colors import Normalize
    from shapely.geometry import Point
    from shapely.ops import nearest_points
    half = crop_size_um / 2; center = capture_area.centroid
    min_x, min_y, max_x, _ = capture_area.bounds
    top_left = nearest_points(capture_area.boundary, Point(min_x, min_y))[0]
    top_right = nearest_points(capture_area.boundary, Point(max_x, min_y))[0]
    crops = {
        "center": [
            ("center", (center.x-half, center.y-half,
                        center.x+half, center.y+half)),
            ("center", (center.x+100-half, center.y+100-half,
                        center.x+100+half, center.y+100+half)),
        ],
        "top_left": [
            ("top_left", (top_left.x-.25*crop_size_um,
                          top_left.y-.25*crop_size_um,
                          top_left.x+.75*crop_size_um,
                          top_left.y+.75*crop_size_um)),
            ("top_right", (top_right.x-.75*crop_size_um,
                           top_right.y-.25*crop_size_um,
                           top_right.x+.25*crop_size_um,
                           top_right.y+.75*crop_size_um)),
        ],
    }
    metrics_ds = ds.dataset(voronoi_dir, format="parquet", partitioning="hive"); edges_ds = ds.dataset(edges_dir, format="parquet", partitioning="hive")
    labels = {"n_neighbors": "Delaunay neighbors", "mean_neighbor_distance_um": "Mean neighbor distance (µm)",
              "voronoi_area_um2": "Voronoi area (µm²)", "equivalent_voronoi_diameter_um": "Equivalent Voronoi diameter (µm)",
              "spot_fraction_of_voronoi": "Spot fraction of Voronoi cell",
              "total_counts": "Recovered UMI counts"}
    radius = {"Nova-ST": .15, "Stereo-Seq": .11}[technology]
    for _, candidates in crops.items():
        frame = None
        for crop_name, (x0, y0, x1, y1) in candidates:
            filt = (ds.field("x_um") >= x0) & (ds.field("x_um") <= x1) & (ds.field("y_um") >= y0) & (ds.field("y_um") <= y1)
            frame = metrics_ds.to_table(filter=filt).to_pandas()
            if not frame.empty:
                break
        if frame.empty: continue
        cells = shapely.from_wkb(frame["geometry"].to_numpy()); segments = [np.asarray(poly.exterior.coords) for cell in cells for poly in shapely.get_parts(cell) if poly.geom_type == "Polygon"]
        edge_filter = (((ds.field("source_x_um") >= x0) & (ds.field("source_x_um") <= x1) & (ds.field("source_y_um") >= y0) & (ds.field("source_y_um") <= y1)) |
                       ((ds.field("target_x_um") >= x0) & (ds.field("target_x_um") <= x1) & (ds.field("target_y_um") >= y0) & (ds.field("target_y_um") <= y1)))
        edge_frame = edges_ds.to_table(filter=edge_filter).to_pandas()
        edge_segments = np.stack((edge_frame[["source_x_um", "source_y_um"]].to_numpy(), edge_frame[["target_x_um", "target_y_um"]].to_numpy()), axis=1) if len(edge_frame) else np.empty((0, 2, 2))
        xy = frame[["x_um", "y_um"]].to_numpy()
        for column, label in labels.items():
            values = frame[column].to_numpy(float); finite = values[np.isfinite(values)]
            if column == "total_counts":
                vmin, vmax = 0.0, np.nanpercentile(finite, 99)
            else:
                vmin, vmax = np.nanpercentile(finite, [1, 99])
            vmax = vmax if vmax > vmin else vmin + 1
            cmap = plt.get_cmap("viridis").copy()
            display_values = values.copy()
            if column == "total_counts":
                cmap.set_under("black")
                display_values[display_values == 0] = -max(1.0, float(vmax)) * 1e-6
            figure, axis = plt.subplots(figsize=(7, 7)); axis.add_collection(LineCollection(segments, colors=".35", linewidths=.35, alpha=.7, zorder=1))
            spots = EllipseCollection(widths=2*radius, heights=2*radius, angles=0, units="xy", offsets=xy, transOffset=axis.transData,
                                      cmap=cmap, norm=Normalize(vmin, vmax), edgecolors="black", linewidths=.15, zorder=3)
            spots.set_array(display_values); axis.add_collection(spots)
            if len(edge_segments): axis.add_collection(LineCollection(edge_segments, colors="#d1495b", linewidths=.5, alpha=.8, zorder=5))
            if tissue_area is not None:
                tissue_crop = shapely.intersection(tissue_area, box(x0, y0, x1, y1))
                for polygon in shapely.get_parts(tissue_crop):
                    if polygon.geom_type == "Polygon":
                        boundary = np.asarray(polygon.exterior.coords)
                        axis.plot(boundary[:, 0], boundary[:, 1], color="red", linewidth=1.0, zorder=6)
                        for ring in polygon.interiors:
                            boundary = np.asarray(ring.coords)
                            axis.plot(boundary[:, 0], boundary[:, 1], color="red", linewidth=1.0, zorder=6)
            figure.colorbar(spots, ax=axis, shrink=.8, label=label); axis.set(xlim=(x0,x1), ylim=(y1,y0), xlabel="x (µm)", ylabel="y (µm)", title=f"{label} — {crop_name} crop", aspect="equal")
            figure.tight_layout(); figure.savefig(output_dir / f"{column}_spatial_{crop_name}.png", dpi=dpi, bbox_inches="tight"); plt.close(figure)


def run_tiled_sample(sample, capture_area, output_dir, max_neighbors, max_distance_um,
                     processing_tile_size_um=500.0, processing_halo_um=25.0,
                     coverage_max_um=50.0, coverage_steps=51,
                     spatial_bin_size_um=200.0, crop_size_um=120.0, dpi=300,
                     saturation_percent=99.9, tissue_area=None,
                     umi_observations=None, umi_mapping_tolerance_um=0.05):
    if max_neighbors is not None:
        raise ValueError("The scalable tiled path currently requires max_neighbors=None")
    if max_distance_um is None:
        raise ValueError("The scalable tiled path requires a finite max_distance_um")
    output_dir.mkdir(parents=True, exist_ok=True)
    configuration = {
        "dataset_id": sample.get("dataset_id"),
        "technology": sample["technology"],
        "max_neighbors": max_neighbors,
        "max_distance_um": max_distance_um,
        "processing_tile_size_um": processing_tile_size_um,
        "processing_halo_um": processing_halo_um,
        "umi_mapping_method": (
            "exact_identifier" if sample["technology"] == "Nova-ST"
            else "coordinate"
        ),
    }
    configuration_file = output_dir / "tiled_configuration.json"
    if configuration_file.exists():
        previous = json.loads(configuration_file.read_text())
        if previous != configuration:
            raise ValueError(
                "The existing tiled output was created with different processing "
                f"arguments. Use a new output directory or restore: {previous}"
            )
    else:
        configuration_file.write_text(json.dumps(configuration, indent=2) + "\n")
    stage = output_dir / "position_tiles"
    start = perf_counter()
    stage_summary = stage_positions(sample, capture_area, stage, processing_tile_size_um)
    had_existing_tiles = any((output_dir / "tile_markers").glob("*.json"))
    metrics, voronoi, edges = process_tiles(sample, capture_area, stage, output_dir,
                                  processing_tile_size_um, processing_halo_um, max_distance_um,
                                  tissue_area=tissue_area,
                                  umi_observations=umi_observations,
                                  umi_mapping_tolerance_um=umi_mapping_tolerance_um)
    if sample["technology"] == "Nova-ST" and umi_observations is not None:
        umi_marker = output_dir / "novast_umi_30bp_mapping.json"
        if not umi_marker.exists():
            if had_existing_tiles:
                remap_novast_umi_counts(
                    metrics, voronoi, umi_observations, output_dir
                )
            else:
                umi_marker.write_text(json.dumps({
                    "identifier_bases": 30,
                    "mapping_created_during_tile_processing": True,
                }, indent=2) + "\n")
    summary = summarize_dataset(metrics, capture_area, output_dir)
    coverage_file = output_dir / "spatial_coverage_curve.csv"
    if coverage_file.exists():
        print(
            f"  Reusing existing spatial coverage: {coverage_file}",
            flush=True,
        )
    else:
        streamed_coverage(voronoi, capture_area, output_dir, sample["technology"],
                          np.linspace(0, coverage_max_um, coverage_steps),
                          saturation_percent=saturation_percent)
    streaming_distribution_plots(metrics, output_dir, dpi=dpi)
    streaming_binned_maps(metrics, capture_area, tissue_area, output_dir, spatial_bin_size_um, dpi=dpi)
    streaming_spatial_crops(voronoi, edges, capture_area, tissue_area, output_dir, crop_size_um, sample["technology"], dpi=dpi)
    summary.update(stage_summary=stage_summary, processing_tile_size_um=processing_tile_size_um,
                   processing_halo_um=processing_halo_um, elapsed_seconds=perf_counter() - start)
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
