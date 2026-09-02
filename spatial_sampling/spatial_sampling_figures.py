from pathlib import Path
import math

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection, LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, to_rgba
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, MultipleLocator, NullLocator, StrMethodFormatter
import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy.ndimage import distance_transform_edt
import shapely
from matplotlib.ticker import LogLocator

from spatial_sampling_metrics import SPOT_CONFIG, get_capture_area, get_tissue_area


TECHNOLOGY_COLORS = {
    "BMKMANU S1000": "#000000", "BMKMANU S3000": "#CC79A7",
    "Curio Seeker": "#9B9B9B", "Nova-ST": "#D55E00",
    "Stereo-Seq": "#0072B2", "Visium (FFPE)": "#E69F00",
    "Visium HD (FFPE)": "#009E73", "Visium HD (FF)": "#F0E442",
    "Visium HD 3'": "#56B4E9",
}

MOUSE_BRAIN_CROP_ORIGINS_PX = {
    "Curio Seeker": (2500, 2500),
    "Nova-ST": (6090, 2960), "BMKMANU S1000": (9880, 10300),
    "BMKMANU S3000": (10300, 10100), "Stereo-Seq": (11885, 11875),
}


def _result_dir(sample, results_root):
    name = sample["technology"] + (f"_{sample['replicate']}" if sample.get("replicate") else "")
    return Path(results_root) / sample["sample_type"] / name


def _dataset(path):
    return ds.dataset(path, format="parquet", partitioning="hive")


def _read_metrics(sample, results_root, columns, bounds=None, sample_size=None):
    root = _result_dir(sample, results_root)
    path = root / "spot_metrics"
    if not path.exists():
        path = root / "spot_metrics.parquet"
    filt = None
    if bounds is not None:
        x0, y0, x1, y1 = bounds
        filt = ((ds.field("x_um") >= x0) & (ds.field("x_um") <= x1) &
                (ds.field("y_um") >= y0) & (ds.field("y_um") <= y1))
    dataset = _dataset(path)
    if sample_size is None:
        return dataset.to_table(columns=columns, filter=filt).to_pandas()
    total = dataset.count_rows(filter=filt)
    stride = max(1, math.ceil(total / sample_size))
    frames=[]; offset=0
    for batch in dataset.scanner(columns=columns, filter=filt, batch_size=250_000).to_batches():
        frame=batch.to_pandas(); take=np.arange((-offset) % stride,len(frame),stride)
        frames.append(frame.iloc[take]); offset += len(frame)
    return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=columns)


def _outline(axis, geometry, color, linewidth=1.2, zorder=6):
    for polygon in shapely.get_parts(geometry):
        if polygon.geom_type != "Polygon":
            continue
        xy = np.asarray(polygon.exterior.coords)
        axis.plot(xy[:, 0], xy[:, 1], color=color, linewidth=linewidth,
                  zorder=zorder)
        for ring in polygon.interiors:
            xy = np.asarray(ring.coords)
            axis.plot(xy[:, 0], xy[:, 1], color=color, linewidth=linewidth,
                      zorder=zorder)


def _scale_bar(axis, length_um=None, color="black", linewidth=2.5):
    x0, x1 = axis.get_xlim(); y0, y1 = axis.get_ylim()
    width = abs(x1 - x0)
    if length_um is None:
        length_um = width / 4
    x = min(x0, x1) + .06 * width
    y = max(y0, y1) - .07 * abs(y1 - y0)
    axis.plot(
        [x, x + length_um], [y, y], color=color, linewidth=linewidth,
        solid_capstyle="butt", zorder=100, clip_on=False,
    )


def _sample_geometry(sample):
    import spatialdata as sd
    sdata = sd.read_zarr(sample["sdata_path"])
    capture = get_capture_area(sdata, sample)
    return sdata, capture, get_tissue_area(sdata, sample, capture)


def export_coverage_map(sample, results_root, output_file, pixel_size_um=2.0,
                        maximum_distance_um=50.0, dpi=300,
                        microns_per_inch=2000.0, cividis_start=0.15,
                        black_and_white=False,
                        sample_geometry=None):
    _, capture, tissue = (
        _sample_geometry(sample) if sample_geometry is None else sample_geometry
    )
    tissue_in_capture = (
        shapely.intersection(capture, tissue) if tissue is not None else capture
    )
    minx, miny, maxx, maxy = capture.bounds
    nx = math.ceil((maxx-minx)/pixel_size_um) + 1
    ny = math.ceil((maxy-miny)/pixel_size_um) + 1
    occupied = np.zeros((ny, nx), dtype=bool)
    dataset = _dataset((_result_dir(sample, results_root) / "spot_metrics")
                       if (_result_dir(sample, results_root) / "spot_metrics").exists()
                       else (_result_dir(sample, results_root) / "spot_metrics.parquet"))
    for batch in dataset.scanner(columns=["x_um", "y_um"], batch_size=500_000).to_batches():
        frame = batch.to_pandas()
        ix = np.clip(np.rint((frame.x_um-minx)/pixel_size_um).astype(int), 0, nx-1)
        iy = np.clip(np.rint((maxy-frame.y_um)/pixel_size_um).astype(int), 0, ny-1)
        occupied[iy, ix] = True
    radius = SPOT_CONFIG[sample["technology"]]["size_um"] / 2
    distance = np.maximum(distance_transform_edt(~occupied) * pixel_size_um - radius, 0)
    # Rasterization scales with polygon complexity and avoids millions of
    # individual point-in-polygon tests for detailed tissue annotations.
    inside = rasterize(
        [(capture, 1)],
        out_shape=(ny, nx),
        transform=from_origin(minx, maxy, pixel_size_um, pixel_size_um),
        fill=0,
        dtype="uint8",
    ).astype(bool)
    distance[~inside] = np.nan
    if black_and_white:
        display_distance = np.where(np.isnan(distance), np.nan, distance >= 2.0)
        cmap = ListedColormap(["black", "white"])
        vmin, vmax = 0, 1
    else:
        display_distance = distance.copy()
        display_distance[display_distance <= 2] = -max(1.0, maximum_distance_um) * 1e-6
        if not 0 <= cividis_start < 1:
            raise ValueError("cividis_start must be at least 0 and less than 1")
        base_cmap = plt.get_cmap("cividis")
        cmap = LinearSegmentedColormap.from_list(
            "light_cividis",
            base_cmap(np.linspace(cividis_start, 1.0, 256)),
        )
        cmap.set_under("black")
        vmin, vmax = 0, maximum_distance_um
    cmap.set_bad((0, 0, 0, 0))

    # The map axes use a fixed physical scale. At the default 2,000 µm/inch,
    # every 1 mm scale bar is exactly 0.5 inch in every exported PNG.
    map_width_in = (maxx - minx) / microns_per_inch
    map_height_in = (maxy - miny) / microns_per_inch
    colorbar_width_in, colorbar_gap_in, colorbar_margin_in = 0.18, 0.12, 0.85
    figure_width = map_width_in if black_and_white else (
        map_width_in + colorbar_gap_in + colorbar_width_in + colorbar_margin_in
    )
    fig = plt.figure(figsize=(figure_width, map_height_in))
    ax = fig.add_axes([0, 0, map_width_in / figure_width, 1])
    if not black_and_white:
        colorbar_axis = fig.add_axes([
            (map_width_in + colorbar_gap_in) / figure_width,
            0.08,
            colorbar_width_in / figure_width,
            0.84,
        ])
    image = ax.imshow(
        display_distance, extent=(minx, maxx, miny, maxy), origin="upper",
        cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )
    if tissue is not None and sample["technology"] != "Curio Seeker":
        _outline(ax, tissue_in_capture, "red")
    _outline(ax, capture, "cornflowerblue", linewidth=2, zorder=7)
    ax.set(xlim=(minx, maxx), ylim=(maxy, miny), aspect="equal")
    ax.set_axis_off()
    _scale_bar(ax, length_um=1000.0, color="white")
    if not black_and_white:
        fig.colorbar(
            image, cax=colorbar_axis, label="Distance to closest spot edge (µm)"
        )
    fig.savefig(output_file, dpi=dpi, transparent=True)
    plt.close(fig)


def _paper_crop(sample, capture, crop_size_um, sdata):
    if "paper_crop_bounds_um" in sample:
        return tuple(sample["paper_crop_bounds_um"])
    origin = MOUSE_BRAIN_CROP_ORIGINS_PX.get(sample["technology"])
    if sample["sample_type"] == "Mouse brain" and origin is not None:
        x0, y0 = np.asarray(origin) * sample["microns_per_pixel"]
    elif sample["sample_type"] == "Mouse brain" and sample["technology"] == "Curio Seeker":
        shape = sdata.labels[sample.get("labels_layer", f"{sample['dataset_id']}_labels")].shape
        x0, y0 = shape[0] / 2 * sample["microns_per_pixel"], shape[1] / 2 * sample["microns_per_pixel"]
    elif sample["sample_type"] == "Mouse brain" and sample.get("image_layer") in sdata.images:
        shape = sdata.images[sample["image_layer"]]["scale0"]["image"].shape
        x0, y0 = shape[2] / 2 * sample["microns_per_pixel"], shape[1] / 2 * sample["microns_per_pixel"]
    else:
        center = capture.centroid; x0 = center.x-crop_size_um/2; y0 = center.y-crop_size_um/2
    return x0, y0, x0+crop_size_um, y0+crop_size_um


def export_umi_spatial_crop(sample, results_root, output_file, crop_size_um=120,
                            dpi=300, color_limits=None, separate_colorbar=True,
                            detail_inset_zoom=10.0,
                            binary_umi=False,
                            sample_geometry=None):
    sdata, capture, tissue = (
        _sample_geometry(sample) if sample_geometry is None else sample_geometry
    )
    bounds = _paper_crop(sample, capture, crop_size_um, sdata)
    root = _result_dir(sample, results_root)
    vor_path = root / "voronoi_cells"
    if not vor_path.exists(): vor_path = root / "voronoi_cells.parquet"
    x0,y0,x1,y1 = bounds
    size = SPOT_CONFIG[sample["technology"]]["size_um"]
    if sample["technology"] == "Curio Seeker":
        size *= .75
    spot_margin = size / 2
    # Include centres just outside the crop when their physical spot footprint
    # overlaps it; the axes clip the glyph at the crop boundary.
    x_filter_min, x_filter_max = x0 - spot_margin, x1 + spot_margin
    y_filter_min, y_filter_max = y0 - spot_margin, y1 + spot_margin
    filt = ((ds.field("x_um") >= x_filter_min) & (ds.field("x_um") <= x_filter_max) &
            (ds.field("y_um") >= y_filter_min) & (ds.field("y_um") <= y_filter_max))
    frame = _dataset(vor_path).to_table(filter=filt).to_pandas()
    if frame.empty:
        raise ValueError(f"The selected paper crop contains no spots: {bounds}")
    geometry = shapely.from_wkb(frame.geometry.to_numpy()) if frame.geometry.dtype == object and isinstance(frame.geometry.iloc[0], bytes) else frame.geometry.to_numpy()
    segments = [np.asarray(p.exterior.coords) for g in geometry for p in shapely.get_parts(g) if p.geom_type == "Polygon"]
    edge_path = root / "delaunay_edges"
    if not edge_path.exists(): edge_path = root / "delaunay_edges.parquet"
    edge_dataset = _dataset(edge_path)
    if {"source_x_um", "source_y_um", "target_x_um", "target_y_um"}.issubset(
        edge_dataset.schema.names
    ):
        # Disk-backed tiled outputs store edge endpoint coordinates directly.
        ef = (((ds.field("source_x_um") >= x_filter_min) & (ds.field("source_x_um") <= x_filter_max) &
               (ds.field("source_y_um") >= y_filter_min) & (ds.field("source_y_um") <= y_filter_max)) |
              ((ds.field("target_x_um") >= x_filter_min) & (ds.field("target_x_um") <= x_filter_max) &
               (ds.field("target_y_um") >= y_filter_min) & (ds.field("target_y_um") <= y_filter_max)))
        edges = edge_dataset.to_table(filter=ef).to_pandas()
    else:
        # Regular outputs store integer row references into spot_metrics.
        edge_indices = edge_dataset.to_table(
            columns=["source_index", "target_index", "distance_um"]
        ).to_pandas()
        metric_path = root / "spot_metrics.parquet"
        coordinates = _dataset(metric_path).to_table(
            columns=["x_um", "y_um"]
        ).to_pandas().to_numpy(float)
        source = edge_indices["source_index"].to_numpy(int)
        target = edge_indices["target_index"].to_numpy(int)
        source_xy, target_xy = coordinates[source], coordinates[target]
        visible = (
            ((source_xy[:, 0] >= x_filter_min) & (source_xy[:, 0] <= x_filter_max) &
             (source_xy[:, 1] >= y_filter_min) & (source_xy[:, 1] <= y_filter_max)) |
            ((target_xy[:, 0] >= x_filter_min) & (target_xy[:, 0] <= x_filter_max) &
             (target_xy[:, 1] >= y_filter_min) & (target_xy[:, 1] <= y_filter_max))
        )
        edges = edge_indices.loc[visible].copy()
        edges["source_x_um"] = source_xy[visible, 0]
        edges["source_y_um"] = source_xy[visible, 1]
        edges["target_x_um"] = target_xy[visible, 0]
        edges["target_y_um"] = target_xy[visible, 1]
    values = frame.total_counts.to_numpy(float)
    if binary_umi:
        display_values = (values >= 1).astype(float)
        norm = Normalize(0, 1)
        cmap = ListedColormap(["#E2F4FA", "#ECD8A2"])
    else:
        if color_limits is None:
            color_limits = (0.0, float(np.nanpercentile(values, 99)))
        else:
            color_limits = (0.0, color_limits[1])
        if color_limits[1] <= color_limits[0]:
            color_limits = (color_limits[0], color_limits[0]+1)
        norm = Normalize(*color_limits)
        cmap = plt.get_cmap("viridis").copy(); cmap.set_under("black")
        display_values = values.copy()
        display_values[display_values == 0] = (
            -max(1.0, float(color_limits[1])) * 1e-6
        )
    xy = frame[["x_um","y_um"]].to_numpy()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    if binary_umi:
        has_umi = values >= 1
        face_colors = np.tile(to_rgba("#E2F4FA"), (len(values), 1))
        edge_colors = np.tile(to_rgba("#57C0E1"), (len(values), 1))
        face_colors[has_umi] = to_rgba("#ECD8A2")
        edge_colors[has_umi] = to_rgba("#F5A509")
        face_colors[:, 3] = .9
    else:
        mapped_colors = cmap(norm(display_values))
        face_colors = mapped_colors.copy(); face_colors[:, 3] = .9
        edge_colors = mapped_colors.copy(); edge_colors[:, 3] = 1.0
    lines = (
        np.stack((edges[["source_x_um","source_y_um"]],
                  edges[["target_x_um","target_y_um"]]), axis=1)
        if len(edges) else []
    )

    def add_spatial_layers(axis, line_scale=1.0, show_geometry=True,
                           neighbors_on_top=False,
                           voronoi_on_top=False,
                           neighbor_linewidth=.80):
        if show_geometry:
            axis.add_collection(LineCollection(
                segments, colors="#292929",
                linewidths=.45 * line_scale,
                zorder=3.5 if voronoi_on_top else 3,
            ))
            if len(lines):
                axis.add_collection(LineCollection(
                    lines, colors="#CC79A7", linewidths=neighbor_linewidth,
                    alpha=1, zorder=4 if neighbors_on_top else 2,
                ))
        if SPOT_CONFIG[sample["technology"]]["shape"] == "circle":
            spots = EllipseCollection(
                size, size, 0, units="xy", offsets=xy,
                transOffset=axis.transData, facecolors=face_colors,
                edgecolors=edge_colors, linewidths=.80 * line_scale, zorder=3,
            )
        else:
            h = size / 2
            spot_polygons = [
                [(x-h,y-h), (x+h,y-h), (x+h,y+h), (x-h,y+h)]
                for x, y in xy
            ]
            spots = PolyCollection(
                spot_polygons, facecolors=face_colors, edgecolors=edge_colors,
                linewidths=.80 * line_scale, zorder=3,
            )
        axis.add_collection(spots)

    technology = sample["technology"]
    use_detail_inset = (
        technology in {"Nova-ST", "Stereo-Seq"}
        or technology.startswith("Visium HD")
    ) and detail_inset_zoom is not None
    add_spatial_layers(ax, show_geometry=not use_detail_inset)
    if tissue is not None and sample["technology"] != "Curio Seeker":
        # Plot the original boundary and let the axes clip it. Intersecting the
        # polygon with the crop would manufacture a red rectangular crop frame.
        _outline(ax, tissue, "red", linewidth=1.0, zorder=5)
    ax.set(xlim=(x0,x1), ylim=(y1,y0), aspect="equal"); ax.set_axis_off()
    _scale_bar(ax, linewidth=7.5)
    if use_detail_inset:
        if detail_inset_zoom <= 1:
            raise ValueError("detail_inset_zoom must be greater than 1 or None")
        center_x, center_y = (x0 + x1) / 2, (y0 + y1) / 2
        inset_width = (x1 - x0) / detail_inset_zoom
        inset_height = (y1 - y0) / detail_inset_zoom
        inset_bounds = (
            center_x - inset_width / 2, center_y - inset_height / 2,
            center_x + inset_width / 2, center_y + inset_height / 2,
        )
        ix0, iy0, ix1, iy1 = inset_bounds
        ax.add_patch(Rectangle(
            (ix0, iy0), inset_width, inset_height, fill=False,
            edgecolor="black", linewidth=1.0, zorder=20,
        ))
        inset_ax = fig.add_axes([.59, .04, .39, .39], facecolor="white")
        # savefig(transparent=True) makes axes patches transparent. A regular
        # artist preserves the white inset background while leaving the rest
        # of the exported panel transparent.
        inset_ax.add_patch(Rectangle(
            (0, 0), 1, 1, transform=inset_ax.transAxes,
            facecolor="white", edgecolor="none", zorder=-100,
            clip_on=False,
        ))
        add_spatial_layers(
            inset_ax, line_scale=1,
            neighbors_on_top=technology.startswith("Visium HD"),
            voronoi_on_top=technology.startswith("Visium HD"),
            neighbor_linewidth=.65,
        )
        inset_ax.set(
            xlim=(ix0, ix1), ylim=(iy1, iy0), aspect="equal",
            xticks=[], yticks=[],
        )
        for spine in inset_ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.2)
    if not separate_colorbar:
        colorbar = fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap), ax=ax,
            label="Recovered UMI counts", shrink=.75,
            ticks=[0, 1] if binary_umi else None,
        )
        if binary_umi:
            colorbar.ax.set_yticklabels(["0", "≥1"])
    fig.savefig(output_file, dpi=dpi, transparent=True, pad_inches=0)
    plt.close(fig)
    if separate_colorbar:
        cbfile = Path(output_file).with_name(Path(output_file).stem + "_colorbar" + Path(output_file).suffix)
        fig, ax = plt.subplots(figsize=(1.2,3.2)); ax.set_axis_off()
        colorbar = fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap), ax=ax,
            label="Recovered UMI counts", fraction=.8,
            ticks=[0, 1] if binary_umi else None,
        )
        if binary_umi:
            colorbar.ax.set_yticklabels(["0", "≥1"])
        fig.savefig(cbfile, dpi=dpi, bbox_inches="tight", transparent=True); plt.close(fig)


def _ecdf(values):
    x=np.sort(np.asarray(values,float)); x=x[np.isfinite(x)]
    return x, np.arange(1,len(x)+1)/len(x)


def _label(sample):
    return sample["technology"] + (f" {sample['replicate']}" if sample.get("replicate") else "")


def _zero_marker_x_positions(samples, spacing=0.35):
    """Place Visium HD zero markers symmetrically around the y-axis."""
    positions = {id(sample): 0.0 for sample in samples}
    hd_samples = [
        sample for sample in samples
        if sample["technology"].startswith("Visium HD")
    ]
    center = (len(hd_samples) - 1) / 2
    for index, sample in enumerate(hd_samples):
        positions[id(sample)] = (index - center) * spacing
    return positions


def _format_percent_axis(axis, upper=103):
    """Format a linear percentage axis with labels every 10 percentage points."""
    axis.set_ylim(0, upper)
    axis.yaxis.set_major_locator(MultipleLocator(10))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))


def _add_log_subdivision_grid(axis, minimum=0.1):
    """Draw vertical lines at 1–9 times every power of ten."""
    maximum = axis.get_xlim()[1]
    if not np.isfinite(maximum) or maximum < minimum:
        return
    first_exponent = math.floor(math.log10(minimum))
    last_exponent = math.ceil(math.log10(maximum))
    for exponent in range(first_exponent, last_exponent + 1):
        decade = 10.0 ** exponent
        for multiplier in range(1, 10):
            value = multiplier * decade
            if minimum <= value <= maximum:
                axis.axvline(
                    value, color="lightgrey",
                    linewidth=1.1 if multiplier == 1 else .45,
                    zorder=0,
                )


def _quantile_summary(values, prefix):
    """Summarize a finite, representative sample of a distribution."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {f"{prefix}_{name}": np.nan for name in (
            "median", "q1", "q3", "iqr", "p95", "p99"
        )}
    q1, median, q3, p95, p99 = np.quantile(
        values, [.25, .50, .75, .95, .99]
    )
    return {
        f"{prefix}_median": median,
        f"{prefix}_q1": q1,
        f"{prefix}_q3": q3,
        f"{prefix}_iqr": q3 - q1,
        f"{prefix}_p95": p95,
        f"{prefix}_p99": p99,
    }


def _coverage_value(coverage, distance_um):
    """Read an exact saved coverage value, including a saturated 100% tail."""
    x = coverage.iloc[:, 0].to_numpy(float)
    y = coverage.iloc[:, 1].to_numpy(float)
    match = np.isclose(x, distance_um, atol=1e-9)
    if np.any(match):
        return float(y[np.flatnonzero(match)[0]])
    if len(x) and x[-1] < distance_um and y[-1] >= 100:
        return 100.0
    return np.nan


def export_sample_summary_statistics(samples, results_root, output_file,
                                     sample_size=1_000_000):
    rows = []
    for sample in samples:
        root = _result_dir(sample, results_root)
        distance_frame = _read_metrics(
            sample, results_root,
            ["mean_neighbor_distance_um", "equivalent_voronoi_diameter_um"],
            sample_size=sample_size,
        )
        coverage = pd.read_csv(root / "spatial_coverage_curve.csv")
        metric_path = root / "spot_metrics"
        if not metric_path.exists():
            metric_path = root / "spot_metrics.parquet"
        metric_dataset = _dataset(metric_path)
        tissue_filter = ds.field("in_tissue") == True
        tissue_count = metric_dataset.count_rows(filter=tissue_filter)
        tissue_sum = 0.0
        tissue_zero = 0
        umi_sample = []
        stride = max(1, math.ceil(tissue_count / sample_size))
        offset = 0
        for batch in metric_dataset.scanner(
                columns=["total_counts"], filter=tissue_filter,
                batch_size=250_000).to_batches():
            values = batch.column(0).to_numpy(zero_copy_only=False)
            tissue_sum += float(np.nansum(values))
            tissue_zero += int(np.count_nonzero(values == 0))
            take = np.arange((-offset) % stride, len(values), stride)
            umi_sample.append(values[take])
            offset += len(values)
        umi_values = (
            np.concatenate(umi_sample) if umi_sample else np.empty(0)
        )
        row = {
            "sample_type": sample["sample_type"],
            "technology": sample["technology"],
            "replicate": sample.get("replicate"),
            "dataset_id": sample["dataset_id"],
            "distance_statistics_sampled_spots": len(distance_frame),
            "tissue_spots": tissue_count,
            "umi_statistics_sampled_spots": len(umi_values),
            "functional_area_percent": _coverage_value(coverage, 0),
            "spatial_coverage_2um_percent": _coverage_value(coverage, 2),
            "spatial_coverage_10um_percent": _coverage_value(coverage, 10),
            "spatial_coverage_20um_percent": _coverage_value(coverage, 20),
            "tissue_umi_sparsity_percent": (
                100 * tissue_zero / tissue_count if tissue_count else np.nan
            ),
            "tissue_umi_mean": (
                tissue_sum / tissue_count if tissue_count else np.nan
            ),
        }
        row.update(_quantile_summary(
            distance_frame["mean_neighbor_distance_um"],
            "mean_neighbor_distance_um",
        ))
        row.update(_quantile_summary(
            distance_frame["equivalent_voronoi_diameter_um"],
            "equivalent_voronoi_diameter_um",
        ))
        umi_summary = _quantile_summary(umi_values, "tissue_umi")
        # Only the median and IQR were requested for UMI counts, but retaining
        # Q1 and Q3 makes the reported interval unambiguous.
        for key in ("tissue_umi_median", "tissue_umi_q1",
                    "tissue_umi_q3", "tissue_umi_iqr"):
            row[key] = umi_summary[key]
        rows.append(row)
    summary = pd.DataFrame(rows)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_file, index=False)
    return summary


def export_tissue_combined_plots(samples, results_root, output_dir, dpi=300):
    output_dir=Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    reference_samples = [
        sample for sample in samples if sample["sample_type"] == "Mouse brain"
    ]
    if not reference_samples:
        raise ValueError(
            "Mouse brain samples are required to define Panel D/E x-axis limits"
        )

    reference_distance_max = 0.0
    reference_umi_max = 0.0
    for sample in reference_samples:
        reference_distances = _read_metrics(
            sample, results_root,
            ["mean_neighbor_distance_um", "equivalent_voronoi_diameter_um"],
            sample_size=1_000_000,
        )
        for column in ("mean_neighbor_distance_um",
                       "equivalent_voronoi_diameter_um"):
            values = reference_distances[column].to_numpy(float)
            values = values[np.isfinite(values)]
            if len(values):
                reference_distance_max = max(
                    reference_distance_max, float(values.max())
                )
        reference_umi = _read_metrics(
            sample, results_root, ["total_counts", "in_tissue"],
            sample_size=1_000_000,
        )
        values = reference_umi.loc[
            reference_umi.in_tissue, "total_counts"
        ].to_numpy(float)
        values = values[np.isfinite(values)]
        if len(values):
            reference_umi_max = max(reference_umi_max, float(values.max()))
    reference_distance_max = max(reference_distance_max, .1)
    reference_umi_max = max(reference_umi_max, 1.0)

    for tissue in sorted({s["sample_type"] for s in samples}):
        group=[s for s in samples if s["sample_type"]==tissue]
        tissue_output = output_dir / tissue
        tissue_output.mkdir(parents=True, exist_ok=True)
        zero_x = _zero_marker_x_positions(group)

        fig,ax=plt.subplots(figsize=(6,4.5))
        for s in group:
            c=TECHNOLOGY_COLORS[s["technology"]]; d=pd.read_csv(_result_dir(s,results_root)/"spatial_coverage_curve.csv")
            # Older saturated output files can stop at the first value rounded
            # to 100%. Explicitly retain the flat 100% tail through 50 µm.
            if d.iloc[-1, 1] >= 100 and d.iloc[-1, 0] < 50:
                differences = np.diff(d.iloc[:, 0].to_numpy(float))
                positive_steps = differences[differences > 0]
                step = float(np.median(positive_steps)) if len(positive_steps) else 1.0
                extra_x = np.arange(d.iloc[-1, 0] + step, 50 + step / 2, step)
                if len(extra_x):
                    d = pd.concat([
                        d,
                        pd.DataFrame({d.columns[0]: extra_x,
                                      d.columns[1]: np.full(len(extra_x), 100.0)}),
                    ], ignore_index=True)
            # Panel C uses a common 2 µm sampling interval even when a saved
            # coverage curve contains additional intermediate distances.
            distances = d.iloc[:, 0].to_numpy(float)
            d = d.loc[np.isclose(np.mod(distances, 2.0), 0.0, atol=1e-9)]
            ax.plot(d.iloc[:,0],d.iloc[:,1],color=c,label=_label(s)); y=d.iloc[0,1]
            marker_x = zero_x[id(s)]
            ax.scatter(marker_x,y,facecolor="white",edgecolor=c,zorder=20,
                       clip_on=False)
        ax.set(xlabel="Distance from spot edge (µm)",
               ylabel="Capture area covered (%)",xlim=(0,35))
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(NullLocator())
        _format_percent_axis(ax, upper=103)
        ax.margins(x=0)
        ax.set_axisbelow(True)
        ax.grid(which="major", axis="both", color="lightgrey", linewidth=.7)
        fig.tight_layout()
        fig.savefig(tissue_output/"C_coverage.png",dpi=dpi,bbox_inches="tight")
        plt.close(fig)

        distance_data = {}
        distance_x_values = []
        for s in group:
            frame=_read_metrics(s,results_root,["mean_neighbor_distance_um","equivalent_voronoi_diameter_um"],sample_size=1_000_000)
            for col,ls in (("mean_neighbor_distance_um","-"),("equivalent_voronoi_diameter_um","--")):
                x,y=_ecdf(frame[col])
                distance_data[(id(s), col)] = (x, y, ls)
                if len(x): distance_x_values.append(x)
        common_distance_x = (
            np.unique(np.concatenate(distance_x_values))
            if distance_x_values else np.empty(0)
        )
        common_distance_x = np.unique(np.append(
            common_distance_x, reference_distance_max
        ))

        fig,ax=plt.subplots(figsize=(6,4.5))
        for s in group:
            c=TECHNOLOGY_COLORS[s["technology"]]
            for col in ("mean_neighbor_distance_um", "equivalent_voronoi_diameter_um"):
                x, y, ls = distance_data[(id(s), col)]
                if len(x):
                    extra_x = common_distance_x[common_distance_x > x[-1]]
                    if len(extra_x):
                        x = np.concatenate((x, extra_x))
                        y = np.concatenate((y, np.ones(len(extra_x))))
                ax.plot(x,y*100,color=c,ls=ls,
                        label=_label(s) if ls=="-" else None)
        ax.set(xlabel="Distance (µm)",ylabel="Cumulative percentage (%)")
        ax.set_xscale("log")
        ax.set_xlim(.1, reference_distance_max)
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_minor_locator(
            LogLocator(base=10, subs=range(2, 10))
        )
        ax.tick_params(axis="x", which="minor", length=0, labelbottom=False)
        ax.margins(x=0)
        _format_percent_axis(ax, upper=103)
        ax.set_axisbelow(True)
        ax.grid(which="major", axis="y", color="lightgrey", linewidth=1.1)
        _add_log_subdivision_grid(ax)
        fig.tight_layout()
        fig.savefig(tissue_output/"D_distance_cdfs.png",dpi=dpi,bbox_inches="tight")
        plt.close(fig)

        umi_data = {}
        umi_x_values = []
        for s in group:
            frame=_read_metrics(s,results_root,["total_counts","in_tissue"],sample_size=1_000_000)
            values=frame.loc[frame.in_tissue,"total_counts"]
            x,y=_ecdf(values)
            umi_data[id(s)] = (x, y, values)
            if len(x): umi_x_values.append(x)
        common_umi_x = (
            np.unique(np.concatenate(umi_x_values))
            if umi_x_values else np.empty(0)
        )
        common_umi_x = np.unique(np.append(common_umi_x, reference_umi_max))

        fig,ax=plt.subplots(figsize=(6,4.5))
        for s in group:
            c=TECHNOLOGY_COLORS[s["technology"]]
            x, y, values = umi_data[id(s)]
            if len(x):
                extra_x = common_umi_x[common_umi_x > x[-1]]
                if len(extra_x):
                    x = np.concatenate((x, extra_x))
                    y = np.concatenate((y, np.ones(len(extra_x))))
            ax.step(x,y*100,where="post",color=c,label=_label(s))
            zero=float(np.mean(np.asarray(values)==0)) * 100
            ax.scatter(0,zero,facecolor="white",edgecolor=c,zorder=20,
                       clip_on=False)
        ax.set(xlabel="Recovered UMI counts per spot",
               ylabel="Cumulative percentage (%)",xlim=(0,None))
        # A symmetric-log scale preserves the meaningful x=0 markers while
        # displaying all positive UMI counts logarithmically.
        ax.set_xscale("symlog", linthresh=1, linscale=1)
        ax.set_xlim(0, reference_umi_max)
        ax.xaxis.set_minor_locator(NullLocator())
        ax.margins(x=0)
        _format_percent_axis(ax, upper=103)
        ax.set_axisbelow(True)
        ax.grid(which="major", axis="y", color="lightgrey", linewidth=1.1)
        _add_log_subdivision_grid(ax, minimum=1.0)
        fig.tight_layout()
        fig.savefig(tissue_output/"E_tissue_umi_cdf.png",dpi=dpi,bbox_inches="tight")
        plt.close(fig)


def export_all_panels(samples, results_root, output_dir, **kwargs):
    output_dir=Path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    for sample in samples:
        name=(_label(sample)+"_"+sample["dataset_id"]).replace(" ","_").replace("/","-")
        tissue_output = output_dir / sample["sample_type"]
        tissue_output.mkdir(parents=True, exist_ok=True)
        sample_geometry = _sample_geometry(sample)
        export_coverage_map(sample,results_root,tissue_output/f"{name}_A_coverage_map.png",pixel_size_um=kwargs.get("coverage_pixel_size_um",2),maximum_distance_um=kwargs.get("coverage_vmax_um",50),dpi=kwargs.get("dpi",300),microns_per_inch=kwargs.get("panel_a_microns_per_inch",2000),cividis_start=kwargs.get("panel_a_cividis_start",.15),black_and_white=kwargs.get("panel_a_black_and_white",False),sample_geometry=sample_geometry)
        export_umi_spatial_crop(sample,results_root,tissue_output/f"{name}_B_umi_crop.png",crop_size_um=kwargs.get("crop_size_um",120),dpi=kwargs.get("dpi",300),detail_inset_zoom=kwargs.get("panel_b_detail_inset_zoom",10),binary_umi=kwargs.get("panel_b_binary_umi",False),sample_geometry=sample_geometry)
    export_tissue_combined_plots(samples,results_root,output_dir,dpi=kwargs.get("dpi",300))
    export_sample_summary_statistics(
        samples, results_root,
        output_dir / "spatial_sampling_summary_statistics.csv",
        sample_size=kwargs.get("statistics_sample_size", 1_000_000),
    )
