from __future__ import annotations

import os
import csv
from dataclasses import dataclass
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import spatialdata_plot  # noqa: F401
import xarray as xr
from matplotlib import font_manager as fm
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.patches import PathPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

import sys
sys.path.append(os.path.abspath('..'))
from utils.plotting import best_datatree_scale

if TYPE_CHECKING:
    from diffusion_utils import DiffusionBinData, GeneBinData

MICROMETER_UNIT = "\u00B5m"


@dataclass
class DiffusionPlotSettings:
    tissue_outline: bool = True
    boundary: bool = True
    bin_gene_counts: bool = True
    bins_and_closest_points: bool = True
    distances_linearized: bool = True
    histogram: bool = True
    boxplot: bool = True


def get_image_figure_parameters(sdata, HE_image, dpi=200, max_dim_for_plots=10000):
    assert isinstance(
        sdata[HE_image], xr.DataTree
    ), f"Expected xarray.DataTree for sdata[{HE_image!r}], got {type(sdata[HE_image])}"

    best_scale = best_datatree_scale(sdata[HE_image], max_dim_for_plots)
    print(f"best scale: {best_scale}")

    figsize_factor = dpi
    figsize = (
        sdata.images[HE_image][best_scale]["image"].shape[2] / figsize_factor,
        sdata.images[HE_image][best_scale]["image"].shape[1] / figsize_factor,
    )
    print(f"figsize: {figsize}")

    return best_scale, dpi, figsize_factor, figsize


def get_labels_figure_parameters(sdata, labels_layer, dpi=200, max_dim_for_plots=10000):
    assert isinstance(
        sdata.labels[labels_layer], xr.DataArray
    ), f"Expected xarray.DataArray for sdata.labels[{labels_layer!r}], got {type(sdata.labels[labels_layer])}"

    figsize_factor = dpi
    figsize = (
        sdata.labels[labels_layer].shape[1] / figsize_factor,
        sdata.labels[labels_layer].shape[0] / figsize_factor,
    )
    print(f"figsize: {figsize}")

    return dpi, figsize_factor, figsize


def save_or_show_figure(fig: plt.Figure, output_path: str) -> None:
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(
        output_path,
        bbox_inches="tight",
        pad_inches=max(fig.get_size_inches()) / 100,
        facecolor="white",
        transparent=False,
    )
    print(f"{output_path} saved")
    plt.close(fig)


def pick_scalebar_length_um(
    visible_width_px: float,
    pixel_size_um: float,
    target_fraction: float = 0.2,
    nice_lengths_um: np.ndarray | None = None,
) -> float:
    assert visible_width_px > 0, "visible_width_px must be positive."
    assert pixel_size_um > 0, "pixel_size_um must be positive."
    assert 0 < target_fraction < 1, "target_fraction must be between 0 and 1."

    if nice_lengths_um is None:
        nice_lengths_um = np.asarray([20, 50, 100, 200, 500, 1000, 2000], dtype=float)
    assert np.all(nice_lengths_um > 0), "All scale bar candidates must be positive."

    visible_width_um = visible_width_px * pixel_size_um
    target_um = visible_width_um * target_fraction
    candidates = nice_lengths_um[nice_lengths_um <= target_um]
    if candidates.size == 0:
        return float(nice_lengths_um[0])
    return float(candidates[-1])


def add_scalebar(
    ax: plt.Axes,
    length_px: float,
    label: str,
    loc: str = "lower right",
    color: str = "black",
    font_size: int = 8,
) -> None:
    assert length_px > 0, "length_px must be positive."
    scalebar = AnchoredSizeBar(
        transform=ax.transData,
        size=length_px,
        label=label,
        loc=loc,
        pad=0.5,
        borderpad=1.0,
        sep=4,
        color=color,
        frameon=True,
        size_vertical=max(1.0, length_px * 0.01),
        fontproperties=fm.FontProperties(size=font_size),
    )
    scalebar.patch.set_facecolor("white")
    scalebar.patch.set_edgecolor("black")
    scalebar.patch.set_alpha(0.8)
    ax.add_artist(scalebar)


def draw_labeled_rectangles(
    ax: plt.Axes,
    rectangles: list[tuple[float, float, float, float]] | None,
) -> None:
    if not rectangles:
        return

    for i, (xmin, ymin, xmax, ymax) in enumerate(rectangles, start=1):
        rect = Rectangle(
            (xmin, ymin),
            xmax - xmin,
            ymax - ymin,
            fill=False,
            edgecolor="yellow",
            linewidth=2,
        )
        ax.add_patch(rect)
        ax.text(
            xmin,
            ymin,
            f"crop{i}",
            color="yellow",
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
            bbox={"facecolor": "black", "alpha": 0.6, "pad": 1.5, "edgecolor": "none"},
        )


def plot_bins_colored_by_feature(
    sdata,
    HE_image: str | None,
    technology: str,
    dataset_id: str,
    level_info,
    feature: str,
    output_path: str,
    add_scale_bar: bool | None = None,
    pixel_size_um: float | None = None,
    rectangles: list[tuple[float, float, float, float]] | None = None,
):
    # color can be for example 'total_counts' or the name of a gene.

    title = f"{feature} ({technology} {dataset_id} {level_info.name})"
    should_add_scale_bar = (pixel_size_um is not None) if add_scale_bar is None else add_scale_bar

    if HE_image:
        best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        sdata.pl.render_images(element=HE_image, scale=best_scale).pl.render_labels(
            element=level_info.labels_layer,
            fill_alpha=0.5,
            contour_px=5,
            outline_alpha=0,
            color=feature,
            table_name=level_info.table_name,
        ).pl.show(
            coordinate_systems=dataset_id,
            title=title,
            figsize=figsize,
            colorbar=False,
            dpi=dpi,
            ax=ax,
        )
        if should_add_scale_bar:
            assert pixel_size_um is not None and pixel_size_um > 0, (
                "pixel_size_um must be set to a positive number when add_scale_bar=True."
            )
            x0, x1 = ax.get_xlim()
            visible_width_px = abs(x1 - x0)
            bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
            bar_px = bar_um / pixel_size_um
            add_scalebar(ax=ax, length_px=bar_px, label=f"{int(bar_um)} {MICROMETER_UNIT}")
        draw_labeled_rectangles(ax=ax, rectangles=rectangles)
        save_or_show_figure(fig, output_path)

    else:
        dpi, _, figsize = get_labels_figure_parameters(sdata, level_info.labels_layer)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        sdata.pl.render_labels(
            element=level_info.labels_layer,
            fill_alpha=0.5,
            contour_px=5,
            outline_alpha=0,
            color=feature,
            table_name=level_info.table_name,
        ).pl.show(
            coordinate_systems=dataset_id,
            title=title,
            figsize=figsize,
            colorbar=False,
            dpi=dpi,
            ax=ax,
        )
        if should_add_scale_bar:
            assert pixel_size_um is not None and pixel_size_um > 0, (
                "pixel_size_um must be set to a positive number when add_scale_bar=True."
            )
            x0, x1 = ax.get_xlim()
            visible_width_px = abs(x1 - x0)
            bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
            bar_px = bar_um / pixel_size_um
            add_scalebar(ax=ax, length_px=bar_px, label=f"{int(bar_um)} {MICROMETER_UNIT}")
        draw_labeled_rectangles(ax=ax, rectangles=rectangles)
        save_or_show_figure(fig, output_path)


def plot_tissue_outline_on_HE_image(
    sdata,
    HE_image: str,
    tissue_mask: str,
    technology: str,
    dataset_id: str,
    output_path: str,
    fill_alpha=0.3,
    outline_width=5,
    bounds: tuple[float, float, float, float] | None = None,
    add_scale_bar: bool = False,
    pixel_size_um: float | None = None,
    rectangles: list[tuple[float, float, float, float]] | None = None,
):
    # Plot H&E image with tissue outline on top.
    # Optionally restrict the plot to bounds = (xmin, ymin, xmax, ymax).

    best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

    title = f"{technology} {dataset_id}: H&E with tissue outline"
    if bounds is not None:
        xmin, ymin, xmax, ymax = bounds
        title = f"{title}\n({xmin} ≤ x ≤ {xmax}, {ymin} ≤ y ≤ {ymax})"

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    sdata.pl.render_images(element=HE_image, scale=best_scale).pl.render_shapes(
        element=tissue_mask,
        color="#1f77b4",
        fill_alpha=fill_alpha,
        outline_width=outline_width,
        outline_alpha=1,
        outline_color="red",
    ).pl.show(
        coordinate_systems=dataset_id,
        title=title,
        figsize=figsize,
        colorbar=False,
        dpi=dpi,
        ax=ax,
    )
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    if bounds is not None:
        # Restrict plot to the requested bounds and show the same bounds on axes.
        xmin, ymin, xmax, ymax = bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
    if add_scale_bar:
        assert pixel_size_um is not None and pixel_size_um > 0, (
            "pixel_size_um must be set to a positive number when add_scale_bar=True."
        )
        x0, x1 = ax.get_xlim()
        visible_width_px = abs(x1 - x0)
        bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
        bar_px = bar_um / pixel_size_um
        add_scalebar(ax=ax, length_px=bar_px, label=f"{int(bar_um)} {MICROMETER_UNIT}")
    draw_labeled_rectangles(ax=ax, rectangles=rectangles)

    save_or_show_figure(fig, output_path)


def plot_bin_outlines(crop_rect, bins_center_coordinates, bin_outlines_gdf, title: str):
    # crop_rect = [xmin, ymin, xmax, ymax]
    # bins_center_coordinates = numpy array of nx2
    # bin_outlines_gdf: geopandas dataframe with outlines of the bins
    rect = Rectangle(
        (crop_rect[0], crop_rect[1]),  # lower-left corner
        crop_rect[2] - crop_rect[0],  # width
        crop_rect[3] - crop_rect[1],  # height
        fill=False,
        edgecolor="green",
        linewidth=1,
    )
    _, ax = plt.subplots(figsize=(8, 8))
    bin_outlines_gdf.plot(ax=ax, facecolor="none")
    ax.set_title(title)
    ax.scatter(x=bins_center_coordinates[:, 0], y=bins_center_coordinates[:, 1], color="red", s=1)
    ax.plot
    ax.add_patch(rect)


def plot_boundary_geoseries(
    boundary: BaseGeometry,
    output_path: str,
    figsize: tuple[float, float] = (6, 6),
):
    fig, ax = plt.subplots(figsize=figsize)

    gpd.GeoSeries([boundary]).plot(ax=ax, figsize=figsize)
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.invert_yaxis()

    save_or_show_figure(fig, output_path)


def _assert_polygon_ring_closed(ring, ring_name: str) -> None:
    coords = list(ring.coords)
    assert len(coords) >= 4, f"{ring_name} must have at least 4 points (including closure point)."
    assert coords[0] == coords[-1], f"{ring_name} is not closed."


def _polygon_to_path(polygon: Polygon) -> Path:
    vertices: list[tuple[float, float]] = []
    codes: list[int] = []

    def add_ring(ring) -> None:
        ring_coords = np.asarray(ring.coords, dtype=float)
        n = ring_coords.shape[0]
        vertices.extend(ring_coords.tolist())
        codes.extend([Path.MOVETO] + [Path.LINETO] * (n - 2) + [Path.CLOSEPOLY])

    add_ring(polygon.exterior)
    for interior in polygon.interiors:
        add_ring(interior)

    return Path(np.asarray(vertices, dtype=float), codes)


def _ring_signed_area(coords: np.ndarray | list[tuple[float, float]]) -> float:
    ring = np.asarray(coords, dtype=float)
    if ring.ndim != 2 or ring.shape[0] < 3:
        return 0.0
    if np.allclose(ring[0], ring[-1]):
        ring = ring[:-1]
    if ring.shape[0] < 3:
        return 0.0
    x = ring[:, 0]
    y = ring[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _ring_as_array(
    polygon: Polygon,
    *,
    tolerance: float = 100.0,
) -> np.ndarray:
    """
    Return the exterior ring used for orientation markers.
    The tolerance is in pixels, larger values yield coarser polygon approximations.
    """
    if tolerance > 0:
        simplified = polygon.simplify(tolerance, preserve_topology=True)
        if isinstance(simplified, Polygon):
            ring = np.asarray(simplified.exterior.coords, dtype=float)
            if ring.shape[0] >= 4:
                return ring
    return np.asarray(polygon.exterior.coords, dtype=float)


def _add_polygon_orientation_arrows(
    ax: plt.Axes,
    polygon: Polygon,
    color: tuple[float, ...],
    base_spacing_px: float = 250.0,
) -> None:
    """
    Draw orientation arrows along the polygon boundary.

    The direction is derived from a simplified polygon to avoid local jitter.

    Code is for debugging / testing purposes only.
    """
    coords = np.asarray(polygon.exterior.coords, dtype=float)
    if coords.shape[0] < 4:
        return

    ring = coords
    if np.allclose(ring[0], ring[-1]):
        ring = ring[:-1]
    if ring.shape[0] < 3:
        return

    simplified_ring = _ring_as_array(polygon)
    if np.allclose(simplified_ring[0], simplified_ring[-1]):
        simplified_ring = simplified_ring[:-1]
    if simplified_ring.shape[0] < 3:
        simplified_ring = ring
    is_simplified_ccw = _ring_signed_area(simplified_ring) > 0.0

    seg_starts = ring
    seg_ends = np.vstack([ring[1:], ring[0]])
    seg_vec = seg_ends - seg_starts
    seg_len = np.hypot(seg_vec[:, 0], seg_vec[:, 1])
    seg_len = np.where(seg_len > 0, seg_len, 1e-12)
    cumulative = np.concatenate(([0.0], np.cumsum(seg_len)))
    perimeter = float(cumulative[-1])
    if perimeter <= 0:
        return

    simp_starts = simplified_ring
    simp_ends = np.vstack([simplified_ring[1:], simplified_ring[0]])
    simp_vec = simp_ends - simp_starts
    simp_len = np.hypot(simp_vec[:, 0], simp_vec[:, 1])
    simp_len = np.where(simp_len > 0, simp_len, 1e-12)
    simp_cumulative = np.concatenate(([0.0], np.cumsum(simp_len)))
    simp_perimeter = float(simp_cumulative[-1])
    if simp_perimeter <= 0:
        return

    x0, y0, x1, y1 = polygon.bounds
    short_side = max(1e-12, min(x1 - x0, y1 - y0))
    arrow_length = max(6.0, 0.5 * short_side / 20.0)
    head_length = arrow_length / 3.0

    n_arrows = max(2, int(np.ceil(perimeter / base_spacing_px)))
    if n_arrows > 12:
        n_arrows = 12

    fractions = np.linspace(0.0, 1.0, n_arrows + 2, endpoint=True)[1:-1]
    for frac in fractions:
        target = frac * perimeter
        seg_idx = int(np.searchsorted(cumulative, target, side="right") - 1)
        seg_idx = np.clip(seg_idx, 0, len(seg_vec) - 1)
        if seg_len[seg_idx] <= 0:
            continue

        t = (target - cumulative[seg_idx]) / seg_len[seg_idx]
        p0 = seg_starts[seg_idx] * (1.0 - t) + seg_ends[seg_idx] * t
        simp_target = frac * simp_perimeter
        simp_seg_idx = int(np.searchsorted(simp_cumulative, simp_target, side="right") - 1)
        simp_seg_idx = np.clip(simp_seg_idx, 0, len(simp_vec) - 1)
        if simp_len[simp_seg_idx] <= 0:
            continue
        direction = simp_vec[simp_seg_idx] / simp_len[simp_seg_idx]
        if not is_simplified_ccw:
            direction = -direction

        start = p0 - 0.5 * arrow_length * direction
        ax.arrow(
            start[0],
            start[1],
            arrow_length * direction[0],
            arrow_length * direction[1],
            length_includes_head=True,
            width=0.0,
            head_width=0.25 * arrow_length,
            head_length=head_length,
            fc=color,
            ec=color,
            linewidth=1.2,
            zorder=3,
            clip_on=False,
        )


def plot_shapely_geometries(
    sdata,
    tissue_mask: str | None,
    output_path: str,
    figsize: tuple[float, float] = (10, 10),
    dpi: int | None = 600,
    title: str | None = None,
):
    # Code for debugging / testing purpose only.

    if tissue_mask is None or tissue_mask not in sdata.shapes:
        return
    gdf = sdata.shapes[tissue_mask]
    geometries = gdf.geometry

    assert len(geometries) > 0, "geometries must not be empty."
    for i, geom in enumerate(geometries):
        assert isinstance(geom, (Polygon, MultiPolygon)), (
            f"Geometry {i} must be a Polygon or MultiPolygon, got {geom.geom_type}."
        )

        polygons = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
        for j, poly in enumerate(polygons):
            _assert_polygon_ring_closed(poly.exterior, f"Geometry {i}, polygon {j} exterior")
            for k, interior in enumerate(poly.interiors):
                _assert_polygon_ring_closed(interior, f"Geometry {i}, polygon {j} hole {k}")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    cmap = plt.get_cmap("tab20", len(geometries))

    for i, geom in enumerate(geometries):
        color = cmap(i)
        polygons = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
        for poly in polygons:
            patch = PathPatch(
                _polygon_to_path(poly),
                facecolor=color,
                edgecolor=color,
                linewidth=1.5,
                alpha=0.6,
            )
            ax.add_patch(patch)
            _add_polygon_orientation_arrows(ax, poly, color=color)

    ax.autoscale()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if title is not None:
        ax.set_title(title)

    save_or_show_figure(fig, output_path)


def plot_bin_gene_counts_and_contour(
    gene_bin_data: "GeneBinData",
    contour,
    output_path: str,
    title: str,
    figsize=(10, 10),
    line_width=1,
    marker_size=20,
    dpi: int | None = None,
    add_scale_bar: bool = False,
    pixel_size_um: float | None = None,
) -> None:
    # Plot the bin centers, colored by the number of transcripts detected in that bin,
    # together with a given contour over it.
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_title(title)
    sc = ax.scatter(
        gene_bin_data.bin_centers[:, 0],
        gene_bin_data.bin_centers[:, 1],
        c=gene_bin_data.transcript_counts,
        s=marker_size,
    )
    cbar = fig.colorbar(sc, ax=ax, shrink=0.4)
    cbar.set_label("Counts")
    gpd.GeoSeries([contour]).plot(
        ax=ax,
        color="red",
        linewidth=line_width,
    )
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.invert_yaxis()
    if add_scale_bar:
        assert pixel_size_um is not None and pixel_size_um > 0, (
            "pixel_size_um must be set to a positive number when add_scale_bar=True."
        )
        x0, x1 = ax.get_xlim()
        visible_width_px = abs(x1 - x0)
        bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
        bar_px = bar_um / pixel_size_um
        add_scalebar(ax=ax, length_px=bar_px, label=f"{int(bar_um)} {MICROMETER_UNIT}")

    save_or_show_figure(fig, output_path)


def plot_bins_and_closest_points_on_contour(
    data: "DiffusionBinData",
    contour,
    output_path: str,
    title: str,
    line_width=1,
    figsize: tuple[float, float] = (6, 6),
    dpi: int | None = None,
):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot the bin centers
    ax.scatter(
        data.gene_bins.bin_centers[:, 0],
        data.gene_bins.bin_centers[:, 1],
        c="blue",
        s=1,
        label="Bin centers",
    )

    # Line segments connecting each bin to its closest point on the contour
    segments = np.stack([data.gene_bins.bin_centers, data.distance_data.closest_xy], axis=1)
    lc = LineCollection(
        segments,
        array=data.distance_data.closest_distances,  # color line by distance
        linewidths=line_width,
    )
    ax.add_collection(lc)

    # Plot the contour
    line_xy = np.array(contour.coords)
    ax.plot(line_xy[:, 0], line_xy[:, 1], "-", c="red", lw=line_width, label="Contour")

    ax.set_aspect("equal")
    # Place legend outside to the right to avoid overlapping the plot.
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    fig.subplots_adjust(right=0.78)
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.set_title(title)
    ax.invert_yaxis()

    save_or_show_figure(fig, output_path)


def plot_distances_linearized(
    data: "DiffusionBinData",
    output_path: str,
    title: str,
    pixel_size_um: float | None = None,
    figsize: tuple[float, float] = (10, 2.5),
    dpi: int | None = None,
) -> None:
    # Scatter plot of bins: along-boundary position vs distance to boundary.
    x = data.distance_data.closest_along_boundary
    y = data.distance_data.closest_distances
    counts = np.asarray(data.gene_bins.transcript_counts, dtype=float)
    positive_counts = counts[counts > 0]
    norm = None
    if positive_counts.size > 0:
        norm = LogNorm(vmin=float(np.min(positive_counts)), vmax=float(np.max(positive_counts)))

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    sc = ax.scatter(
        x,
        y,
        s=3.0,
        alpha=0.9,
        c=counts,
        cmap="viridis",
        norm=norm,
        marker=".",
        edgecolors="none",
        linewidths=0,
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=1.2)  # WAS: 0.45
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label("Transcript count per bin")
    ax.axhline(0, color="red", linewidth=1, linestyle=":", zorder=3)
    ax.set_xlabel("Normalized position along boundary")
    ax.set_ylabel("Distance to boundary [px]")

    if pixel_size_um is not None and pixel_size_um > 0:
        secax = ax.secondary_yaxis(
            "right",
            functions=(
                lambda y_px: y_px * pixel_size_um,
                lambda y_um: y_um / pixel_size_um,
            ),
        )
        secax.yaxis.set_label_position("right")
        secax.yaxis.tick_right()
        secax.set_ylabel(f"Distance to boundary [{MICROMETER_UNIT}]", labelpad=14)

    ax.set_title(title)

    save_or_show_figure(fig, output_path)


def plot_distance_histogram(
    data: "DiffusionBinData",
    output_path: str,
    unit: str = "px",
    pixel_size_um: float | None = None,
    title: str | None = None,
) -> None:
    distances = np.asarray(data.distance_data.closest_distances, dtype=float)
    assert unit in {"px", "um"}, "unit must be 'px' or 'um'."
    assert unit == "px" or (pixel_size_um is not None and pixel_size_um > 0), (
        "pixel_size_um must be set to a positive number when unit='um'."
    )
    if unit == "um":
        distances = distances * pixel_size_um
        unit_label = MICROMETER_UNIT
        bin_width = 5.0
    else:
        unit_label = "px"
        # Use a ~5 µm bin width in pixel space when conversion is available.
        bin_width = 5.0 / pixel_size_um if pixel_size_um is not None and pixel_size_um > 0 else 5.0

    fig, ax = plt.subplots(figsize=(6, 6))
    dmin = float(np.min(distances))
    dmax = float(np.max(distances))
    bins = np.arange(dmin, dmax + bin_width, bin_width)
    if bins.size < 2:
        bins = 10
    ax.hist(distances, bins=bins)
    ax.axvline(0, color="black", linewidth=1, linestyle=":", zorder=3)
    ax.set_xlabel(f"Distance to boundary [{unit_label}]")
    ax.set_ylabel("Count")
    ax.set_title(title or f"Distance to boundary histogram [{unit_label}]")
    save_or_show_figure(fig, output_path)


def combined_distance_histogram(
    crop_specs: list[tuple[str, str, float, str]],
    output_path: str,
    plots_root: str,
    title: str | None = None,
    figsize = (12, 6),
    align_peaks: bool = False,
) -> None:
    """
    Plot overlaid distance histograms (as lines) for multiple crops, in micrometers.

    If `align_peaks` is True, each curve is shifted so its histogram peak is placed
    at x=0.

    Each entry in crop_specs is (technology, dataset_id, pixel_size_um, crop_name).
    """
    data_per_label = _read_crop_distance_data(crop_specs, plots_root)
    bins = _distance_hist_bins(data_per_label, bin_width=5.0)
    curves = _distance_hist_curves(data_per_label, bins=bins, align_peaks=align_peaks, density=False)

    fig, ax = plt.subplots(figsize=figsize)
    for label, centers, counts in curves:
        ax.plot(centers, counts, linewidth=1.5, label=label)

    ax.axvline(0, color="black", linewidth=1, linestyle=":", zorder=3)
    x_label = (
        f"Peak-aligned distance to boundary [{MICROMETER_UNIT}]"
        if align_peaks
        else f"Distance to boundary [{MICROMETER_UNIT}]"
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel("Count")
    ax.set_title(title or f"Combined distance histogram")
    ax.legend(fontsize="small")

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(
        output_path,
        bbox_inches="tight",
        pad_inches=max(fig.get_size_inches()) / 100,
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)


def _read_crop_distance_data(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
) -> list[tuple[str, np.ndarray]]:
    assert len(crop_specs) > 0, "crop_specs must contain at least one entry."

    data_per_label: list[tuple[str, np.ndarray]] = []

    for technology, dataset_id, pixel_size_um, crop_name in crop_specs:
        assert pixel_size_um > 0, "pixel_size_um must be > 0."
        csv_path = os.path.join(
            plots_root,
            technology,
            dataset_id,
            "analysis",
            crop_name,
            "distances.csv",
        )
        assert os.path.exists(csv_path), f"Missing file: {csv_path}"

        with open(csv_path, "r", encoding="utf-8", newline="") as csv_file:
            reader = csv.DictReader(csv_file)
            assert reader.fieldnames is not None and "closest_distance_px" in reader.fieldnames, (
                f"{csv_path} must contain a 'closest_distance_px' column"
            )
            distances_px = [float(row["closest_distance_px"]) for row in reader if row.get("closest_distance_px")]

        if len(distances_px) == 0:
            continue

        distances_um = np.asarray(distances_px, dtype=float) * pixel_size_um
        label = f"{technology}, {dataset_id}, {crop_name}"
        data_per_label.append((label, distances_um))

    if len(data_per_label) == 0:
        raise ValueError("No distance rows found in any selected crop.")

    return data_per_label


def _distance_hist_bins(data_per_label: list[tuple[str, np.ndarray]], bin_width: float = 5.0) -> np.ndarray:
    concat = np.concatenate([distances_um for _, distances_um in data_per_label])
    dmin = float(np.min(concat))
    dmax = float(np.max(concat))
    bins = np.arange(dmin, dmax + bin_width, bin_width)
    if bins.size < 2:
        bins = 10
    return bins


def _distance_hist_curves(
    data_per_label: list[tuple[str, np.ndarray]],
    bins: np.ndarray,
    align_peaks: bool = False,
    density: bool = False,
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    curves: list[tuple[str, np.ndarray, np.ndarray]] = []
    for label, distances_um in data_per_label:
        counts, edges = np.histogram(distances_um, bins=bins, density=density)
        centers = 0.5 * (edges[:-1] + edges[1:])
        if align_peaks:
            peak_idx = int(np.argmax(counts))
            centers = centers - centers[peak_idx]
        curves.append((label, centers, counts))
    return curves

def plot_distance_boxplot(
    data: "DiffusionBinData",
    output_path: str,
    unit: str = "px",
    pixel_size_um: float | None = None,
    title: str | None = None,
) -> None:
    distances = np.asarray(data.distance_data.closest_distances, dtype=float)
    assert unit in {"px", "um"}, "unit must be 'px' or 'um'."
    assert unit == "px" or (pixel_size_um is not None and pixel_size_um > 0), (
        "pixel_size_um must be set to a positive number when unit='um'."
    )
    if unit == "um":
        distances = distances * pixel_size_um
        unit_label = MICROMETER_UNIT
    else:
        unit_label = "px"

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.boxplot(distances)
    ax.set_ylabel(f"Distance to boundary [{unit_label}]")
    ax.set_title(title or f"Distance to boundary boxplot [{unit_label}]")
    save_or_show_figure(fig, output_path)


def ridgeline_distance_histogram(
    crop_specs: list[tuple[str, str, float, str]],
    output_path: str,
    plots_root: str,
    title: str | None = None,
    figsize=(14, 8),
    align_peaks: bool = False,
    overlap: float = 0.6,
    fill_alpha: float = 0.45,
) -> None:
    """
    Plot distance histograms as a ridgeline (joy) plot for multiple crops.

    If `align_peaks` is True, each curve is shifted so its histogram peak is placed
    at x=0.
    """
    assert len(crop_specs) > 0, "crop_specs must contain at least one entry."
    data_per_label = _read_crop_distance_data(crop_specs, plots_root)
    bins = _distance_hist_bins(data_per_label, bin_width=5.0)
    all_curves = _distance_hist_curves(data_per_label, bins=bins, align_peaks=align_peaks, density=True)

    max_height = max(float(np.max(counts)) for _, _, counts in all_curves)
    if max_height <= 0:
        raise ValueError("No positive histogram values found.")

    normalized_curves = [
        (label, centers, counts / max_height)
        for label, centers, counts in all_curves
    ]

    cmap = plt.get_cmap("viridis")
    step = 1.0 - overlap
    assert 0 <= overlap < 1, "overlap must be in [0, 1)."

    fig, ax = plt.subplots(figsize=figsize)
    num_curves = len(normalized_curves)
    for idx, (label, centers, counts) in enumerate(normalized_curves):
        baseline = (num_curves - 1 - idx) * step
        color = cmap((idx + 1) / max(len(normalized_curves), 1))
        filled = baseline + counts
        ax.plot(centers, filled, color=color, linewidth=1.2, label=label)
        ax.fill_between(centers, baseline, filled, color=color, alpha=fill_alpha, linewidth=0)

    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", length=0)
    ax.set_ylabel("")
    ax.legend(fontsize="small", bbox_to_anchor=(1.01, 1), loc="upper left")

    ax.set_xlabel(
        f"Peak-aligned distance to boundary [{MICROMETER_UNIT}]"
        if align_peaks
        else f"Distance to boundary [{MICROMETER_UNIT}]"
    )
    ax.set_title(title or "Distance ridgeline histogram")
    ax.axvline(0, color="black", linewidth=1, linestyle=":", zorder=3)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(
        output_path,
        bbox_inches="tight",
        pad_inches=max(fig.get_size_inches()) / 100,
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)
