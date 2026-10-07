from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from time import perf_counter
from typing import TYPE_CHECKING, Sequence

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import shapely
import spatialdata_plot  # noqa: F401
import xarray as xr
from matplotlib import font_manager as fm
from matplotlib.colors import LogNorm, Normalize
from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.patches import PathPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from shapely.geometry import MultiPoint, MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry

import sys
sys.path.append(os.path.abspath('..'))
from utils.plotting import best_datatree_scale
from colors import COLOR_PALETTE

if TYPE_CHECKING:
    from diffusion_utils import DiffusionBinData, GeneBinData

MICROMETER_UNIT = "\u00B5m"
KDE_GRID_STEP_UM = 0.25
KDE_TAIL_BANDWIDTHS = 3.0
KDE_BANDWIDTH_FALLBACK_UM = 1.0


def _estimate_bin_outline_radius(points_xy: np.ndarray) -> float | None:
    points_xy = np.asarray(points_xy, dtype=float).reshape(-1, 2)
    if points_xy.shape[0] < 2:
        return None

    candidate_steps: list[float] = []
    for axis in (0, 1):
        ordered = points_xy[np.argsort(points_xy[:, axis], kind="mergesort")]
        diffs = np.linalg.norm(np.diff(ordered, axis=0), axis=1)
        diffs = diffs[np.isfinite(diffs) & (diffs > 1e-8)]
        if diffs.size > 0:
            candidate_steps.append(float(np.percentile(diffs, 5.0)))

    if len(candidate_steps) == 0:
        return None

    step = min(candidate_steps)
    if not np.isfinite(step) or step <= 0.0:
        return None
    return 0.55 * step


def _bin_region_geometry(bin_centers: np.ndarray) -> BaseGeometry:
    bin_centers = np.asarray(bin_centers, dtype=float).reshape(-1, 2)
    if bin_centers.shape[0] == 0:
        raise ValueError("bin_centers must contain at least one point.")

    point_cloud = MultiPoint(bin_centers)
    if point_cloud.is_empty:
        raise ValueError("Failed to build a point cloud from bin centers.")

    if hasattr(shapely, "concave_hull") and bin_centers.shape[0] >= 4:
        region = shapely.concave_hull(point_cloud, ratio=0.02, allow_holes=False)
    else:
        region = point_cloud.convex_hull

    if region is None or region.is_empty:
        raise ValueError("Failed to build analysis-region geometry from bin centers.")

    radius = _estimate_bin_outline_radius(bin_centers)
    if radius is not None:
        region = region.buffer(radius, join_style="round")
        if region is None or region.is_empty:
            raise ValueError("Failed to pad analysis-region geometry.")
    return region


def _plot_bin_region_boundary(
    ax: plt.Axes,
    bin_centers: np.ndarray,
    edgecolor: str,
    linewidth: float,
    zorder: float,
) -> None:
    region = _bin_region_geometry(bin_centers)
    region_boundary = region.boundary
    if region_boundary is None or region_boundary.is_empty:
        raise ValueError("Failed to build analysis-region boundary from bin centers.")

    gpd.GeoSeries([region_boundary]).plot(
        ax=ax,
        color=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )


def _plot_geometry_boundary(
    ax: plt.Axes,
    geometry: BaseGeometry,
    edgecolor: str,
    linewidth: float,
    zorder: float,
) -> None:
    geometry_boundary = geometry.boundary
    if geometry_boundary is not None and not geometry_boundary.is_empty:
        gpd.GeoSeries([geometry_boundary]).plot(
            ax=ax,
            color=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )


@dataclass
class DiffusionPlotSettings:
    tissue_outline: bool = True
    boundary: bool = True
    bin_gene_counts: bool = True
    full_tissue_bin_gene_counts: bool = True
    bins_and_closest_points: bool = True
    distances_linearized: bool = True
    histogram: bool = True
    histogram_bin_widths_um: list[float] = field(default_factory=lambda: [5.0])
    max_plot_distance_um: float = 100.0
    boxplot: bool = False


@dataclass(frozen=True)
class BoundarySectionOverlay:
    """One analyzed boundary section and its distance zone, for whole-tissue plots.

    An explicit `analysis_region` preserves interior rings for closed sections;
    otherwise the zone outline is estimated from `analysis_bin_centers`.
    """
    boundary_section: BaseGeometry
    analysis_region: BaseGeometry | None = None
    analysis_bin_centers: np.ndarray | None = None
    approx_points: np.ndarray | None = None
    snapped_points: np.ndarray | None = None


def full_tissue_plot_bounds(
    full_tissue_contour: BaseGeometry,
    bin_centers: np.ndarray,
    section_overlays: Sequence[BoundarySectionOverlay] = (),
    margin_fraction: float = 0.05,
) -> tuple[float, float, float, float]:
    """(xmin, ymin, xmax, ymax) of everything drawn on a whole-tissue plot.

    Padded like matplotlib's default autoscale margins, so that every
    whole-tissue figure shows the same extent.
    """
    extents = [np.asarray(full_tissue_contour.bounds, dtype=float)]
    points = np.asarray(bin_centers, dtype=float).reshape(-1, 2)
    if points.size > 0:
        extents.append(np.concatenate([points.min(axis=0), points.max(axis=0)]))
    for overlay in section_overlays:
        for geometry in (overlay.boundary_section, overlay.analysis_region):
            if geometry is not None and not geometry.is_empty:
                extents.append(np.asarray(geometry.bounds, dtype=float))
        if overlay.analysis_bin_centers is not None and np.size(overlay.analysis_bin_centers) > 0:
            centers = np.asarray(overlay.analysis_bin_centers, dtype=float).reshape(-1, 2)
            extents.append(np.concatenate([centers.min(axis=0), centers.max(axis=0)]))
    extents = np.vstack(extents)
    xmin, ymin = extents[:, :2].min(axis=0)
    xmax, ymax = extents[:, 2:].max(axis=0)
    dx = (xmax - xmin) * margin_fraction
    dy = (ymax - ymin) * margin_fraction
    return float(xmin - dx), float(ymin - dy), float(xmax + dx), float(ymax + dy)


@dataclass(frozen=True)
class DistanceDistribution:
    label: str
    distances_um: np.ndarray
    weights: np.ndarray
    technology: str | None = None
    dataset_id: str | None = None


@dataclass(frozen=True)
class DistanceHistogramHWHM:
    label: str
    left_hwhm_um: float
    right_hwhm_um: float
    peak_position_um: float
    technology: str | None = None


@dataclass(frozen=True)
class DistanceKDEPeakMedianDistance:
    label: str
    left_median_distance_um: float
    right_median_distance_um: float
    peak_position_um: float
    technology: str | None = None


@dataclass(frozen=True)
class DistanceKDECurve:
    label: str
    x_um: np.ndarray
    density: np.ndarray
    peak_position_um: float
    peak_density: float
    alignment_shift_um: float
    technology: str | None = None


@dataclass(frozen=True)
class DistanceDistributionSummary:
    label: str
    quantile_signature_um: dict[str, float]
    cdf_distance_um: dict[str, float]
    side_mass_fraction: dict[str, float]


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


def get_labels_figure_parameters(sdata, labels_layer, dpi=200):
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


def is_single_channel_image(sdata, image_element: str, scale: str) -> bool:
    image = sdata.images[image_element][scale]["image"]
    if "c" in image.dims:
        return image.sizes["c"] == 1
    return image.ndim == 2 or (image.ndim >= 3 and image.shape[0] == 1)


def render_background_image_in_grayscale(ax: plt.Axes) -> None:
    images = ax.get_images()
    if not images:
        return

    background_image = images[0]
    image_array = background_image.get_array()
    image_ndim = getattr(image_array, "ndim", None)
    if image_ndim == 2:
        background_image.set_cmap("gray_r")
        return
    if image_ndim == 3 and image_array.shape[-1] in (3, 4):
        rgb_image = np.asarray(image_array)[..., :3]
        grayscale_image = np.tensordot(rgb_image, [0.2989, 0.5870, 0.1140], axes=([-1], [0]))
        background_image.set_data(grayscale_image)
        background_image.set_cmap("gray")


def save_or_show_figure(fig: plt.Figure, output_path: str, dpi: int = 300) -> None:
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=max(fig.get_size_inches()) / 100,
        facecolor="white",
        transparent=False,
    )
    print(f"Saved {output_path}")
    plt.close(fig)


def pick_scalebar_length_um(
    visible_width_px: float,
    pixel_size_um: float,
    target_fraction: float = 0.2,
    nice_lengths_um: np.ndarray | None = None,
) -> float:
    assert visible_width_px > 0, "visible_width_px must be positive."
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
    frameon: bool = True,
    size_vertical_scale: float = 1.0,
) -> None:
    assert length_px > 0, "length_px must be positive."
    assert size_vertical_scale > 0, "size_vertical_scale must be positive."
    size_vertical = max(1.0, length_px * 0.01) * size_vertical_scale
    scalebar = AnchoredSizeBar(
        transform=ax.transData,
        size=length_px,
        label=label,
        loc=loc,
        pad=0.5,
        borderpad=1.0,
        sep=4,
        color=color,
        frameon=frameon,
        size_vertical=size_vertical,
        fontproperties=fm.FontProperties(size=font_size),
    )
    if frameon:
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


def draw_boundary_section_overlays(
    ax: plt.Axes,
    overlay_region: BaseGeometry | None = None,
    overlay_bin_centers: np.ndarray | None = None,
    overlay_contour: BaseGeometry | None = None,
    approx_points: np.ndarray | None = None,
    snapped_points: np.ndarray | None = None,
    show_approx_points: bool = False,
) -> None:
    if overlay_region is not None and not overlay_region.is_empty:
        _plot_geometry_boundary(
            ax=ax,
            geometry=overlay_region,
            edgecolor="#f39c12",
            linewidth=1.5,
            zorder=6,
        )

    if overlay_bin_centers is not None and np.size(overlay_bin_centers) > 0:
        _plot_bin_region_boundary(
            ax=ax,
            bin_centers=overlay_bin_centers,
            edgecolor="#f39c12",
            linewidth=1.5,
            zorder=6,
        )

    if overlay_contour is not None and not overlay_contour.is_empty:
        gpd.GeoSeries([overlay_contour]).plot(
            ax=ax,
            color="red",
            linewidth=2.0,
            zorder=7,
        )

    if show_approx_points and approx_points is not None and approx_points.size > 0:
        approx_points = np.asarray(approx_points, dtype=float).reshape(-1, 2)
        ax.scatter(
            approx_points[:, 0],
            approx_points[:, 1],
            c="black",
            marker="x",
            s=45,
            linewidths=1.6,
            zorder=8,
        )

    if snapped_points is not None and snapped_points.size > 0:
        snapped_points = np.asarray(snapped_points, dtype=float).reshape(-1, 2)
        ax.scatter(
            snapped_points[:, 0],
            snapped_points[:, 1],
            c="#ff00aa",
            marker="o",
            s=28,
            edgecolors="black",
            linewidths=0.7,
            zorder=9,
        )


def plot_bins_colored_by_feature(
    sdata,
    HE_image: str | None,
    technology: str,
    dataset_id: str,
    pixel_size_um: float,
    level_info,
    feature: str,
    output_path: str,
    rectangles: list[tuple[float, float, float, float]] | None = None,
):
    def _style_colorbar_for_labels(ax: plt.Axes) -> None:
        # helper function to shrink the vertical size of the colorbar,
        # and to add a label.
        fig = ax.figure
        if len(fig.axes) < 2:
            return
        cbar_ax = fig.axes[-1]
        pos = cbar_ax.get_position()
        new_h = pos.height * 0.5
        new_y = pos.y0 + (pos.height - new_h) / 2.0
        cbar_ax.set_position([pos.x0, new_y, pos.width, new_h])
        cbar_ax.set_ylabel(f"{feature} transcript count")

    title = f"{technology}"
    if HE_image:
        best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        sdata.pl.render_images(
            element=HE_image,
            scale=best_scale,
            colorbar=False,
        ).pl.render_labels(
            element=level_info.labels_layer,
            fill_alpha=0.5,
            contour_px=5,
            outline_alpha=0,
            color=feature,
            table_name=level_info.table_name,
            colorbar_params={"fraction": 0.0375, "pad": 0.015},
        ).pl.show(
            coordinate_systems=dataset_id,
            title=title,
            figsize=figsize,
            colorbar=True,
            dpi=dpi,
            ax=ax,
        )
        render_background_image_in_grayscale(ax)
        _style_colorbar_for_labels(ax)
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
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
            colorbar=True,
            dpi=dpi,
            ax=ax,
        )
        _style_colorbar_for_labels(ax)
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
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
    pixel_size_um: float,
    output_path: str,
    fill_alpha=0.3,
    outline_width=5,
    bounds: tuple[float, float, float, float] | None = None,
    rectangles: list[tuple[float, float, float, float]] | None = None,
    section_overlays: Sequence[BoundarySectionOverlay] = (),
    show_approx_points: bool = False,
):
    # Plot H&E image with tissue outline on top.
    # Optionally restrict the plot to bounds = (xmin, ymin, xmax, ymax).

    best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

    title = technology

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    sdata.pl.render_images(element=HE_image, scale=best_scale).pl.render_shapes(
        element=tissue_mask,
        color="#1f77b4",
        fill_alpha=fill_alpha,
        outline_width=outline_width,
        outline_alpha=(1.0, 0.0),
        outline_color="red",
    ).pl.show(
        coordinate_systems=dataset_id,
        title=title,
        figsize=figsize,
        colorbar=False,
        dpi=dpi,
        ax=ax,
    )
    if is_single_channel_image(sdata, HE_image, best_scale):
        render_background_image_in_grayscale(ax)
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    if bounds is not None:
        # Restrict plot to the requested bounds and show the same bounds on axes.
        xmin, ymin, xmax, ymax = bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
    x0, x1 = ax.get_xlim()
    visible_width_px = abs(x1 - x0)
    bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
    bar_px = bar_um / pixel_size_um
    add_scalebar(ax=ax, length_px=bar_px, label=f"{int(bar_um)} {MICROMETER_UNIT}")
    draw_labeled_rectangles(ax=ax, rectangles=rectangles)
    for overlay in section_overlays:
        draw_boundary_section_overlays(
            ax=ax,
            overlay_region=overlay.analysis_region,
            overlay_bin_centers=(
                overlay.analysis_bin_centers if overlay.analysis_region is None else None
            ),
            overlay_contour=overlay.boundary_section,
            approx_points=overlay.approx_points,
            snapped_points=overlay.snapped_points,
            show_approx_points=show_approx_points,
        )

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

        # Ensure arrows follow the traversal direction of the original exterior ring
        # (from the first point through subsequent vertices), regardless of CW/CCW orientation.
        local_direction = seg_vec[seg_idx] / seg_len[seg_idx]
        if float(np.dot(direction, local_direction)) < 0.0:
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
    """Plot the full mask and each geometry with matching full-mask axis limits."""
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

    def draw_geometry(ax, geom, color):
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

    for i, geom in enumerate(geometries):
        draw_geometry(ax, geom, cmap(i))

    ax.autoscale()
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.invert_yaxis()
    if title is not None:
        ax.set_title(title)

    # Capture the combined plot's limits, including padding and inverted y-axis.
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    save_or_show_figure(fig, output_path)


def plot_bin_gene_counts_and_contour(
    gene_bin_data: "GeneBinData",
    contour,
    output_path: str,
    pixel_size_um: float,
    title: str,
    figsize=(10, 10),
    line_width=1,
    marker_size=0.5,
    show_axes: bool = False,
    dpi: int | None = None,
) -> None:
    """Scatter plot non-zero gene bin counts, overlaid with the given contour."""
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_title(title)
    counts = np.asarray(gene_bin_data.transcript_counts)
    nonzero_mask = counts > 0
    sc = ax.scatter(
        gene_bin_data.bin_centers[nonzero_mask, 0],
        gene_bin_data.bin_centers[nonzero_mask, 1],
        c=counts[nonzero_mask],
        s=marker_size,
        marker="o",
        linewidths=0,
        edgecolors="none",
    )

    if np.any(nonzero_mask):
        cbar = fig.colorbar(sc, ax=ax, shrink=0.4)
        cbar.set_label(f"{gene_bin_data.gene_name} transcript count")

    gpd.GeoSeries([contour]).plot(
        ax=ax,
        color="red",
        linewidth=line_width,
    )
    if show_axes:
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labelleft=False,
        )
    ax.invert_yaxis()
    x0, x1 = ax.get_xlim()
    visible_width_px = abs(x1 - x0)
    bar_um = pick_scalebar_length_um(visible_width_px=visible_width_px, pixel_size_um=pixel_size_um)
    bar_px = bar_um / pixel_size_um
    add_scalebar(
        ax=ax,
        length_px=bar_px,
        label=f"{int(bar_um)} {MICROMETER_UNIT}",
        loc="lower left",
        frameon=False,
        size_vertical_scale=18.0,
    )

    save_or_show_figure(fig, output_path)


def plot_full_tissue_bin_gene_counts_with_highlighted_sections(
    gene_bin_data: "GeneBinData",
    full_tissue_contour,
    section_overlays: Sequence[BoundarySectionOverlay],
    output_path: str,
    pixel_size_um: float,
    title: str,
    figsize=(20, 20),
    contour_line_width=3,
    marker_size=1,
    show_axes: bool = False,
    dpi: int | None = None,
    bounds: tuple[float, float, float, float] | None = None,
) -> None:
    """Plot gene bins, the whole tissue contour (red) and all analysis zones (cyan).

    `bounds` = (xmin, ymin, xmax, ymax) defaults to `full_tissue_plot_bounds`.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_title(title)
    shared_line_width = contour_line_width
    full_tissue_zorder = 1
    highlighted_analysis_region_zorder = 2
    transcripts_zorder = 3
    gpd.GeoSeries([full_tissue_contour]).plot(
        ax=ax,
        color="red",
        linewidth=shared_line_width,
        zorder=full_tissue_zorder,
    )
    counts = np.asarray(gene_bin_data.transcript_counts)
    nonzero_mask = counts > 0
    sc = ax.scatter(
        gene_bin_data.bin_centers[nonzero_mask, 0],
        gene_bin_data.bin_centers[nonzero_mask, 1],
        c=counts[nonzero_mask],
        s=marker_size,
        marker="o",
        linewidths=0,
        edgecolors="none",
        zorder=transcripts_zorder,
    )

    if np.any(nonzero_mask):
        cbar = fig.colorbar(sc, ax=ax, shrink=0.4)
        cbar.set_label(f"{gene_bin_data.gene_name} transcript count")

    for overlay in section_overlays:
        if overlay.analysis_region is not None:
            _plot_geometry_boundary(
                ax=ax,
                geometry=overlay.analysis_region,
                edgecolor="#00BFFF",
                linewidth=shared_line_width,
                zorder=highlighted_analysis_region_zorder,
            )
        elif (
            overlay.analysis_bin_centers is not None
            and np.size(overlay.analysis_bin_centers) > 0
        ):
            _plot_bin_region_boundary(
                ax=ax,
                bin_centers=overlay.analysis_bin_centers,
                edgecolor="#00BFFF",
                linewidth=shared_line_width,
                zorder=highlighted_analysis_region_zorder,
            )

    if show_axes:
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labelleft=False,
        )
    if bounds is None:
        bounds = full_tissue_plot_bounds(
            full_tissue_contour, gene_bin_data.bin_centers[nonzero_mask], section_overlays,
        )
    xmin, ymin, xmax, ymax = bounds
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    bar_um = 1000.0
    bar_px = bar_um / pixel_size_um
    add_scalebar(
        ax=ax,
        length_px=bar_px,
        label="",
        loc="lower left",
        frameon=False,
        size_vertical_scale=18.0,
    )

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
    _plot_start = perf_counter()
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Plot bin centers colored by transcript count (same style as linearized-distance plot).
    counts = np.asarray(data.gene_bins.transcript_counts, dtype=float)
    positive_counts = counts[counts > 0]
    norm = None
    if positive_counts.size > 0:
        norm = LogNorm(vmin=float(np.min(positive_counts)), vmax=float(np.max(positive_counts)))

    # Draw low-count bins first and high-count bins last so high counts appear on top.
    draw_order = np.argsort(counts, kind="stable")
    centers_sorted = data.gene_bins.bin_centers[draw_order]
    counts_sorted = counts[draw_order]
    closest_xy_sorted = data.distance_data.closest_xy[draw_order]
    distances_sorted = np.asarray(data.distance_data.closest_distances, dtype=float)[draw_order]

    cmap = plt.get_cmap("viridis")
    if norm is not None:
        color_norm = norm
    elif counts.size > 0:
        vmin = float(np.min(counts))
        vmax = float(np.max(counts))
        if np.isclose(vmin, vmax):
            vmax = vmin + 1.0
        color_norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        color_norm = Normalize(vmin=0.0, vmax=1.0)
    line_colors = cmap(color_norm(counts_sorted))

    sm = plt.cm.ScalarMappable(norm=color_norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(f"{data.gene_bins.gene_name} transcript count")

    # Draw each bin pair in order: line first, then the matching dot.
    # This gives true interleaving depth control for each bin/line pair.
    for i, (center_xy, closest_xy, distance, color) in enumerate(
        zip(centers_sorted, closest_xy_sorted, distances_sorted, line_colors)
    ):
        line_z = 2 * i
        dot_z = line_z + 1
        line_style = ":" if distance < 0 else "-"
        line_kwargs = {}
        if distance < 0:
            line_kwargs["dash_capstyle"] = "round"

        ax.plot(
            [center_xy[0], closest_xy[0]],
            [center_xy[1], closest_xy[1]],
            linestyle=line_style,
            color=color,
            linewidth=line_width,
            zorder=line_z,
            **line_kwargs,
        )
        ax.scatter(
            center_xy[0],
            center_xy[1],
            c=[color],
            s=1,
            zorder=dot_z,
        )

    # Plot the contour
    line_xy = np.array(contour.coords)
    ax.plot(line_xy[:, 0], line_xy[:, 1], "-", c="red", lw=line_width, label="Boundary", zorder=1_000_000)

    ax.set_aspect("equal")
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.set_title(title)
    ax.invert_yaxis()

    save_or_show_figure(fig, output_path)
    print(f"plot_bins_and_closest_points_on_contour took {perf_counter() - _plot_start:.1f}s")


def plot_distances_linearized(
    data: "DiffusionBinData",
    output_path: str,
    title: str,
    pixel_size_um: float,
    figsize: tuple[float, float] = (10, 2.5),
    dpi: int | None = None,
) -> None:
    # Scatter plot of bins: along-boundary position vs distance to boundary.
    x = data.distance_data.closest_along_boundary
    y = np.asarray(data.distance_data.closest_distances, dtype=float) * pixel_size_um
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
        zorder=2,
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label(f"{data.gene_bins.gene_name} transcript count")
    ax.plot([0.0, 1.0], [0.0, 0.0], color="red", linewidth=1, linestyle="-", zorder=0)
    ax.set_xlabel("Normalized position along boundary")
    ax.set_ylabel(f"Distance to boundary [{MICROMETER_UNIT}]")

    ax.set_title(title)

    save_or_show_figure(fig, output_path)


def plot_distance_histogram(
    data: DistanceDistribution,
    output_path: str,
    gene_name: str,
    max_distance_um: float,
    title: str | None = None,
    bin_width_um: float = 5.0,
) -> None:
    kde_color = "#1f77b4"
    distances = np.asarray(data.distances_um, dtype=float)
    weights = np.asarray(data.weights, dtype=float)
    distances, weights = _sanitize_distance_weights(distances, weights, label="distance histogram")

    fig, ax = plt.subplots(figsize=(6, 6))
    bins = _distance_hist_bins(
        [DistanceDistribution("distance histogram", distances, weights)],
        bin_width=bin_width_um,
    )
    bin_centers = _distance_hist_bin_centers_from_edges(bins)
    bin_widths = np.diff(bins)
    histogram_counts, _ = np.histogram(distances, bins=bins, weights=weights, density=False)
    density_curve = _distance_kde_curves(
        [DistanceDistribution("distance histogram", distances, weights)],
        align_peaks=False,
    )[0]

    bar_container = ax.bar(
        bin_centers,
        histogram_counts,
        width=bin_widths,
        align="center",
        color="#d9d9d9",
        edgecolor="#666666",
        linewidth=0.8,
        alpha=0.9,
        label="Transcript-weighted histogram",
    )
    ax.axvline(
        density_curve.peak_position_um,
        color=kde_color,
        linewidth=1.2,
        linestyle=":",
        zorder=4,
    )
    # Add inside tissue/outside tissue labels on the histogram plots
    label_bbox = {
        "boxstyle": "square,pad=0.25",
        "facecolor": "white",
        "edgecolor": "black",
        "linewidth": 1.0,
    }
    ax.text(
        0.25,
        0.95,
        "Outside tissue",
        transform=ax.transAxes,
        ha="center",
        va="top",
        bbox=label_bbox,
        zorder=4,
    )
    ax.text(
        0.75,
        0.95,
        "Inside tissue",
        transform=ax.transAxes,
        ha="center",
        va="top",
        bbox=label_bbox,
        zorder=4,
    )
    ax.set_xlabel(f"Distance to boundary [{MICROMETER_UNIT}]")
    ax.set_ylabel(f"{gene_name} transcript count")
    ax.set_title(title)
    density_ax = ax.twinx()
    kde_line, = density_ax.plot(
        density_curve.x_um,
        density_curve.density,
        color=kde_color,
        linewidth=1.2,
        label="Transcript-weighted KDE",
    )
    peak_marker = density_ax.scatter(
        [density_curve.peak_position_um],
        [density_curve.peak_density],
        s=28,
        color=kde_color,
        edgecolors="white",
        linewidths=0.8,
        zorder=5,
        label=f"Detected peak (d={density_curve.peak_position_um:.2f} {MICROMETER_UNIT})",
    )
    if max_distance_um <= 0:
        raise ValueError("max_distance_um must be positive.")
    ax.set_xlim(-max_distance_um, max_distance_um)
    density_ax.set_ylabel(f"{gene_name} weighted KDE density")
    density_ax.tick_params(right=True, labelright=True)
    density_ax.spines["right"].set_visible(True)

    handles = [bar_container, kde_line, peak_marker]
    labels = [handle.get_label() for handle in handles]
    ax.legend(handles, labels, fontsize="small", loc="upper left", bbox_to_anchor=(1.16, 1.0), borderaxespad=0.0)
    save_or_show_figure(fig, output_path)


def plot_distance_kde_combined(
    crop_specs: list[tuple[str, str, float, str]],
    output_path: str,
    plots_root: str,
    gene_name: str,
    max_distance_um: float | None = None,
    title: str | None = None,
    figsize = (12, 6),
    align_peaks: bool = False,
    bin_width_um: float = 5.0,
) -> None:
    """
    Plot one normalized signed-distance KDE per experiment, pooling its crops.

    If `align_peaks` is True, each KDE curve is shifted so its detected density peak
    is placed at x=0.
    Each entry in `crop_specs` is `(technology, dataset_id, pixel_size_um, crop_name)`.
    Histograms are computed in micrometers with shared bins across experiments, and
    the overlaid density curves are true weighted KDEs evaluated on the raw
    distance values.

    Args:
        crop_specs: Crop descriptors whose signed distances/counts are pooled per experiment.
        output_path: Destination path for the rendered figure.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.
        gene_name: Name used in the y-axis label.
        max_distance_um: Optional symmetric x-limit in micrometers.
        title: Optional plot title.
        figsize: Matplotlib figure size.
        align_peaks: Whether to shift each KDE curve so its mode is at zero.
        bin_width_um: Histogram bin width in micrometers.

    Raises:
        ValueError: If `max_distance_um` is provided but not positive.
    """
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    bins = _distance_hist_bins(
        data_per_label,
        bin_width=bin_width_um,
    )
    curves = _distance_kde_curves(
        data_per_label,
        align_peaks=align_peaks,
    )

    fig, ax = plt.subplots(figsize=figsize)
    for distribution, curve in zip(data_per_label, curves, strict=True):
        color = COLOR_PALETTE.get(distribution.technology)
        label_distances_um = distribution.distances_um
        label_weights = distribution.weights
        histogram_density, _ = np.histogram(label_distances_um, bins=bins, weights=label_weights, density=True)
        histogram_centers_um = _distance_hist_bin_centers_from_edges(bins)
        if align_peaks:
            histogram_centers_um = histogram_centers_um - curve.alignment_shift_um
        ax.step(
            histogram_centers_um,
            histogram_density,
            where="mid",
            linewidth=0.9,
            alpha=0.3,
            color=color,
        )
        ax.plot(
            curve.x_um,
            curve.density,
            linewidth=1.8,
            label=curve.label,
            color=color,
        )
        ax.scatter(
            [curve.peak_position_um],
            [curve.peak_density],
            s=24,
            color=color,
            edgecolors="white",
            linewidths=0.8,
            zorder=4,
        )

    ax.axvline(0, color="black", linewidth=1, linestyle=":", zorder=3)
    x_label = (
        f"Peak-aligned distance to boundary [{MICROMETER_UNIT}]"
        if align_peaks
        else f"Distance to boundary [{MICROMETER_UNIT}]"
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(f"Normalized {gene_name} expression")
    ax.set_title(title)
    ax.legend(fontsize="small")
    if max_distance_um is not None:
        if max_distance_um <= 0:
            raise ValueError("max_distance_um must be positive.")
        ax.set_xlim(-max_distance_um, max_distance_um)

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


def read_experiment_distance_data(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
) -> list[DistanceDistribution]:
    """Pool configured crops by (technology, dataset_id), in first-seen order.

    Sections are assumed to select disjoint bins: choose boundary sections far
    enough apart for the analysis distance limit. No overlap detection or
    deduplication is performed. Concatenate signed distances and transcript counts
    before normalization, KDE peak detection/alignment, or summary statistics.
    Only explicitly listed crops are read; leftover crop folders are ignored.
    """
    if not crop_specs:
        raise ValueError("crop_specs must contain at least one entry.")

    grouped: dict[tuple[str, str], list[tuple[np.ndarray, np.ndarray]]] = {}
    for technology, dataset_id, pixel_size_um, crop_name in crop_specs:
        chunks = grouped.setdefault((technology, dataset_id), [])
        npz_path = os.path.join(
            plots_root, technology, dataset_id, "analysis", crop_name, "diffusion.npz",
        )
        if not os.path.exists(npz_path):
            raise FileNotFoundError(
                f"Missing diffusion results file: {npz_path}. Run compute_diffusion() first."
            )
        with np.load(npz_path) as data:
            distances_px = np.asarray(data["filtered_closest_distances"], dtype=float)
            counts = np.asarray(data["filtered_transcript_counts"], dtype=float)
        if distances_px.ndim != 1 or counts.shape != distances_px.shape:
            raise ValueError(f"Distance and transcript-count arrays must be matching 1D arrays: {npz_path}")
        if distances_px.size:
            chunks.append((distances_px * pixel_size_um, counts))

    technology_counts: dict[str, int] = {}
    for technology, _ in grouped:
        technology_counts[technology] = technology_counts.get(technology, 0) + 1

    distributions: list[DistanceDistribution] = []
    for (technology, dataset_id), chunks in grouped.items():
        if not chunks:
            raise ValueError(f"No distance rows found for experiment {technology!r}, {dataset_id!r}.")
        label = technology if technology_counts[technology] == 1 else f"{technology} ({dataset_id})"
        distributions.append(DistanceDistribution(
            label=label,
            distances_um=np.concatenate([distances for distances, _ in chunks]),
            weights=np.concatenate([counts for _, counts in chunks]),
            technology=technology,
            dataset_id=dataset_id,
        ))
    return distributions


def _distance_hist_bins(
    data_per_label: list[DistanceDistribution],
    bin_width: float = 5.0,
) -> np.ndarray:
    if not np.isfinite(bin_width) or bin_width <= 0:
        raise ValueError("bin_width must be a positive finite value.")
    concat = np.concatenate([distribution.distances_um for distribution in data_per_label])
    max_abs = float(np.max(np.abs(concat)))
    if not np.isfinite(max_abs):
        raise ValueError("Histogram distances must contain finite values.")
    half_bin = 0.5 * bin_width
    n_half = int(np.ceil(max_abs / bin_width))
    centers = np.arange(-n_half, n_half + 1, dtype=float) * bin_width
    return np.concatenate([centers - half_bin, [centers[-1] + half_bin]])


def _distance_hist_bin_centers_from_edges(bins: np.ndarray) -> np.ndarray:
    bins = np.asarray(bins, dtype=float)
    if bins.ndim != 1 or bins.size < 2:
        raise ValueError("bins must be a 1D array with at least two edges.")
    return 0.5 * (bins[:-1] + bins[1:])


def _distance_kde_grid(
    max_abs_um: float,
    max_bandwidth_um: float = 0.0,
    step_um: float = KDE_GRID_STEP_UM,
    tail_bandwidths: float = KDE_TAIL_BANDWIDTHS,
) -> np.ndarray:
    if not np.isfinite(max_abs_um) or max_abs_um < 0:
        raise ValueError("max_abs_um must be finite and non-negative.")
    if not np.isfinite(step_um) or step_um <= 0:
        raise ValueError("step_um must be finite and positive.")

    limit_um = max_abs_um + tail_bandwidths * max(0.0, max_bandwidth_um)
    n_half = max(1, int(np.ceil(limit_um / step_um)))
    return np.arange(-n_half, n_half + 1, dtype=float) * step_um


def _weighted_mean_and_std(
    values: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0:
        raise ValueError("weights must sum to a positive finite value.")
    mean = float(np.sum(weights * values) / total_weight)
    variance = float(np.sum(weights * (values - mean) ** 2) / total_weight)
    variance = max(variance, 0.0)
    return mean, float(np.sqrt(variance))


def _effective_sample_size(weights: np.ndarray) -> float:
    total_weight = float(np.sum(weights))
    sumsq_weight = float(np.sum(weights * weights))
    if not np.isfinite(total_weight) or not np.isfinite(sumsq_weight) or total_weight <= 0 or sumsq_weight <= 0:
        raise ValueError("weights must be positive and finite.")
    return (total_weight * total_weight) / sumsq_weight


def _weighted_kde_bandwidth_um(
    values_um: np.ndarray,
    weights: np.ndarray,
    fallback_scale_um: float = KDE_BANDWIDTH_FALLBACK_UM,
) -> float:
    _, weighted_std = _weighted_mean_and_std(values_um, weights)
    q25_um, q75_um = _weighted_quantiles(values_um, weights, probs=(0.25, 0.75))
    weighted_iqr_um = float(q75_um - q25_um)
    robust_scale_um = weighted_iqr_um / 1.349 if weighted_iqr_um > 0 else np.nan

    scale_candidates = [weighted_std]
    if np.isfinite(robust_scale_um) and robust_scale_um > 0:
        scale_candidates.append(robust_scale_um)
    positive_scales = [candidate for candidate in scale_candidates if np.isfinite(candidate) and candidate > 0]
    scale_um = min(positive_scales) if positive_scales else np.nan

    n_eff = _effective_sample_size(weights)
    if np.isfinite(scale_um) and scale_um > 0 and np.isfinite(n_eff) and n_eff > 1:
        bandwidth_um = 0.9 * scale_um * (n_eff ** (-1.0 / 5.0))
        if np.isfinite(bandwidth_um) and bandwidth_um > 0:
            return float(bandwidth_um)

    data_range_um = float(np.max(values_um) - np.min(values_um))
    fallback_candidates = [fallback_scale_um, data_range_um / 20.0, 1e-6]
    return float(max(candidate for candidate in fallback_candidates if np.isfinite(candidate) and candidate > 0))


def _shared_kde_bandwidth_um(
    data_per_label: list[DistanceDistribution],
) -> float:
    if len(data_per_label) == 0:
        raise ValueError("data_per_label must contain at least one dataset.")

    pooled_distances: list[np.ndarray] = []
    pooled_weights: list[np.ndarray] = []
    for distribution in data_per_label:
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        pooled_distances.append(distances_um)
        pooled_weights.append(weights)

    return _weighted_kde_bandwidth_um(
        np.concatenate(pooled_distances),
        np.concatenate(pooled_weights),
    )


def _evaluate_weighted_gaussian_kde(
    x_um: np.ndarray,
    values_um: np.ndarray,
    weights: np.ndarray,
    bandwidth_um: float,
    value_chunk_size: int = 4096,
) -> np.ndarray:
    if bandwidth_um <= 0 or not np.isfinite(bandwidth_um):
        raise ValueError("bandwidth_um must be finite and positive.")
    total_weight = float(np.sum(weights))
    if total_weight <= 0 or not np.isfinite(total_weight):
        raise ValueError("weights must sum to a positive finite value.")

    density = np.zeros_like(x_um, dtype=float)
    normalization = total_weight * bandwidth_um * np.sqrt(2.0 * np.pi)
    for start in range(0, values_um.size, value_chunk_size):
        stop = min(start + value_chunk_size, values_um.size)
        value_chunk = values_um[start:stop]
        weight_chunk = weights[start:stop]
        scaled_delta = (x_um[:, None] - value_chunk[None, :]) / bandwidth_um
        density += np.exp(-0.5 * scaled_delta * scaled_delta) @ weight_chunk
    return density / normalization


def _detect_peak_from_curve(
    x_um: np.ndarray,
    density_values: np.ndarray,
) -> tuple[float, float]:
    if x_um.ndim != 1 or density_values.ndim != 1 or x_um.size != density_values.size:
        raise ValueError("x_um and density_values must be 1D arrays with equal length.")
    if x_um.size == 0:
        raise ValueError("Density curve must contain at least one sample.")

    peak_density = float(np.max(density_values))
    if not np.isfinite(peak_density) or peak_density <= 0:
        raise ValueError("Density peak must be finite and positive.")
    peak_indices = np.flatnonzero(np.isclose(density_values, peak_density))
    if peak_indices.size == 0:
        raise ValueError("No peak index found in density curve.")
    peak_position_um = float(np.mean(x_um[peak_indices]))
    return peak_position_um, peak_density


def _distance_kde_curves(
    data_per_label: list[DistanceDistribution],
    align_peaks: bool = False,
    x_um: np.ndarray | None = None,
) -> list[DistanceKDECurve]:
    sanitized_data: list[DistanceDistribution] = []
    max_abs_um = 0.0
    for distribution in data_per_label:
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        sanitized_data.append(replace(distribution, distances_um=distances_um, weights=weights))
        max_abs_um = max(max_abs_um, float(np.max(np.abs(distances_um))))

    shared_bandwidth_um = _shared_kde_bandwidth_um(sanitized_data)
    if len(sanitized_data) > 1:
        labels_str = ", ".join(distribution.label for distribution in sanitized_data)
        print(f"Using shared KDE bandwidth: {shared_bandwidth_um:.3f} {MICROMETER_UNIT} for {labels_str}")
    bandwidths_um = [shared_bandwidth_um] * len(sanitized_data)

    if x_um is None:
        x_um = _distance_kde_grid(
            max_abs_um=max_abs_um,
            max_bandwidth_um=max(bandwidths_um, default=0.0),
        )
    x_um = np.asarray(x_um, dtype=float)

    curves: list[DistanceKDECurve] = []
    for distribution, bandwidth_um in zip(sanitized_data, bandwidths_um, strict=True):
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        density_values = _evaluate_weighted_gaussian_kde(x_um, distances_um, weights, bandwidth_um=bandwidth_um)
        peak_position_um, peak_density = _detect_peak_from_curve(x_um, density_values)
        alignment_shift_um = peak_position_um
        curve_x_um = x_um
        if align_peaks:
            curve_x_um = x_um - peak_position_um
            peak_position_um = 0.0
        curves.append(
            DistanceKDECurve(
                label=label,
                x_um=curve_x_um,
                density=density_values,
                peak_position_um=peak_position_um,
                peak_density=peak_density,
                alignment_shift_um=alignment_shift_um,
                technology=distribution.technology,
            )
        )
    return curves


def _distance_hist_curves(
    data_per_label: list[DistanceDistribution],
    bins: np.ndarray,
    align_peaks: bool = False,
    density: bool = False,
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    curves: list[tuple[str, np.ndarray, np.ndarray]] = []
    for distribution in data_per_label:
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        counts, edges = np.histogram(distances_um, bins=bins, weights=weights, density=density)
        centers = _distance_hist_bin_centers_from_edges(edges)
        if align_peaks:
            peak_idx = int(np.argmax(counts))
            centers = centers - centers[peak_idx]
        curves.append((label, centers, counts))
    return curves


def _distance_kde_curves_from_crop_specs(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
    align_peaks: bool = False,
) -> list[DistanceKDECurve]:
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    return _distance_kde_curves(
        data_per_label,
        align_peaks=align_peaks,
    )


def _distance_hist_curves_from_crop_specs(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
    bin_width_um: float = 5.0,
    align_peaks: bool = False,
    density: bool = True,
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    bins = _distance_hist_bins(
        data_per_label,
        bin_width=bin_width_um,
    )
    return _distance_hist_curves(data_per_label, bins=bins, align_peaks=align_peaks, density=density)


def _linear_interpolate_x_at_y(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    target_y: float,
) -> float:
    if np.isclose(y1, y0):
        return 0.5 * (x0 + x1)
    t = (target_y - y0) / (y1 - y0)
    return x0 + t * (x1 - x0)


def _find_hwhm_from_curve(
    x_um: np.ndarray,
    y_values: np.ndarray,
) -> tuple[float, float, float]:
    """
    Compute half-width at half-maximum (HWHM) on both sides of a 1D peak curve.

    The peak height is defined as the global maximum of `y_values`. Half-maximum is
    `0.5 * peak_count`. Peak samples are identified as all x-locations whose value is
    `np.isclose` to the global maximum.

    If multiple samples share the maximum:
    - `peak_left_idx` is the first (leftmost) max index.
    - `peak_right_idx` is the last (rightmost) max index.
    - `peak_x` is the mean of the corresponding x positions.

    Left and right half-maximum crossing positions are found by scanning outward
    from the peak edges and taking the first interval that brackets `half_max`.
    Crossing x-coordinates are linearly interpolated between neighboring x positions.
    If no crossing exists on a side, that side's HWHM is returned as `np.nan`.

    Parameters
    ----------
    x_um : np.ndarray
        1D array of x positions.
    y_values : np.ndarray
        1D array of curve values at `x_um` (same length as `x_um`).

    Returns
    -------
    tuple[float, float, float]
        `(peak_x, left_hwhm, right_hwhm)`, where:
        - `peak_x` is the peak x-position as defined above,
        - `left_hwhm` is `peak_x - left_halfmax_x` (or `np.nan` if not found),
        - `right_hwhm` is `right_halfmax_x - peak_x` (or `np.nan` if not found).
    """
    if x_um.ndim != 1 or y_values.ndim != 1 or x_um.size != y_values.size:
        raise ValueError("x_um and y_values must be 1D arrays with equal length.")
    if x_um.size == 0:
        raise ValueError("Curve must contain at least one sample.")

    peak_count = float(np.max(y_values))
    if not np.isfinite(peak_count) or peak_count <= 0:
        raise ValueError("Curve peak must be finite and positive.")
    half_max = 0.5 * peak_count

    peak_indices = np.flatnonzero(np.isclose(y_values, peak_count))
    if peak_indices.size == 0:
        raise ValueError("No peak index found in curve.")
    peak_left_idx = int(peak_indices[0])
    peak_right_idx = int(peak_indices[-1])
    peak_x = float(np.mean(x_um[peak_indices]))

    left_x = np.nan
    for i in range(peak_left_idx - 1, -1, -1):
        if y_values[i] <= half_max <= y_values[i + 1]:
            left_x = float(
                _linear_interpolate_x_at_y(
                    float(x_um[i]),
                    float(y_values[i]),
                    float(x_um[i + 1]),
                    float(y_values[i + 1]),
                    half_max,
                )
            )
            break

    right_x = np.nan
    for i in range(peak_right_idx + 1, y_values.size):
        if y_values[i] <= half_max <= y_values[i - 1]:
            right_x = float(
                _linear_interpolate_x_at_y(
                    float(x_um[i - 1]),
                    float(y_values[i - 1]),
                    float(x_um[i]),
                    float(y_values[i]),
                    half_max,
                )
            )
            break

    left_hwhm = float(peak_x - left_x) if np.isfinite(left_x) else np.nan
    right_hwhm = float(right_x - peak_x) if np.isfinite(right_x) else np.nan
    return peak_x, left_hwhm, right_hwhm


def calculate_distance_histogram_hwhm(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
    align_peaks: bool = False,
) -> list[DistanceHistogramHWHM]:
    """
    Calculate left/right HWHM values from pooled experiment signed-distance KDE curves.

    The left HWHM corresponds to diffusion towards the outside of the tissue (negative distances),
    the right HWHM towards the inside of the tissue (positive distances).
    Each entry in `crop_specs` is `(technology, dataset_id, pixel_size_um, crop_name)`.
    KDE curves are evaluated on a shared symmetric x-grid anchored at zero,
    consistent with the combined and ridgeline visualizations.

    Args:
        crop_specs: Crop descriptors whose signed distances/counts are pooled per experiment.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.
        align_peaks: Whether to shift each KDE curve so its mode is at zero before
            HWHM extraction.

    Returns:
        List of `DistanceHistogramHWHM` rows (one row per experiment) containing
        peak position and left/right HWHM in micrometers.
    """
    rows: list[DistanceHistogramHWHM] = []
    curves = _distance_kde_curves_from_crop_specs(
        crop_specs=crop_specs,
        plots_root=plots_root,
        align_peaks=align_peaks,
    )

    for curve in curves:
        peak_x, left_hwhm, right_hwhm = _find_hwhm_from_curve(curve.x_um, curve.density)
        rows.append(
            DistanceHistogramHWHM(
                label=curve.label,
                left_hwhm_um=left_hwhm,
                right_hwhm_um=right_hwhm,
                peak_position_um=peak_x,
                technology=curve.technology,
            )
        )

    return rows


def plot_distance_histogram_hwhm_bars(
    hwhm_rows: list[DistanceHistogramHWHM],
    output_path: str,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
) -> None:
    """
    Plot left/right HWHM per experiment as bidirectional horizontal bars.

    Left HWHM is drawn towards negative x; right HWHM towards positive x.

    Args:
        hwhm_rows: HWHM summary rows produced by `calculate_distance_histogram_hwhm`.
        output_path: Destination path for the rendered figure.
        title: Optional plot title.
        figsize: Optional figure size. If omitted, a height is chosen from row count.

    Raises:
        ValueError: If input is empty or no rows contain finite non-negative left/right HWHM.
    """
    if len(hwhm_rows) == 0:
        raise ValueError("hwhm_rows must contain at least one row.")

    labels = [row.label for row in hwhm_rows]
    left_hwhm = np.asarray([row.left_hwhm_um for row in hwhm_rows], dtype=float)
    right_hwhm = np.asarray([row.right_hwhm_um for row in hwhm_rows], dtype=float)
    valid = np.isfinite(left_hwhm) & np.isfinite(right_hwhm) & (left_hwhm >= 0) & (right_hwhm >= 0)
    if not np.any(valid):
        raise ValueError("No finite non-negative left/right HWHM pairs available to plot.")

    dropped = int(np.sum(~valid))
    if dropped > 0:
        print(f"Dropping {dropped} rows with invalid HWHM values before plotting.")

    keep_idx = np.flatnonzero(valid)
    labels = [labels[i] for i in keep_idx]
    left_hwhm = left_hwhm[keep_idx]
    right_hwhm = right_hwhm[keep_idx]

    if figsize is None:
        figsize = (12, max(4.5, 0.45 * len(labels) + 1.5))

    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=figsize)
    tech_colors = [
        COLOR_PALETTE.get(hwhm_rows[i].technology or hwhm_rows[i].label, "#7f7f7f")
        for i in keep_idx
    ]
    missing_color_labels = [label for label, color in zip(labels, tech_colors) if color == "#7f7f7f"]
    assert len(missing_color_labels) == 0

    # Keep technology identity in color; use hatch style to encode inside/outside tissue.
    ax.barh(
        y,
        -left_hwhm,
        color=tech_colors,
        alpha=0.9,
        hatch="///",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.barh(
        y,
        right_hwhm,
        color=tech_colors,
        alpha=0.9,
        hatch="",
        edgecolor="black",
        linewidth=0.5,
    )

    max_abs = float(np.max(np.concatenate([left_hwhm, right_hwhm])))
    margin = max(0.5, 0.1 * max_abs)
    ax.set_xlim(-(max_abs + margin), max_abs + margin)
    ax.axvline(0.0, color="black", linewidth=1, linestyle=":", zorder=3)
    ax.grid(axis="x", alpha=0.25)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"HWHM [{MICROMETER_UNIT}]")
    ax.set_title(title or "Distance histogram half-width at half-maximum")
    side_legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="///", linewidth=0.8),
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="", linewidth=0.8),
    ]
    ax.legend(
        side_legend_handles,
        ["Outside tissue", "Inside tissue"],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )

    save_or_show_figure(fig, output_path)


def calculate_distance_to_kde_peak_median(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
) -> list[DistanceKDEPeakMedianDistance]:
    """
    Calculate weighted median distances to the detected KDE peak on both sides.

    Distances are first shifted so the pooled experiment KDE peak is at zero. The left
    median summarizes distances on the negative side of that peak and the right
    median summarizes distances on the positive side.

    Args:
        crop_specs: Crop descriptors whose signed distances/counts are pooled per experiment.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.

    Returns:
        One row per experiment containing the detected peak position together with
        the weighted median distance to that peak on the left and right side.
    """
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    curves = _distance_kde_curves(data_per_label, align_peaks=False)

    rows: list[DistanceKDEPeakMedianDistance] = []
    for distribution, curve in zip(data_per_label, curves, strict=True):
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        if label != curve.label:
            raise ValueError(f"Label mismatch while computing KDE-peak medians: {label!r} != {curve.label!r}")

        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        left_distances_um, left_weights, right_distances_um, right_weights = _split_peak_aligned_distance_sides(
            distances_um,
            weights,
            peak_position_um=curve.peak_position_um,
        )

        left_median_distance_um = (
            float(_weighted_quantiles(left_distances_um, left_weights, probs=(0.50,))[0])
            if left_distances_um.size > 0
            else np.nan
        )
        right_median_distance_um = (
            float(_weighted_quantiles(right_distances_um, right_weights, probs=(0.50,))[0])
            if right_distances_um.size > 0
            else np.nan
        )

        rows.append(
            DistanceKDEPeakMedianDistance(
                label=label,
                left_median_distance_um=left_median_distance_um,
                right_median_distance_um=right_median_distance_um,
                peak_position_um=curve.peak_position_um,
                technology=distribution.technology,
            )
        )

    return rows


def plot_distance_to_kde_peak_median_bars(
    median_rows: list[DistanceKDEPeakMedianDistance],
    output_path: str,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
) -> None:
    """
    Plot left/right weighted median distance to the detected KDE peak per experiment.

    Left-of-peak medians (~outside tissue) are drawn towards negative x;
    right-of-peak medians (~inside tissue) are drawn towards positive x.

    Args:
        median_rows: Summary rows produced by `calculate_distance_to_kde_peak_median`.
        output_path: Destination path for the rendered figure.
        title: Optional plot title.
        figsize: Optional figure size. If omitted, a height is chosen from row count.

    Raises:
        ValueError: If input is empty or no rows contain finite non-negative left/right medians.
    """
    if len(median_rows) == 0:
        raise ValueError("median_rows must contain at least one row.")

    labels = [row.label for row in median_rows]
    left_median = np.asarray([row.left_median_distance_um for row in median_rows], dtype=float)
    right_median = np.asarray([row.right_median_distance_um for row in median_rows], dtype=float)
    valid = np.isfinite(left_median) & np.isfinite(right_median) & (left_median >= 0) & (right_median >= 0)
    if not np.any(valid):
        raise ValueError("No finite non-negative left/right KDE-peak medians available to plot.")

    dropped = int(np.sum(~valid))
    if dropped > 0:
        print(f"Dropping {dropped} rows with invalid KDE-peak median values before plotting.")

    keep_idx = np.flatnonzero(valid)
    labels = [labels[i] for i in keep_idx]
    left_median = left_median[keep_idx]
    right_median = right_median[keep_idx]

    if figsize is None:
        figsize = (12, max(4.5, 0.45 * len(labels) + 1.5))

    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=figsize)
    tech_colors = [
        COLOR_PALETTE.get(median_rows[i].technology or median_rows[i].label, "#7f7f7f")
        for i in keep_idx
    ]
    missing_color_labels = [label for label, color in zip(labels, tech_colors) if color == "#7f7f7f"]
    assert len(missing_color_labels) == 0

    ax.barh(
        y,
        -left_median,
        color=tech_colors,
        alpha=0.9,
        hatch="///",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.barh(
        y,
        right_median,
        color=tech_colors,
        alpha=0.9,
        hatch="",
        edgecolor="black",
        linewidth=0.5,
    )

    max_abs = float(np.max(np.concatenate([left_median, right_median])))
    margin = max(0.5, 0.1 * max_abs)
    ax.set_xlim(-(max_abs + margin), max_abs + margin)
    ax.axvline(0.0, color="black", linewidth=1, linestyle=":", zorder=3)
    ax.grid(axis="x", alpha=0.25)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"Median distance [{MICROMETER_UNIT}]")
    ax.set_title(title or "Median distance to KDE peak")
    side_legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="///", linewidth=0.8),
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="", linewidth=0.8),
    ]
    ax.legend(
        side_legend_handles,
        ["Outside tissue", "Inside tissue"], # Left of KDE peak, Right of KDE peak
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )

    save_or_show_figure(fig, output_path)


def _weighted_quantiles(
    values: np.ndarray,
    weights: np.ndarray,
    probs: Sequence[float],
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    probs = np.asarray(probs, dtype=float)
    if values.ndim != 1 or weights.ndim != 1 or values.size != weights.size:
        raise ValueError("values and weights must be 1D arrays of equal length.")
    if probs.ndim != 1 or probs.size == 0:
        raise ValueError("probs must be a non-empty 1D sequence.")
    if np.any(probs < 0) or np.any(probs > 1):
        raise ValueError("All probabilities must be in [0, 1].")

    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    wsum = float(np.sum(w))
    if not np.isfinite(wsum) or wsum <= 0:
        raise ValueError("weights must sum to a positive finite value.")
    cdf = np.cumsum(w) / wsum
    return np.interp(probs, cdf, v, left=v[0], right=v[-1])


def _split_signed_distance_sides(
    distances_um: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    outside_mask = distances_um <= 0
    inside_mask = distances_um > 0
    outside_distances_um = -distances_um[outside_mask]
    outside_weights = weights[outside_mask]
    inside_distances_um = distances_um[inside_mask]
    inside_weights = weights[inside_mask]
    return outside_distances_um, outside_weights, inside_distances_um, inside_weights


def _split_peak_aligned_distance_sides(
    distances_um: np.ndarray,
    weights: np.ndarray,
    peak_position_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    shifted_distances_um = np.asarray(distances_um, dtype=float) - float(peak_position_um)
    left_mask = shifted_distances_um <= 0
    right_mask = shifted_distances_um > 0
    left_distances_um = -shifted_distances_um[left_mask]
    left_weights = weights[left_mask]
    right_distances_um = shifted_distances_um[right_mask]
    right_weights = weights[right_mask]
    return left_distances_um, left_weights, right_distances_um, right_weights


def summarize_distance_distributions(
    crop_specs: list[tuple[str, str, float, str]],
    plots_root: str,
    quantiles: Sequence[float] = (0.10, 0.25, 0.50, 0.75, 0.90),
    cdf_levels: Sequence[float] = (0.25, 0.50, 0.75, 0.90),
) -> list[DistanceDistributionSummary]:
    """
    Summarize outside/inside signed-distance distributions for each experiment, pooling its crops.

    Distances are split by sign into outside (`d <= 0`, reported as absolute distance)
    and inside (`d > 0`) components, both weighted by transcript counts.

    Args:
        crop_specs: Crop descriptors as `(technology, dataset_id, pixel_size_um, crop_name)`.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.
        quantiles: Probability levels used to report weighted distance quantiles per side.
        cdf_levels: CDF levels used to report distance-at-CDF summaries per side.
    Returns:
        One `DistanceDistributionSummary` per experiment with quantile signatures,
        distance-at-CDF summaries, and side mass fractions.
    """
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    quantiles = tuple(float(q) for q in quantiles)
    cdf_levels = tuple(float(q) for q in cdf_levels)

    rows: list[DistanceDistributionSummary] = []
    for distribution in data_per_label:
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        outside_distances_um, outside_weights, inside_distances_um, inside_weights = _split_signed_distance_sides(
            distances_um, weights
        )

        total = float(np.sum(weights))
        outside_total = float(np.sum(outside_weights))
        inside_total = float(np.sum(inside_weights))

        quantile_map: dict[str, float] = {}
        if outside_distances_um.size > 0:
            outside_quantile_values = _weighted_quantiles(outside_distances_um, outside_weights, quantiles)
            quantile_map.update(
                {
                    f"outside_q{int(round(q * 100)):02d}_um": float(val)
                    for q, val in zip(quantiles, outside_quantile_values, strict=True)
                }
            )
        else:
            quantile_map.update({f"outside_q{int(round(q * 100)):02d}_um": np.nan for q in quantiles})
        if inside_distances_um.size > 0:
            inside_quantile_values = _weighted_quantiles(inside_distances_um, inside_weights, quantiles)
            quantile_map.update(
                {
                    f"inside_q{int(round(q * 100)):02d}_um": float(val)
                    for q, val in zip(quantiles, inside_quantile_values, strict=True)
                }
            )
        else:
            quantile_map.update({f"inside_q{int(round(q * 100)):02d}_um": np.nan for q in quantiles})

        cdf_map: dict[str, float] = {}
        if outside_distances_um.size > 0:
            outside_cdf_values = _weighted_quantiles(outside_distances_um, outside_weights, cdf_levels)
            cdf_map.update(
                {
                    f"outside_d_at_cdf_{int(round(q * 100)):02d}_um": float(val)
                    for q, val in zip(cdf_levels, outside_cdf_values, strict=True)
                }
            )
        else:
            cdf_map.update({f"outside_d_at_cdf_{int(round(q * 100)):02d}_um": np.nan for q in cdf_levels})
        if inside_distances_um.size > 0:
            inside_cdf_values = _weighted_quantiles(inside_distances_um, inside_weights, cdf_levels)
            cdf_map.update(
                {
                    f"inside_d_at_cdf_{int(round(q * 100)):02d}_um": float(val)
                    for q, val in zip(cdf_levels, inside_cdf_values, strict=True)
                }
            )
        else:
            cdf_map.update({f"inside_d_at_cdf_{int(round(q * 100)):02d}_um": np.nan for q in cdf_levels})

        side_fraction = {
            "outside_fraction": float(outside_total / total) if total > 0 else np.nan,
            "inside_fraction": float(inside_total / total) if total > 0 else np.nan,
        }

        rows.append(
            DistanceDistributionSummary(
                label=label,
                quantile_signature_um=quantile_map,
                cdf_distance_um=cdf_map,
                side_mass_fraction=side_fraction,
            )
        )
    return rows


def plot_distance_cdf(
    crop_specs: list[tuple[str, str, float, str]],
    output_path: str,
    plots_root: str,
    max_distance_um: float | None = None,
    title: str | None = None,
    origin: str = "boundary",
) -> None:
    """
    Plot pooled experiment transcript-count CDFs using boundary- or peak-based origin.

    Two side-by-side panels show non-negative distances on either side of the
    chosen origin. With `origin="boundary"`, panels correspond to outside and
    inside tissue. With `origin="peak"`, distances are first shifted by the
    detected KDE peak used for peak-aligned histogram/ridgeline plots, and the
    panels correspond to left and right of that peak.

    Args:
        crop_specs: Crop descriptors as `(technology, dataset_id, pixel_size_um, crop_name)`.
        output_path: Destination path for the rendered figure.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.
        max_distance_um: Optional x-limit (0 to `max_distance_um`) for both panels.
        title: Optional figure title.
        origin: Distance origin used for the CDFs. Must be `"boundary"` or `"peak"`.

    Raises:
        ValueError: If `origin` is not one of the supported values.
        ValueError: If `max_distance_um` is provided but not positive.
    """
    data_per_label = read_experiment_distance_data(crop_specs, plots_root)
    if origin not in {"boundary", "peak"}:
        raise ValueError("origin must be 'boundary' or 'peak'.")

    peak_positions_by_label: dict[str, float] = {}
    if origin == "peak":
        curves = _distance_kde_curves(data_per_label, align_peaks=False)
        peak_positions_by_label = {curve.label: curve.alignment_shift_um for curve in curves}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True, gridspec_kw={"wspace": 0.05})
    ax_left, ax_right = axes

    for distribution in data_per_label:
        label = distribution.label
        distances_um, weights = distribution.distances_um, distribution.weights
        distances_um, weights = _sanitize_distance_weights(distances_um, weights, label=label)
        if origin == "boundary":
            left_distances_um, left_weights, right_distances_um, right_weights = _split_signed_distance_sides(
                distances_um,
                weights,
            )
        else:
            if label not in peak_positions_by_label:
                raise ValueError(f"Missing KDE peak position for label {label!r}.")
            left_distances_um, left_weights, right_distances_um, right_weights = _split_peak_aligned_distance_sides(
                distances_um,
                weights,
                peak_position_um=peak_positions_by_label[label],
            )
        color = COLOR_PALETTE.get(distribution.technology, "#7f7f7f")

        for ax, side_distances_um, side_weights in [
            (ax_left, left_distances_um, left_weights),
            (ax_right, right_distances_um, right_weights),
        ]:
            if side_distances_um.size == 0:
                continue
            order = np.argsort(side_distances_um)
            x = side_distances_um[order]
            sorted_weights = side_weights[order]
            cdf = np.cumsum(sorted_weights)
            cdf = cdf / cdf[-1]
            ax.step(x, cdf, where="post", linewidth=1.6, color=color, label=label)

    if max_distance_um is not None:
        if max_distance_um <= 0:
            raise ValueError("max_distance_um must be positive.")
        ax_left.set_xlim(0, max_distance_um)
        ax_right.set_xlim(0, max_distance_um)
    ax_left.set_ylim(0, 1.02)
    ax_left.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_left.grid(axis="y", alpha=0.25)
    ax_right.grid(axis="y", alpha=0.25)
    if origin == "boundary":
        x_label = f"Distance to boundary [{MICROMETER_UNIT}]"
        left_title = "Outside tissue"
        right_title = "Inside tissue"
    else:
        x_label = f"Distance to KDE peak [{MICROMETER_UNIT}]"
        left_title = "Outside tissue" # "Left of peak"
        right_title = "Inside tissue" # "Right of peak"
    ax_left.set_xlabel(x_label)
    ax_right.set_xlabel(x_label)
    ax_left.set_ylabel("Transcript count CDF")
    ax_left.set_title(left_title)
    ax_right.set_title(right_title)
    handles, labels = ax_right.get_legend_handles_labels()
    if len(handles) > 0:
        ax_right.legend(handles, labels, fontsize="small", ncol=1, loc="lower right")
    if title:
        fig.suptitle(title)
    save_or_show_figure(fig, output_path)


def plot_distance_boxplot(
    data: DistanceDistribution,
    output_path: str,
    title: str | None = None,
) -> None:
    distances = np.asarray(data.distances_um, dtype=float)
    weights = np.asarray(data.weights, dtype=float)
    distances, weights = _sanitize_distance_weights(distances, weights, label="distance boxplot")
    if not np.all(np.isclose(weights, np.rint(weights))):
        raise ValueError("Transcript counts must be integer to expand-by-count.")
    repeats = np.rint(weights).astype(int)
    expanded_distances = np.repeat(distances, repeats)
    if expanded_distances.size == 0:
        raise ValueError("No transcript-expanded distances available for boxplot.")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.boxplot(expanded_distances)
    ax.set_ylabel(f"Distance to boundary [{MICROMETER_UNIT}]")
    ax.set_title(title)
    save_or_show_figure(fig, output_path)


def _sanitize_distance_weights(
    distances: np.ndarray,
    weights: np.ndarray,
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    distances = np.asarray(distances, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if distances.shape[0] != weights.shape[0]:
        raise ValueError(f"Distance and weight arrays must align for {label}.")
    valid = np.isfinite(distances) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        raise ValueError(f"No valid positive weights found for {label}.")
    return distances[valid], weights[valid]


def plot_distance_kde_ridgeline(
    crop_specs: list[tuple[str, str, float, str]],
    output_path: str,
    plots_root: str,
    gene_name: str,
    max_distance_um: float | None = None,
    title: str | None = None,
    figsize=(14, 8),
    align_peaks: bool = False,
    overlap: float = 0.6,
    fill_alpha: float = 0.45,
) -> None:
    """
    Plot one pooled signed-distance KDE per experiment as a ridgeline (joy) chart.

    If `align_peaks` is True, each curve is shifted so its KDE peak is placed
    at x=0.
    Curves are stacked with configurable overlap and optional filled areas.

    Args:
        crop_specs: Crop descriptors as `(technology, dataset_id, pixel_size_um, crop_name)`.
        output_path: Destination path for the rendered figure.
        plots_root: Root directory containing per-crop precomputed diffusion outputs.
        gene_name: Name used in the y-axis label.
        max_distance_um: Optional symmetric x-limit in micrometers.
        title: Optional plot title.
        figsize: Matplotlib figure size.
        align_peaks: Whether to shift each curve so its mode is at zero.
        overlap: Vertical overlap fraction between neighboring ridges in `[0, 1)`.
        fill_alpha: Fill alpha for ridge areas.

    Raises:
        ValueError: If no positive KDE values are available or
            `max_distance_um` is provided but not positive.
        AssertionError: If `crop_specs` is empty or `overlap` is outside `[0, 1)`.
    """
    assert len(crop_specs) > 0, "crop_specs must contain at least one entry."
    all_curves = _distance_kde_curves_from_crop_specs(
        crop_specs=crop_specs,
        plots_root=plots_root,
        align_peaks=align_peaks,
    )

    max_height = max(float(np.max(curve.density)) for curve in all_curves)
    if max_height <= 0:
        raise ValueError("No positive KDE values found.")

    step = 1.0 - overlap
    assert 0 <= overlap < 1, "overlap must be in [0, 1)."

    fig, ax = plt.subplots(figsize=figsize)
    num_curves = len(all_curves)
    for idx, curve in enumerate(all_curves):
        baseline = (num_curves - 1 - idx) * step
        assert curve.technology in COLOR_PALETTE
        color = COLOR_PALETTE[curve.technology]
        filled = baseline + curve.density
        ax.plot(curve.x_um, filled, color=color, linewidth=1.4, label=curve.label)
        ax.fill_between(curve.x_um, baseline, filled, color=color, alpha=fill_alpha, linewidth=0)
        ax.scatter(
            [curve.peak_position_um],
            [baseline + curve.peak_density],
            s=18,
            color=color,
            edgecolors="white",
            linewidths=0.6,
            zorder=4,
        )

    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", length=0)
    ax.set_ylabel(f"Normalized {gene_name} expression")
    ax.legend(fontsize="small", bbox_to_anchor=(1.01, 1), loc="upper left")

    ax.set_xlabel(
        f"Peak-aligned distance to boundary [{MICROMETER_UNIT}]"
        if align_peaks
        else f"Distance to boundary [{MICROMETER_UNIT}]"
    )
    ax.set_title(title)
    ax.axvline(0, color="black", linewidth=1, linestyle=":", zorder=3)
    if max_distance_um is not None:
        if max_distance_um <= 0:
            raise ValueError("max_distance_um must be positive.")
        ax.set_xlim(-max_distance_um, max_distance_um)

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
