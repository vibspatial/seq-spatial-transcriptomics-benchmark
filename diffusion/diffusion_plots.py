from __future__ import annotations

import os
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
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
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
    ), f"Expected xarray.DataTree for sdata.labels[{labels_layer!r}], got {type(sdata.labels[labels_layer])}"

    figsize_factor = dpi
    figsize = (
        sdata.labels[labels_layer].shape[1] / figsize_factor,
        sdata.labels[labels_layer].shape[0] / figsize_factor,
    )
    print(f"figsize: {figsize}")

    return dpi, figsize_factor, figsize


def save_or_show_figure(fig: plt.Figure, output_folder: str | None, filename: str) -> None:
    if output_folder is not None:
        save_path = os.path.join(output_folder, filename)
        fig.savefig(
            save_path,
            bbox_inches="tight",
            pad_inches=max(fig.get_size_inches()) / 100,
            facecolor="white",
            transparent=False,
        )
        print(f"{save_path} saved")
    else:
        plt.show()
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


def plot_bins_colored_by_feature(
    sdata,
    HE_image: str | None,
    technology: str,
    dataset_id: str,
    level_info,
    feature: str,
    output_folder: str | None,
    filename: str | None = None,
    add_scale_bar: bool | None = None,
    pixel_size_um: float | None = None,
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
        save_or_show_figure(fig, output_folder, filename)

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
        save_or_show_figure(fig, output_folder, filename)


def plot_tissue_outline_on_HE_image(
    sdata,
    HE_image: str,
    tissue_mask: str,
    technology: str,
    dataset_id: str,
    output_folder: str | None,
    filename: str | None = None,
    fill_alpha=0.3,
    outline_width=5,
    bounds: tuple[float, float, float, float] | None = None,
    add_scale_bar: bool = False,
    pixel_size_um: float | None = None,
):
    # Plot H&E image with tissue outline on top.
    # Optionally restrict the plot to bounds = (xmin, ymin, xmax, ymax).

    best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

    title = f"{technology} {dataset_id}: H&E with tissue outline"
    if bounds is not None:
        xmin, ymin, xmax, ymax = bounds
        title = f"{title}\n(crop: {xmin}, {ymin}, {xmax}, {ymax})"

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

    save_or_show_figure(fig, output_folder, filename)


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
    figsize: tuple[float, float] = (6, 6),
    output_folder: str | None = None,
    filename: str | None = None,
):
    fig, ax = plt.subplots(figsize=figsize)

    gpd.GeoSeries([boundary]).plot(ax=ax, figsize=figsize)
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_bin_gene_counts_and_contour(
    gene_bin_data: "GeneBinData",
    contour,
    title: str,
    figsize=(10, 10),
    line_width=1,
    marker_size=20,
    dpi: int | None = None,
    output_folder: str | None = None,
    filename: str | None = None,
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

    save_or_show_figure(fig, output_folder, filename)


def plot_bins_and_closest_points_on_contour(
    data: "DiffusionBinData",
    contour,
    title: str,
    line_width=1,
    output_folder: str | None = None,
    filename: str | None = None,
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

    save_or_show_figure(fig, output_folder, filename)


def plot_distances_linearized(
    data: "DiffusionBinData",
    title: str,
    pixel_size_um: float | None = None,
    output_folder: str | None = None,
    filename: str | None = None,
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

    save_or_show_figure(fig, output_folder, filename)


def plot_distance_histogram(
    data: "DiffusionBinData",
    output_folder: str | None,
    filename: str,
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
    ax.set_xlabel(f"Distance to boundary [{unit_label}]")
    ax.set_ylabel("Count")
    ax.set_title(title or f"Distance to boundary histogram [{unit_label}]")
    save_or_show_figure(fig, output_folder, filename)


def plot_distance_boxplot(
    data: "DiffusionBinData",
    output_folder: str | None,
    filename: str,
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
    save_or_show_figure(fig, output_folder, filename)
