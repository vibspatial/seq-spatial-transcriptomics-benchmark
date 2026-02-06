from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import spatialdata_plot  # noqa: F401
import xarray as xr
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from shapely.geometry.base import BaseGeometry

import sys
sys.path.append(os.path.abspath('..'))
from utils.plotting import best_datatree_scale

if TYPE_CHECKING:
    from diffusion_utils import DiffusionBinData, GeneBinData


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


def plot_bins_colored_by_feature(
    sdata,
    HE_image: str | None,
    technology: str,
    dataset_id: str,
    level_info,
    feature: str,
    output_folder: str | None,
    filename: str | None = None,
):
    # color can be for example 'total_counts' or the name of a gene.

    title = f"{feature} ({technology} {dataset_id} {level_info.name})"

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
    if bounds is not None:
        # Restrict plot to the requested bounds and show the same bounds on axes.
        xmin, ymin, xmax, ymax = bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)

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
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_bin_gene_counts_and_contour(
    gene_bin_data: "GeneBinData",
    contour,
    title: str,
    figsize=(10, 10),
    line_width=1,
    marker_size=20,
    output_folder: str | None = None,
    filename: str | None = None,
) -> None:
    # Plot the bin centers, colored by the number of transcripts detected in that bin,
    # together with a given contour over it.
    fig, ax = plt.subplots(figsize=figsize)
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
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_bins_and_closest_points_on_contour(
    data: "DiffusionBinData",
    contour,
    title: str,
    line_width=1,
    output_folder: str | None = None,
    filename: str | None = None,
):
    fig, ax = plt.subplots(figsize=(6, 6))

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
    ax.legend()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(title)
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_distances_linearized(
    data: "DiffusionBinData",
    title: str,
    output_folder: str | None = None,
    filename: str | None = None,
    figsize: tuple[float, float] = (10, 2.5),
) -> None:
    # Scatter plot of bins: along-boundary position vs distance to boundary.
    sizes = np.asarray(data.gene_bins.transcript_counts, dtype=float)
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        data.distance_data.closest_along_boundary,
        data.distance_data.closest_distances,
        s=sizes,
        alpha=0.7,
    )
    ax.set_xlabel("Normalized position along boundary")
    ax.set_ylabel("Distance to boundary")
    ax.set_title(title)

    save_or_show_figure(fig, output_folder, filename)


def plot_distance_histogram(
    data: "DiffusionBinData",
    output_folder: str | None,
    filename: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist(data.distance_data.closest_distances)  # TODO: add title
    save_or_show_figure(fig, output_folder, filename)


def plot_distance_boxplot(
    data: "DiffusionBinData",
    output_folder: str | None,
    filename: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.boxplot(data.distance_data.closest_distances)  # TODO: add title
    save_or_show_figure(fig, output_folder, filename)
