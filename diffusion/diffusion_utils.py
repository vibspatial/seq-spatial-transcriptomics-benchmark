import os
import sys
sys.path.append(os.path.abspath('..'))
from utils.in_tissue import add_in_tissue, assign_region
from utils.plotting import best_datatree_scale, rel_fontsize
from diffusion import LevelInfo, Experiment

import shapely
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.geometry.base import BaseGeometry
import spatialdata_plot
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
import numpy as np
import geopandas as gpd
import xarray as xr


# from shapely.geometry import box, Polygon

# def crop_polygon_to_bbox(polygon: Polygon, crop_rect):
#     """
#     Crop a Shapely Polygon to a rectangular bounding box.

#     Parameters
#     ----------
#     polygon : shapely.geometry.Polygon
#         Input polygon to be cropped.
#     xmin, ymin, xmax, ymax : float
#         Bounding box coordinates.

#     Returns
#     -------
#     shapely.geometry.base.BaseGeometry
#         Cropped geometry (Polygon, MultiPolygon, or empty geometry).
#     """
#     xmin, ymin, xmax, ymax = crop_rect
#     bbox = box(xmin, ymin, xmax, ymax)
#     return polygon.intersection(bbox)

def clip_polygon_to_bbox_lines(polygon, crop_rect):
    """
    Clip a polygon to a rectangle and return only the generated polylines,
    excluding any parts that lie on the rectangle edges.

    Returns a LineString or MultiLineString (or empty geometry).
    """
    xmin, ymin, xmax, ymax = crop_rect
    bbox = box(xmin, ymin, xmax, ymax)

    # Clip polygon
    clipped = polygon.intersection(bbox)

    if clipped.is_empty:
        return clipped

    # Extract boundary of clipped geometry
    clipped_boundary = clipped.boundary

    # Remove rectangle edges
    result = clipped_boundary.difference(bbox.boundary)

    return result


def _cropped_shape_boundary(
    geom: BaseGeometry,
    crop_rect: tuple[float, float, float, float],
) -> list[shapely.LineString]:
    def _iter_lines_local(geometry: BaseGeometry):
        if geometry.is_empty:
            return
        for part in shapely.get_parts(geometry):
            if part.geom_type == "LineString":
                yield part

    bbox = box(*crop_rect)
    if geom is None or geom.is_empty:
        return []
    if not geom.intersects(bbox):
        return []

    clipped = geom.intersection(bbox)
    if clipped.is_empty:
        return []

    # Remove bbox edges so only interior polygon boundaries remain.
    parts = list(shapely.get_parts(clipped))
    poly_parts = [p for p in parts if p.geom_type in {"Polygon", "MultiPolygon"}]
    if poly_parts:
        boundary_geom = unary_union(poly_parts).boundary
    else:
        boundary_geom = clipped

    if boundary_geom is None or boundary_geom.is_empty:
        return []
    boundary_geom = boundary_geom.difference(bbox.boundary)
    if boundary_geom is None or boundary_geom.is_empty:
        return []

    return list(_iter_lines_local(boundary_geom))


def longest_cropped_shape_contour(
    sdata,
    shapes_layer: str,
    crop_rect: tuple[float, float, float, float],
) -> shapely.LineString | None:
    """
    Return the longest clipped contour (LineString) from a shapes layer.

    If no shape overlaps the crop rectangle, return None.
    """
    gdf = sdata.shapes[shapes_layer]
    if gdf.empty:
        return None

    lines = [line for geom in gdf.geometry for line in _cropped_shape_boundary(geom, crop_rect)]

    if not lines:
        return None
    return max(lines, key=lambda line: line.length)


def get_image_figure_parameters(sdata, HE_image, dpi=200, max_dim_for_plots=10000):
    assert isinstance(sdata[HE_image], xr.DataTree), f"Expected xarray.DataTree for sdata[{HE_image!r}], got {type(sdata[HE_image])}"

    best_scale = best_datatree_scale(sdata[HE_image], max_dim_for_plots)
    print(f'best scale: {best_scale}')

    figsize_factor = dpi
    figsize = (sdata.images[HE_image][best_scale]['image'].shape[2]/figsize_factor,
               sdata.images[HE_image][best_scale]['image'].shape[1]/figsize_factor)
    print(f'figsize: {figsize}')

    return best_scale, dpi, figsize_factor, figsize


def get_labels_figure_parameters(sdata, labels_layer, dpi=200, max_dim_for_plots=10000):
    assert isinstance(sdata.labels[labels_layer], xr.DataArray), f"Expected xarray.DataTree for sdata.labels[{labels_layer!r}], got {type(sdata.labels[labels_layer])}"

    figsize_factor = dpi
    figsize = (sdata.labels[labels_layer].shape[1]/figsize_factor,
               sdata.labels[labels_layer].shape[0]/figsize_factor)
    print(f'figsize: {figsize}')

    return dpi, figsize_factor, figsize


def save_or_show_figure(fig: plt.Figure, output_folder: str | None, filename: str) -> None:
    if output_folder is not None:
        save_path = os.path.join(output_folder, filename)
        fig.savefig(save_path, bbox_inches="tight", pad_inches=max(fig.get_size_inches())/100, facecolor="white", transparent=False)
        print(f"{save_path} saved")
    else:
        plt.show()
    plt.close(fig)


def plot_bins_colored_by_feature(sdata,
                                 HE_image: str | None,
                                 technology: str,
                                 dataset_id: str,
                                 level_info: LevelInfo,
                                 feature: str,
                                 output_folder: str | None):
    
    # color can be for example 'total_counts' or the name of a gene.

    title = f"{feature} ({technology} {dataset_id} {level_info.name})"

    if HE_image:
        best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        sdata.pl.render_images(
            element=HE_image, scale=best_scale
        ).pl.render_labels(
            element=level_info.labels_layer, fill_alpha=0.5, contour_px=5, outline_alpha=0, 
            color=feature, 
            table_name=level_info.table_name,
        ).pl.show(
            coordinate_systems=dataset_id,
            title=title,
            figsize=figsize,
            colorbar=False,
            dpi=dpi,
            ax=ax
        )
        save_or_show_figure(fig, output_folder, f'{technology}_{dataset_id}_{level_info.name}_{feature}.png')

    else:
        dpi, _, figsize = get_labels_figure_parameters(sdata, level_info.labels_layer)

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        sdata.pl.render_labels(
            element=level_info.labels_layer, fill_alpha=0.5, contour_px=5, outline_alpha=0, 
            color=feature, 
            table_name=level_info.table_name,
        ).pl.show(
            coordinate_systems=dataset_id,
            title=title,
            figsize=figsize,
            colorbar=False,
            dpi=dpi,
            ax=ax
        )
        save_or_show_figure(fig, output_folder, f'{technology}_{dataset_id}_{level_info.name}_{feature}.png')



def plot_tissue_outline_on_HE_image(sdata,
                                    HE_image: str,
                                    tissue_mask: str,
                                    technology: str,
                                    dataset_id: str,
                                    output_folder: str | None,
                                    fill_alpha=0.3,
                                    outline_width=5,
                                    bounds: tuple[float, float, float, float] | None = None):
    # Plot H&E image with tissue outline on top.
    # Optionally restrict the plot to bounds = (xmin, ymin, xmax, ymax).

    best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

    title = f"{technology} {dataset_id}: H&E with tissue outline"

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    sdata.pl.render_images(
        element=HE_image, scale=best_scale
    ).pl.render_shapes(
        element=tissue_mask, color='#1f77b4', fill_alpha=fill_alpha, outline_width=outline_width, outline_alpha=1, outline_color='red'
    ).pl.show(
        coordinate_systems=dataset_id,
        title=title,
        figsize=figsize,
        colorbar=False,
        dpi=dpi,
        ax=ax
    )
    if bounds is not None:
        # Restrict plot to the requested bounds and show the same bounds on axes.
        xmin, ymin, xmax, ymax = bounds
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)

    save_or_show_figure(fig, output_folder, f'{technology}_{dataset_id}_tissue.png')


def plot_bin_outlines(crop_rect, bins_center_coordinates, bin_outlines_gdf, title: str):
    # crop_rect = [xmin, ymin, xmax, ymax] 
    # bins_center_coordinates = numpy array of nx2
    # bin_outlines_gdf: geopandas dataframe with outlines of the bins
    # print(type(crop_rect))
    # print(type(bins_center_coordinates))
    # print(type(bin_outlines_gdf))
    rect = Rectangle(
        (crop_rect[0], crop_rect[1]),          # lower-left corner
        crop_rect[2] - crop_rect[0],           # width
        crop_rect[3] - crop_rect[1],           # height
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


def plot_boundary_geoseries(boundary: BaseGeometry, figsize: tuple[float, float] = (6, 6), invert_y: bool = True):
    ax = gpd.GeoSeries([boundary]).plot(figsize=figsize)
    if invert_y:
        ax.set_ylim(ax.get_ylim()[::-1])
    plt.show()


def closest_points_and_distances(points: np.ndarray, line: shapely.LineString) -> tuple[np.ndarray, np.ndarray]:
    points = shapely.points(points)

    _, nearest_on_line = shapely.ops.nearest_points(points, line)

    closest_xy: np.ndarray = shapely.get_coordinates(nearest_on_line)
    distances: np.ndarray = shapely.distance(points, nearest_on_line)

    assert closest_xy.shape[1] == 2
    assert len(distances.shape) == 1

    return closest_xy, distances


def plot_bin_gene_counts_and_contour(bin_centers, transcript_counts, contour, title:str) -> None:
    # Plot the bin centers, colored by the number of transcripts detected in that bin,
    # together with a given contour over it.
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title(title)
    sc = ax.scatter(bin_centers[:,0], bin_centers[:, 1], c=transcript_counts, s=20)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.4)
    cbar.set_label("Counts")
    gpd.GeoSeries([contour]).plot(
        ax=ax,
        color="red",
        linewidth=1
    )


def plot_bins_and_closest_points_on_contour(points_xy, closest_xy, closest_distances, contour, dataset_and_level):

    # Line segments connecting each bin to its closest point on the contour
    segments = np.stack([points_xy, closest_xy], axis=1)  # shape (number of bins, 2, 2)
    lc = LineCollection(
        segments,
        # cmap='viridis',
        array=closest_distances,  # color line by distance
        linewidths=1.5 + 2*closest_distances/closest_distances.max()  # scale line width by distance
    )

    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the contour
    line_xy = np.array(contour.coords)
    ax.plot(line_xy[:, 0], line_xy[:, 1], '-', c="red", lw=1, label='Contour')

    # Plot the bin centers
    ax.scatter(points_xy[:, 0], points_xy[:, 1], c='blue', s=1, label='Bin centers')

    # Draw connecting line segments
    ax.add_collection(lc)

    ax.set_aspect('equal')
    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f'{dataset_and_level}: bins and closest point on contour')

    plt.show()


def get_experiment_by_id(experiments: list[Experiment], dataset_id: str) -> Experiment | None:
    for experiment in experiments:
        if experiment.dataset_id == dataset_id:
            return experiment
    return None
