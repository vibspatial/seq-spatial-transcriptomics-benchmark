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
                                 output_folder: str | None,
                                 filename: str | None = None):
    
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
        save_or_show_figure(fig, output_folder, filename)

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
        save_or_show_figure(fig, output_folder, filename)



def plot_tissue_outline_on_HE_image(sdata,
                                    HE_image: str,
                                    tissue_mask: str,
                                    technology: str,
                                    dataset_id: str,
                                    output_folder: str | None,
                                    filename: str | None = None,
                                    fill_alpha=0.3,
                                    outline_width=5,
                                    bounds: tuple[float, float, float, float] | None = None):
    # Plot H&E image with tissue outline on top.
    # Optionally restrict the plot to bounds = (xmin, ymin, xmax, ymax).

    best_scale, dpi, _, figsize = get_image_figure_parameters(sdata, HE_image)

    title = f"{technology} {dataset_id}: H&E with tissue outline"
    if bounds is not None:
        xmin, ymin, xmax, ymax = bounds
        title = f"{title}\n(crop: {xmin}, {ymin}, {xmax}, {ymax})"

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

    save_or_show_figure(fig, output_folder, filename)


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


def plot_boundary_geoseries(boundary: BaseGeometry,
                            figsize: tuple[float, float] = (6, 6),
                            output_folder: str | None = None,
                            filename: str | None = None):
    fig, ax = plt.subplots(figsize=figsize)

    gpd.GeoSeries([boundary]).plot(ax=ax, figsize=figsize)
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def signed_distance_to_linestring(
    points: np.ndarray,
    line: shapely.LineString,
    nearest_on_line: BaseGeometry | None = None,
    closest_xy: np.ndarray | None = None,
    along_distances: np.ndarray | None = None,
) -> np.ndarray:
    """
    Signed distance to a LineString.
    - Use closest point on interior segments.
    - If the closest segment is the first and the projection is before the start,
      extend the first segment to infinity.
    - If the closest segment is the last and the projection is after the end,
      extend the last segment to infinity.
    """
    points_xy = np.asarray(points, dtype=float)
    if points_xy.ndim != 2 or points_xy.shape[1] != 2:
        raise ValueError("points must be an (N, 2) array.")

    if nearest_on_line is None:
        points_geom = shapely.points(points_xy)
        _, nearest_on_line = shapely.ops.nearest_points(points_geom, line)
    if closest_xy is None:
        closest_xy = shapely.get_coordinates(nearest_on_line)
    if along_distances is None:
        along_distances = shapely.line_locate_point(line, nearest_on_line, normalized=False)

    coords = np.asarray(line.coords, dtype=float)
    if coords.shape[0] < 2:
        raise ValueError("line must have at least two coordinates.")

    seg_starts = coords[:-1]
    seg_ends = coords[1:]
    seg_vec = seg_ends - seg_starts
    seg_len = np.sqrt(np.einsum("ij,ij->i", seg_vec, seg_vec))
    valid = seg_len > 0
    if not np.any(valid):
        raise ValueError("line has no non-degenerate segments.")

    seg_starts = seg_starts[valid]
    seg_vec = seg_vec[valid]
    seg_len = seg_len[valid]
    seg_len2 = seg_len * seg_len

    cum_len = np.concatenate(([0.0], np.cumsum(seg_len)))
    seg_idx = np.searchsorted(cum_len, along_distances, side="right") - 1
    seg_idx = np.clip(seg_idx, 0, seg_vec.shape[0] - 1)

    seg_start_sel = seg_starts[seg_idx]
    seg_vec_sel = seg_vec[seg_idx]
    seg_len_sel = seg_len[seg_idx]
    seg_len2_sel = seg_len2[seg_idx]

    rel = points_xy - seg_start_sel
    t_sel = np.einsum("ij,ij->i", rel, seg_vec_sel) / seg_len2_sel

    use_first_ext = (seg_idx == 0) & (t_sel < 0.0)
    use_last_ext = (seg_idx == (seg_vec.shape[0] - 1)) & (t_sel > 1.0)
    use_ext = use_first_ext | use_last_ext

    seg_dir = seg_vec_sel / seg_len_sel[:, None]

    vec_local = points_xy - closest_xy
    vec_ext = points_xy - seg_start_sel
    vec_for_cross = np.where(use_ext[:, None], vec_ext, vec_local)
    cross = seg_dir[:, 0] * vec_for_cross[:, 1] - seg_dir[:, 1] * vec_for_cross[:, 0]

    dist_local = np.sqrt(np.sum(vec_local * vec_local, axis=1))
    dist_ext = np.abs(cross)
    dist = np.where(use_ext, dist_ext, dist_local)

    return dist * np.sign(cross)


def closest_points_and_distances(
    points: np.ndarray,
    line: shapely.LineString,
    use_signed_distances: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    print("Calculating distances to and along boundary contour...")
    points_xy = np.asarray(points, dtype=float)
    points_geom = shapely.points(points_xy)

    _, nearest_on_line = shapely.ops.nearest_points(points_geom, line)

    closest_xy: np.ndarray = shapely.get_coordinates(nearest_on_line)
    if use_signed_distances:
        along_distances_abs = shapely.line_locate_point(line, nearest_on_line, normalized=False)
        distances = signed_distance_to_linestring(
            points_xy,
            line,
            nearest_on_line=nearest_on_line,
            closest_xy=closest_xy,
            along_distances=along_distances_abs,
        )
    else:
        distances = shapely.distance(points_geom, nearest_on_line)
    along_distances: np.ndarray = shapely.line_locate_point(line, nearest_on_line, normalized=True)

    assert closest_xy.shape[1] == 2
    assert len(distances.shape) == 1
    assert len(along_distances.shape) == 1

    print("...done")
    return closest_xy, distances, along_distances


def plot_bin_gene_counts_and_contour(
    bin_centers,
    transcript_counts,
    contour,
    title: str,
    figsize=(10, 10),
    line_width=1,
    marker_size=20,
    output_folder: str | None = None,
    filename: str | None = None
) -> None:
    # Plot the bin centers, colored by the number of transcripts detected in that bin,
    # together with a given contour over it.
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title)
    sc = ax.scatter(bin_centers[:,0], bin_centers[:, 1], c=transcript_counts, s=marker_size)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.4)
    cbar.set_label("Counts")
    gpd.GeoSeries([contour]).plot(
        ax=ax,
        color="red",
        linewidth=line_width,
    )
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_bins_and_closest_points_on_contour(points_xy, closest_xy, closest_distances, contour, title: str, output_folder: str | None = None, filename: str | None = None):

    # Line segments connecting each bin to its closest point on the contour
    segments = np.stack([points_xy, closest_xy], axis=1)  # shape (number of bins, 2, 2)
    lc = LineCollection(
        segments,
        # cmap='viridis',
        array=closest_distances,  # color line by distance
        linewidths=1
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
    ax.set_title(title)
    ax.invert_yaxis()

    save_or_show_figure(fig, output_folder, filename)


def plot_filtered_bins_scatter(
    closest_along_boundary: np.ndarray,
    closest_distances: np.ndarray,
    transcript_counts: np.ndarray,
    title: str,
    output_folder: str | None = None,
    filename: str | None = None,
    figsize: tuple[float, float] = (10, 2.5),
) -> None:
    # Scatter plot of filtered bins: along-boundary position vs distance to boundary.
    sizes = np.asarray(transcript_counts, dtype=float)
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(closest_along_boundary, closest_distances, s=sizes, alpha=0.7)
    ax.set_xlabel("Closest along boundary (normalized)")
    ax.set_ylabel("Closest distance")
    ax.set_title(title)

    save_or_show_figure(fig, output_folder, filename)


def get_experiment_by_id(experiments: list[Experiment], dataset_id: str) -> Experiment | None:
    for experiment in experiments:
        if experiment.dataset_id == dataset_id:
            return experiment
    return None

