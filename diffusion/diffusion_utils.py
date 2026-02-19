import os
import sys
from dataclasses import dataclass
sys.path.append(os.path.abspath('..'))
from diffusion import Experiment
from diffusion_plots import (
    DiffusionPlotSettings,
    plot_bin_gene_counts_and_contour,
    plot_bins_and_closest_points_on_contour,
    plot_boundary_geoseries,
    plot_distance_boxplot,
    plot_distance_histogram,
    plot_distances_linearized,
    plot_tissue_outline_on_HE_image,
)

import shapely
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.geometry.base import BaseGeometry
import spatialdata as sd
import numpy as np


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
) -> "DistanceData":
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
    assert len(points) == len(distances)
    assert len(points) == len(along_distances)
    assert len(points) == closest_xy.shape[0]

    print("...done")
    return DistanceData(
        closest_xy=closest_xy,
        closest_distances=distances,
        closest_along_boundary=along_distances,
    )


def get_experiment_by_id(experiments: list[Experiment], dataset_id: str) -> Experiment | None:
    for experiment in experiments:
        if experiment.dataset_id == dataset_id:
            return experiment
    return None

@dataclass
class DistanceData:
    closest_xy: np.ndarray
    closest_distances: np.ndarray
    closest_along_boundary: np.ndarray


@dataclass
class GeneBinData:
    bin_centers: np.ndarray
    transcript_counts: np.ndarray


@dataclass
class DiffusionBinData:
    gene_bins: GeneBinData
    distance_data: DistanceData

    def filter(self, mask: np.ndarray) -> "DiffusionBinData":
        return DiffusionBinData(
            gene_bins=GeneBinData(
                bin_centers=self.gene_bins.bin_centers[mask],
                transcript_counts=self.gene_bins.transcript_counts[mask],
            ),
            distance_data=DistanceData(
                closest_xy=self.distance_data.closest_xy[mask],
                closest_distances=self.distance_data.closest_distances[mask],
                closest_along_boundary=self.distance_data.closest_along_boundary[mask],
            ),
        )


def get_gene_bin_data(sdata, table_name: str, gene_name: str) -> GeneBinData:

    # Get the feature column number for the gene we want to check diffusion for
    gene_index: int = list(sdata.tables[table_name].var.index).index(gene_name)
    assert sdata.tables[table_name].var.iloc[gene_index].name == gene_name

    # Get the number of transcripts for the given gene, for each bin.
    transcript_counts: np.ndarray = sdata.tables[table_name].X[:, gene_index].toarray().flatten()

    # Get the (x,y) coordinates of the centers of each bin
    bin_centers: np.ndarray = sdata.tables[table_name].obsm["spatial"]

    # We must have a transcript count for exactly each bin
    assert transcript_counts.shape[0] == bin_centers.shape[0]

    return GeneBinData(
        bin_centers=bin_centers,
        transcript_counts=transcript_counts,
    )


def analyze_diffusion(
    experiment: Experiment,
    level: str,
    gene_name: str,
    max_distance_px: float = 400,
    use_signed_distances: bool = False,
    base_folder: str | None = None,
    plot_settings: DiffusionPlotSettings | None = None,
) -> None:
    if plot_settings is None:
        plot_settings = DiffusionPlotSettings()

    print(f"Analyzing diffusion for {experiment.technology} {experiment.dataset_id} at level {level}")

    sdata = sd.read_zarr(experiment.sdata_path)

    for crop_idx, crop_rect in enumerate(experiment.crop_rects):
        print(f"Using crop rect: {crop_rect}")

        # Create a cropped spatial data
        crop_rect = np.array(crop_rect)
        sdata_crop = sdata.query.bounding_box(
            min_coordinate=crop_rect[:2],
            max_coordinate=crop_rect[2:],
            axes=("x", "y"),
            target_coordinate_system=experiment.dataset_id,
        )

        if base_folder:
            plot_folder = os.path.join(
                base_folder,
                experiment.technology,
                experiment.dataset_id,
                "analysis",
                f"crop{crop_idx + 1}",
            )
            os.makedirs(plot_folder, exist_ok=True)
        else:
            plot_folder = None

        # Plot close-up of cropped area
        if plot_settings.tissue_outline and experiment.he_image:
            plot_tissue_outline_on_HE_image(
                sdata_crop,
                experiment.he_image,
                experiment.tissue_mask,
                experiment.technology,
                experiment.dataset_id,
                plot_folder,
                "tissue.png",
                fill_alpha=0.1,
                outline_width=1,
                bounds=crop_rect,
                add_scale_bar=experiment.pixel_size_um is not None,
                pixel_size_um=experiment.pixel_size_um,
            )

        # Calculate tissue boundary against which we will calculate diffusion distances
        boundary = longest_cropped_shape_contour(sdata_crop, experiment.tissue_mask, crop_rect)
        if plot_settings.boundary:
            # TODO? perhaps add title to plot
            plot_boundary_geoseries(
                boundary,
                output_folder=plot_folder,
                filename="boundary.png",
            )

        gene_bin_data = get_gene_bin_data(sdata_crop, level, gene_name)

        [xmin, ymin, xmax, ymax] = crop_rect
        title = (
            f"{experiment.dataset_id}_{level}: {gene_name} transcript counts per bin\n"
            f"(crop: {xmin}, {ymin}, {xmax}, {ymax})"
        )
        if plot_settings.bin_gene_counts:
            plot_bin_gene_counts_and_contour(
                gene_bin_data,
                boundary,
                title,
                dpi=600,
                output_folder=plot_folder,
                filename=f"counts_{gene_name}.png",
                add_scale_bar=experiment.pixel_size_um is not None,
                pixel_size_um=experiment.pixel_size_um,
            )

        # Calculate closest points on contour and distances
        distance_data = closest_points_and_distances(
            gene_bin_data.bin_centers,
            boundary,
            use_signed_distances=use_signed_distances,
        )
        filtered = DiffusionBinData(
            distance_data=distance_data,
            gene_bins=gene_bin_data,
        )

        # Keep only the bins that are within a certain maximum distance to the contour
        filtered = filtered.filter(np.abs(filtered.distance_data.closest_distances) <= max_distance_px)

        # Keep only the bins with at least one transcript of the gene of interest
        filtered = filtered.filter(filtered.gene_bins.transcript_counts > 0)

        # IMPORTANT FIXME: take number of transcripts in bin into account when analyzing diffusion distances!

        title = (
            f"{level}: bins and closest point on boundary\n"
            f"(crop: {xmin}, {ymin}, {xmax}, {ymax})"
        )
        # TODO: add gene name to title?
        if plot_settings.bins_and_closest_points:
            plot_bins_and_closest_points_on_contour(
                filtered,
                boundary,
                title,
                line_width=0.5,
                dpi=600,
                output_folder=plot_folder,
                filename="bins_closest_point_on_boundary.png",
            )

        title = (
            f"{level}: bin distance to boundary vs position along boundary\n"
            f"(crop: {xmin}, {ymin}, {xmax}, {ymax})"
        )
        if plot_settings.distances_linearized:
            plot_distances_linearized(
                filtered,
                title,
                pixel_size_um=experiment.pixel_size_um,
                dpi=600,
                output_folder=plot_folder,
                filename="distances_linearized.png",
            )

        if plot_settings.histogram:
            plot_distance_histogram(
                filtered,
                plot_folder,
                "histogram_px.png",
                unit="px",
            )
            if experiment.pixel_size_um is not None and experiment.pixel_size_um > 0:
                plot_distance_histogram(
                    filtered,
                    plot_folder,
                    "histogram_um.png",
                    unit="um",
                    pixel_size_um=experiment.pixel_size_um,
                )

        if plot_settings.boxplot:
            plot_distance_boxplot(
                filtered,
                plot_folder,
                "boxplot_px.png",
                unit="px",
            )
            if experiment.pixel_size_um is not None and experiment.pixel_size_um > 0:
                plot_distance_boxplot(
                    filtered,
                    plot_folder,
                    "boxplot_um.png",
                    unit="um",
                    pixel_size_um=experiment.pixel_size_um,
                )
