import os
import sys
from dataclasses import dataclass
from time import perf_counter
sys.path.append(os.path.abspath('..'))
from diffusion import Experiment
from diffusion_plots import (
    DiffusionPlotSettings,
    read_experiment_distance_data,
    plot_bin_gene_counts_and_contour,
    plot_bins_and_closest_points_on_contour,
    plot_boundary_geoseries,
    plot_distance_boxplot,
    plot_distance_histogram,
    plot_distances_linearized,
    plot_full_tissue_bin_gene_counts_with_highlighted_section,
    plot_tissue_outline_on_HE_image,
)

import shapely
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.geometry.polygon import orient
from shapely.ops import substring
from shapely.geometry.base import BaseGeometry
import spatialdata as sd
import numpy as np

# Note about signed distances and contour orientation:
# - points outside the closed tissue annotation contour have a negative distance to it;
# - points inside have a positive distance to it;
# - in a coordinate system with the origin at the top-left corner, x-axis pointing
#   to the right, and y-axis pointing down, the tissue annotation contour is oriented counter-clockwise.


def _signed_ring_area(coords: np.ndarray | list[tuple[float, float]]) -> float:
    """Compute signed area for a single linear ring."""
    ring = np.asarray(coords, dtype=float)
    if ring.ndim != 2 or ring.shape[0] < 3:
        return 0.0
    x = ring[:, 0]
    y = ring[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def signed_area(polygon: Polygon) -> float:
    """
    Return signed polygon area using shapely.orientation sign convention.

    Positive means the oriented exterior ring is CCW, negative means clockwise.
    This matches shapely.geometry.polygon.orient(..., sign=1.0).
    """
    if not isinstance(polygon, Polygon):
        raise TypeError(f"Expected a shapely Polygon, got {type(polygon)!r}.")
    if polygon.is_empty:
        return 0.0
    oriented = orient(polygon, sign=1.0)
    area = _signed_ring_area(oriented.exterior.coords)
    for interior in oriented.interiors:
        area += _signed_ring_area(interior.coords)
    return area


def _iter_lines(geometry: BaseGeometry):
    if geometry is None or geometry.is_empty:
        return
    for part in shapely.get_parts(geometry):
        if part.geom_type == "LineString":
            yield part
        elif part.geom_type == "LinearRing":
            yield LineString(part.coords)


def _shape_boundary_lines(geom: BaseGeometry) -> list[LineString]:
    if geom is None or geom.is_empty:
        return []
    boundary = geom.boundary
    if boundary is None or boundary.is_empty:
        return []
    return list(_iter_lines(boundary))


def longest_shape_contour(
    sdata,
    shapes_layer: str,
) -> LineString | None:
    """Return the longest contour (LineString) from a shapes layer."""
    gdf = sdata.shapes[shapes_layer]
    if gdf.empty:
        return None
    lines = [line for geom in gdf.geometry for line in _shape_boundary_lines(geom)]
    if not lines:
        return None
    return max(lines, key=lambda line: line.length)


def full_shape_boundary_geometry(
    sdata,
    shapes_layer: str,
) -> BaseGeometry | None:
    """Return the complete boundary geometry from a shapes layer."""
    gdf = sdata.shapes[shapes_layer]
    if gdf.empty:
        return None
    lines = [line for geom in gdf.geometry for line in _shape_boundary_lines(geom)]
    if not lines:
        return None
    if len(lines) == 1:
        return lines[0]
    return MultiLineString([line.coords for line in lines])


def _is_closed_linestring(line: LineString, atol: float = 1e-8) -> bool:
    coords = np.asarray(line.coords, dtype=float)
    return coords.shape[0] >= 2 and np.allclose(coords[0], coords[-1], atol=atol)


def _project_point_to_line(point_xy: tuple[float, float], line: LineString) -> tuple[float, np.ndarray]:
    p = Point(float(point_xy[0]), float(point_xy[1]))
    d = float(shapely.line_locate_point(line, p, normalized=False))
    snapped_geom = shapely.line_interpolate_point(line, d, normalized=False)
    snapped_xy = np.asarray(shapely.get_coordinates(snapped_geom), dtype=float).reshape(-1, 2)[0]
    return d, snapped_xy


def _concatenate_lines(line_a: LineString, line_b: LineString) -> LineString:
    coords_a = np.asarray(line_a.coords, dtype=float)
    coords_b = np.asarray(line_b.coords, dtype=float)
    if coords_a.shape[0] == 0:
        return LineString(coords_b)
    if coords_b.shape[0] == 0:
        return LineString(coords_a)
    if np.allclose(coords_a[-1], coords_b[0], atol=1e-8):
        coords = np.vstack([coords_a, coords_b[1:]])
    else:
        coords = np.vstack([coords_a, coords_b])
    return LineString(coords)


def oriented_boundary_section(
    contour: LineString,
    start_distance: float,
    end_distance: float,
) -> LineString:
    """
    Extract boundary arc from start to end following the existing contour orientation.

    For closed contours and end < start, the arc wraps around contour end.
    """
    length = float(contour.length)
    if length <= 0:
        raise ValueError("Contour length must be positive.")
    if _is_closed_linestring(contour):
        start_distance = float(start_distance % length)
        end_distance = float(end_distance % length)
        if np.isclose(start_distance, end_distance):
            raise ValueError("Snapped start and end points are too close on the contour.")
        if end_distance >= start_distance:
            segment = substring(contour, start_distance, end_distance, normalized=False)
            assert isinstance(segment, LineString)
            return segment
        tail = substring(contour, start_distance, length, normalized=False)
        head = substring(contour, 0.0, end_distance, normalized=False)
        assert isinstance(tail, LineString)
        assert isinstance(head, LineString)
        return _concatenate_lines(tail, head)

    start_distance = float(np.clip(start_distance, 0.0, length))
    end_distance = float(np.clip(end_distance, 0.0, length))
    if np.isclose(start_distance, end_distance):
        raise ValueError("Snapped start and end points are too close on the contour.")
    if end_distance < start_distance:
        raise ValueError(
            "Open contour cannot wrap around. Ensure start precedes end in contour orientation."
        )
    segment = substring(contour, start_distance, end_distance, normalized=False)
    assert isinstance(segment, LineString)
    return segment


@dataclass(frozen=True)
class BoundarySectionGeometry:
    approx_start_xy: np.ndarray
    approx_end_xy: np.ndarray
    snapped_start_xy: np.ndarray
    snapped_end_xy: np.ndarray
    boundary_section: LineString
    thickened_section: BaseGeometry
    crop_rect: tuple[float, float, float, float]


def build_boundary_section_geometry(
    contour: LineString,
    approx_start_xy: tuple[float, float],
    approx_end_xy: tuple[float, float],
    max_distance_px: float,
    crop_margin_px: float,
) -> BoundarySectionGeometry:
    start_d, snapped_start = _project_point_to_line(approx_start_xy, contour)
    end_d, snapped_end = _project_point_to_line(approx_end_xy, contour)
    section = oriented_boundary_section(contour, start_d, end_d)
    thickened = section.buffer(float(max_distance_px), cap_style="round", join_style="round")
    if thickened is None or thickened.is_empty:
        raise ValueError("Failed to create thickened contour section.")
    xmin, ymin, xmax, ymax = thickened.bounds
    margin = float(crop_margin_px)
    crop_rect = (xmin - margin, ymin - margin, xmax + margin, ymax + margin)
    return BoundarySectionGeometry(
        approx_start_xy=np.asarray(approx_start_xy, dtype=float),
        approx_end_xy=np.asarray(approx_end_xy, dtype=float),
        snapped_start_xy=snapped_start,
        snapped_end_xy=snapped_end,
        boundary_section=section,
        thickened_section=thickened,
        crop_rect=crop_rect,
    )


def print_geodataframe_geometry_info(sdata, tissue_mask: str | None) -> None:
    """
    Print basic information about the polygons in a GeoDataFrame.
    """

    if tissue_mask is None:
        print("No tissue mask name provided.")
        return

    print(f"Tissue mask: {tissue_mask}")
    if not tissue_mask in sdata.shapes:
        print(f"{tissue_mask} is not a shapes layer!")

    gdf = sdata.shapes[tissue_mask]

    assert hasattr(gdf, "geometry"), "Input must have a GeoDataFrame-like geometry column."
    geometries = gdf.geometry
    print(f"Number of geometries: {len(geometries)}")

    for i, geom in enumerate(geometries):
        print(f"\nGeometry {i}:")
        print(f"  type: {type(geom)}")
        print(f"  object: {geom}")

        assert isinstance(geom, Polygon)

        print(f"  number of holes: {len(geom.interiors)}")

        if geom.is_empty:
            print("  empty geometry")
        else:
            signed_poly_area = signed_area(geom)
            print(f"  signed area of polygon: {signed_poly_area}")
            assert signed_poly_area >= 0.0 
            # Note: we rely on polygon sign to interpret inside/outside tissue;
            # if needed we could normalize orientation to the expected convention
            # with shapely.geometry.polygon.orient(geom, sign=1.0).
            # Expected convention: the coordinates of the product’s exterior ring are 
            # oriented counter-clockwise and the interior rings (holes) clockwise;
            # the signed area of the polygon is nonnegative.

            xmin, ymin, xmax, ymax = geom.bounds
            bbox_area = float((xmax - xmin) * (ymax - ymin))
            print(f"  bbox: {xmin} <= x <= {xmax}, {ymin} <= y <= {ymax}")
            print(f"  bbox unsigned area: {bbox_area}")


def signed_distance_to_linestring(
    points: np.ndarray,
    line: LineString,
    nearest_on_line: BaseGeometry | None = None,
    closest_xy: np.ndarray | None = None,
    along_distances: np.ndarray | None = None,
) -> np.ndarray:
    """
    Calculates the signed distance to a LineString.
    
    Sign convention:
    Points on the left of the linestring
    (when moving from its first to its last point) will have a negative distance,
    points on the right a positive distance.

    Distance:
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

    # Enforce our sign convention for signed distances:
    # negative on the left, positive on the right of the polyline.
    return dist * (-np.sign(cross))


def closest_points_and_distances(
    points: np.ndarray,  # bin centers
    line: LineString,
    use_signed_distances: bool = False,
) -> "DistanceData":
    print("Calculating distances to and along boundary contour...")
    t_total_start = perf_counter()

    t0 = perf_counter()
    points_xy = np.asarray(points, dtype=float)
    t_points_xy = perf_counter() - t0

    t0 = perf_counter()
    points_geom = shapely.points(points_xy)
    t_points_geom = perf_counter() - t0

    print("  Stage: nearest points (started)")
    t0 = perf_counter()
    _, nearest_on_line = shapely.ops.nearest_points(points_geom, line)
    t_nearest = perf_counter() - t0
    print("  Stage: nearest points (done)")

    t0 = perf_counter()
    closest_xy: np.ndarray = shapely.get_coordinates(nearest_on_line)
    t_get_coords = perf_counter() - t0
    if use_signed_distances:
        print("  Stage: distance along curve (started)")
        t0 = perf_counter()
        along_distances_abs = shapely.line_locate_point(line, nearest_on_line, normalized=False)
        t_along_abs = perf_counter() - t0
        print("  Stage: distance along curve (done)")

        t0 = perf_counter()
        distances = signed_distance_to_linestring(
            points_xy,
            line,
            nearest_on_line=nearest_on_line,
            closest_xy=closest_xy,
            along_distances=along_distances_abs,
        )
        t_distances = perf_counter() - t0
    else:
        t_along_abs = 0.0
        t0 = perf_counter()
        distances = shapely.distance(points_geom, nearest_on_line)
        t_distances = perf_counter() - t0

    if use_signed_distances:
        t0 = perf_counter()
        line_length = float(line.length)
        assert line_length > 0, "line length must be positive."
        along_distances: np.ndarray = along_distances_abs / line_length
        t_along_norm = perf_counter() - t0
    else:
        t0 = perf_counter()
        along_distances = shapely.line_locate_point(line, nearest_on_line, normalized=True)
        t_along_norm = perf_counter() - t0

    assert closest_xy.shape[1] == 2
    assert len(distances.shape) == 1
    assert len(along_distances.shape) == 1
    assert len(points) == len(distances)
    assert len(points) == len(along_distances)
    assert len(points) == closest_xy.shape[0]

    t_total = perf_counter() - t_total_start
    print(
        "Timing (s): "
        f"asarray={t_points_xy:.4f}, "
        f"points={t_points_geom:.4f}, "
        f"nearest={t_nearest:.4f}, "
        f"coords={t_get_coords:.4f}, "
        f"along_abs={t_along_abs:.4f}, "
        f"dist={t_distances:.4f}, "
        f"along_norm={t_along_norm:.4f}, "
        f"total={t_total:.4f}"
    )
    print("...done")
    return DistanceData(
        closest_xy=closest_xy,
        closest_distances=distances,
        closest_along_boundary=along_distances,
    )


def mask_excluding_contour_endpoint_closest_points(
    closest_xy: np.ndarray,
    contour: LineString,
    atol: float = 1e-8,
) -> np.ndarray:
    """
    Avoid edge artifact: keep only points whose closest contour point is not the first/last contour coordinate.
    """
    boundary_coords = np.asarray(contour.coords, dtype=float)
    boundary_start = boundary_coords[0]
    boundary_end = boundary_coords[-1]
    closest_is_start = np.all(np.isclose(closest_xy, boundary_start, atol=atol), axis=1)
    closest_is_end = np.all(np.isclose(closest_xy, boundary_end, atol=atol), axis=1)
    return ~(closest_is_start | closest_is_end)


def get_experiment_by_id(experiments: list[Experiment], dataset_id: str) -> Experiment | None:
    for experiment in experiments:
        if experiment.dataset_id == dataset_id:
            return experiment
    return None

@dataclass
class DistanceData:
    closest_xy: np.ndarray  # closest_xy[i] = (x,y) coordinate of the point on the tissue outline closest to the i-th bin; shape = n x 2 (with n the number of bins)
    closest_distances: np.ndarray  # closest_distances[i] = distance from bin i to the closest point on the tissue outline; shape = n x 1
    closest_along_boundary: np.ndarray  # distance measured along the tissue outline contour between closest_xy[i] and the start of the contour; shape = n x 1


@dataclass
class GeneBinData:
    gene_name: str
    bin_centers: np.ndarray  # bin_centers[i] = (x,y) coordinates a the i-th bin; shape = n x 2 (with n the number of bins)
    transcript_counts: np.ndarray  # transcript_counts[i] = number of transcripts of our gene of interest (ptgds) in bin i; shape = n x 1


@dataclass
class DiffusionBinData:
    # Information about transcripts in bins and distance of bins to tissue boundary,
    # for a subset of the bins (the bins in the neighborhood of the tissue boundary,
    # where we analyze diffusion).
    gene_bins: GeneBinData
    distance_data: DistanceData

    def filter(self, mask: np.ndarray) -> "DiffusionBinData":
        return DiffusionBinData(
            gene_bins=GeneBinData(
                gene_name=self.gene_bins.gene_name,
                bin_centers=self.gene_bins.bin_centers[mask],
                transcript_counts=self.gene_bins.transcript_counts[mask],
            ),
            distance_data=DistanceData(
                closest_xy=self.distance_data.closest_xy[mask],
                closest_distances=self.distance_data.closest_distances[mask],
                closest_along_boundary=self.distance_data.closest_along_boundary[mask],
            ),
        )


@dataclass
class DiffusionResults:
    raw_gene_bins: GeneBinData
    filtered: DiffusionBinData
    analysis_region_bin_centers: np.ndarray
    boundary_section: LineString
    crop_rect: np.ndarray
    approx_start_xy: np.ndarray
    approx_end_xy: np.ndarray
    snapped_start_xy: np.ndarray
    snapped_end_xy: np.ndarray
    max_distance_px: float


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
    assert bin_centers.shape[1] == 2  # (x,y)

    return GeneBinData(
        gene_name=gene_name,
        bin_centers=bin_centers,
        transcript_counts=transcript_counts,
    )


def save_diffusion_results(
    output_path: str,
    raw_gene_bins: GeneBinData,
    analysis_region: DiffusionBinData,
    filtered: DiffusionBinData,
    section: BoundarySectionGeometry,
    max_distance_px: float,
) -> None:
    assert raw_gene_bins.gene_name == filtered.gene_bins.gene_name, (
        "raw_gene_bins and filtered.gene_bins must have the same gene_name."
    )
    assert raw_gene_bins.gene_name == analysis_region.gene_bins.gene_name, (
        "raw_gene_bins and analysis_region.gene_bins must have the same gene_name."
    )
    np.savez_compressed(
        output_path,
        gene_name=np.asarray(raw_gene_bins.gene_name),
        raw_bin_centers=np.asarray(raw_gene_bins.bin_centers, dtype=float),
        raw_transcript_counts=np.asarray(raw_gene_bins.transcript_counts, dtype=float),
        analysis_region_bin_centers=np.asarray(analysis_region.gene_bins.bin_centers, dtype=float),
        filtered_bin_centers=np.asarray(filtered.gene_bins.bin_centers, dtype=float),
        filtered_transcript_counts=np.asarray(filtered.gene_bins.transcript_counts, dtype=float),
        filtered_closest_xy=np.asarray(filtered.distance_data.closest_xy, dtype=float),
        filtered_closest_distances=np.asarray(filtered.distance_data.closest_distances, dtype=float),
        filtered_closest_along_boundary=np.asarray(filtered.distance_data.closest_along_boundary, dtype=float),
        boundary_section_xy=np.asarray(section.boundary_section.coords, dtype=float),
        crop_rect=np.asarray(section.crop_rect, dtype=float),
        approx_start_xy=np.asarray(section.approx_start_xy, dtype=float),
        approx_end_xy=np.asarray(section.approx_end_xy, dtype=float),
        snapped_start_xy=np.asarray(section.snapped_start_xy, dtype=float),
        snapped_end_xy=np.asarray(section.snapped_end_xy, dtype=float),
        max_distance_px=np.asarray([float(max_distance_px)], dtype=float),
    )
    print(f"Saved diffusion results to {output_path}")


def load_diffusion_results(output_path: str) -> DiffusionResults:
    with np.load(output_path) as data:
        gene_name = str(np.asarray(data["gene_name"]).item())
        if "analysis_region_bin_centers" not in data.files:
            raise ValueError(
                f"Missing analysis_region_bin_centers in {output_path}. "
                "Regenerate this file with compute_diffusion()."
            )
        return DiffusionResults(
            raw_gene_bins=GeneBinData(
                gene_name=gene_name,
                bin_centers=np.asarray(data["raw_bin_centers"], dtype=float),
                transcript_counts=np.asarray(data["raw_transcript_counts"], dtype=float),
            ),
            filtered=DiffusionBinData(
                gene_bins=GeneBinData(
                    gene_name=gene_name,
                    bin_centers=np.asarray(data["filtered_bin_centers"], dtype=float),
                    transcript_counts=np.asarray(data["filtered_transcript_counts"], dtype=float),
                ),
                distance_data=DistanceData(
                    closest_xy=np.asarray(data["filtered_closest_xy"], dtype=float),
                    closest_distances=np.asarray(data["filtered_closest_distances"], dtype=float),
                    closest_along_boundary=np.asarray(data["filtered_closest_along_boundary"], dtype=float),
                ),
            ),
            analysis_region_bin_centers=np.asarray(data["analysis_region_bin_centers"], dtype=float),
            boundary_section=LineString(np.asarray(data["boundary_section_xy"], dtype=float)),
            crop_rect=np.asarray(data["crop_rect"], dtype=float),
            approx_start_xy=np.asarray(data["approx_start_xy"], dtype=float),
            approx_end_xy=np.asarray(data["approx_end_xy"], dtype=float),
            snapped_start_xy=np.asarray(data["snapped_start_xy"], dtype=float),
            snapped_end_xy=np.asarray(data["snapped_end_xy"], dtype=float),
            max_distance_px=float(np.asarray(data["max_distance_px"], dtype=float).reshape(-1)[0]),
        )


def filter_bins_for_diffusion_analysis(
    diffusion_data: DiffusionBinData,
    boundary_section: LineString,
    max_distance_px: float,
) -> DiffusionBinData:
    """
    Keep bins eligible for diffusion analysis before transcript-count filtering.

    This mirrors compute_diffusion(): stay within the round max-distance neighborhood and
    remove bins whose closest boundary point is exactly one of the section endpoints.
    """
    points_geom = shapely.points(diffusion_data.gene_bins.bin_centers)
    thickened_section = boundary_section.buffer(
        float(max_distance_px),
        cap_style="round",
        join_style="round",
    )
    inside_thickened_mask = np.asarray(
        shapely.covers(thickened_section, points_geom),
        dtype=bool,
    )
    filtered = diffusion_data.filter(inside_thickened_mask)

    # Numerical guard: keep expected max-distance neighborhood only.
    filtered = filtered.filter(
        np.abs(filtered.distance_data.closest_distances) <= (max_distance_px + 1e-6)
    )

    endpoint_mask = mask_excluding_contour_endpoint_closest_points(
        filtered.distance_data.closest_xy,
        boundary_section,
    )
    return filtered.filter(endpoint_mask)


def compute_diffusion(
    sdata_folder: str,
    experiment: Experiment,
    level: str,
    gene_name: str,
    plots_folder: str,
    max_distance_px: float = 400,
    crop_margin_px: float = 0.0,
    use_signed_distances: bool = False,
) -> None:
    print(f"Computing diffusion for {experiment.technology} {experiment.dataset_id} at level {level}")

    sdata = sd.read_zarr(os.path.join(sdata_folder, experiment.sdata_path))

    if experiment.tissue_mask is None:
        raise ValueError("tissue_mask must be set when using boundary_point_pairs.")
    if not experiment.boundary_point_pairs:
        raise ValueError(
            f"No boundary_point_pairs provided for experiment {experiment.dataset_id!r}."
        )

    contour = longest_shape_contour(sdata, experiment.tissue_mask)
    if contour is None or contour.is_empty:
        raise ValueError(f"No tissue contour found in shapes layer {experiment.tissue_mask!r}.")

    for crop_idx, (approx_start, approx_end) in enumerate(experiment.boundary_point_pairs):
        section = build_boundary_section_geometry(
            contour=contour,
            approx_start_xy=approx_start,
            approx_end_xy=approx_end,
            max_distance_px=max_distance_px,
            crop_margin_px=crop_margin_px,
        )
        crop_rect = np.asarray(section.crop_rect, dtype=float)
        print(
            f"Using computed crop{crop_idx + 1} rect: "
            f"{tuple(np.round(crop_rect, 3).tolist())}"
        )

        sdata_crop = sdata.query.bounding_box(
            min_coordinate=crop_rect[:2],
            max_coordinate=crop_rect[2:],
            axes=("x", "y"),
            target_coordinate_system=experiment.dataset_id,
        )

        plot_folder = os.path.join(
            plots_folder,
            experiment.technology,
            experiment.dataset_id,
            "analysis",
            f"crop{crop_idx + 1}",
        )
        os.makedirs(plot_folder, exist_ok=True)

        gene_bin_data = get_gene_bin_data(sdata_crop, level, gene_name)

        distance_data = closest_points_and_distances(
            gene_bin_data.bin_centers,
            section.boundary_section,
            use_signed_distances=use_signed_distances,
        )
        diffusion_data = DiffusionBinData(
            distance_data=distance_data,
            gene_bins=gene_bin_data,
        )
        analysis_region = filter_bins_for_diffusion_analysis(
            diffusion_data,
            section.boundary_section,
            max_distance_px,
        )
        filtered = analysis_region.filter(analysis_region.gene_bins.transcript_counts > 0)

        diffusion_results_output_path = os.path.join(plot_folder, "diffusion.npz")
        save_diffusion_results(
            diffusion_results_output_path,
            raw_gene_bins=gene_bin_data,
            analysis_region=analysis_region,
            filtered=filtered,
            section=section,
            max_distance_px=max_distance_px,
        )


def generate_diffusion_plots(
    sdata_folder: str,
    experiment: Experiment,
    level: str,
    plots_folder: str,
    plot_settings: DiffusionPlotSettings,
) -> None:
    print(f"Generating diffusion plots for {experiment.technology} {experiment.dataset_id} at level {level}")

    sdata = None
    if (
        (plot_settings.tissue_outline and experiment.he_image)
        or plot_settings.full_tissue_bin_gene_counts
    ):
        sdata = sd.read_zarr(os.path.join(sdata_folder, experiment.sdata_path))

    if not experiment.boundary_point_pairs:
        raise ValueError(
            f"No boundary_point_pairs provided for experiment {experiment.dataset_id!r}."
        )

    full_tissue_boundary = None
    full_tissue_gene_bin_data = None
    if plot_settings.full_tissue_bin_gene_counts:
        assert sdata is not None
        if experiment.tissue_mask is None:
            raise ValueError("tissue_mask must be set to generate full-tissue diffusion plots.")
        full_tissue_boundary = full_shape_boundary_geometry(sdata, experiment.tissue_mask)
        if full_tissue_boundary is None or full_tissue_boundary.is_empty:
            raise ValueError(
                f"No tissue boundary found in shapes layer {experiment.tissue_mask!r}."
            )

    for crop_idx, _ in enumerate(experiment.boundary_point_pairs):
        plot_folder = os.path.join(
            plots_folder,
            experiment.technology,
            experiment.dataset_id,
            "analysis",
            f"crop{crop_idx + 1}",
        )
        os.makedirs(plot_folder, exist_ok=True)

        diffusion_results_path = os.path.join(plot_folder, "diffusion.npz")
        if not os.path.exists(diffusion_results_path):
            raise FileNotFoundError(
                f"Missing diffusion results file: {diffusion_results_path}. "
                "Run compute_diffusion() first."
            )
        results = load_diffusion_results(diffusion_results_path)
        crop_rect = np.asarray(results.crop_rect, dtype=float)
        analysis_region_bin_centers = results.analysis_region_bin_centers
        if plot_settings.full_tissue_bin_gene_counts and full_tissue_gene_bin_data is None:
            assert sdata is not None
            full_tissue_gene_bin_data = get_gene_bin_data(
                sdata,
                level,
                results.raw_gene_bins.gene_name,
            )

        if plot_settings.tissue_outline and experiment.he_image:
            assert sdata is not None
            sdata_crop = sdata.query.bounding_box(
                min_coordinate=crop_rect[:2],
                max_coordinate=crop_rect[2:],
                axes=("x", "y"),
                target_coordinate_system=experiment.dataset_id,
            )
            plot_tissue_outline_on_HE_image(
                sdata_crop,
                experiment.he_image,
                experiment.tissue_mask,
                experiment.technology,
                experiment.dataset_id,
                experiment.pixel_size_um,
                os.path.join(plot_folder, "tissue.png"),
                fill_alpha=0.0,
                outline_width=1,
                bounds=tuple(crop_rect),
                overlay_contour=results.boundary_section,
                overlay_bin_centers=analysis_region_bin_centers,
                approx_points=np.vstack([results.approx_start_xy, results.approx_end_xy]),
                snapped_points=np.vstack([results.snapped_start_xy, results.snapped_end_xy]),
            )

        if plot_settings.boundary:
            plot_boundary_geoseries(
                results.boundary_section,
                output_path=os.path.join(plot_folder, "boundary.png"),
            )

        if plot_settings.bin_gene_counts:
            plot_bin_gene_counts_and_contour(
                results.raw_gene_bins,
                results.boundary_section,
                os.path.join(plot_folder, "transcripts_closeup.png"),
                experiment.pixel_size_um,
                f"{experiment.technology}",
                dpi=600,
            )

        if plot_settings.full_tissue_bin_gene_counts:
            assert full_tissue_boundary is not None
            assert full_tissue_gene_bin_data is not None
            plot_full_tissue_bin_gene_counts_with_highlighted_section(
                full_tissue_gene_bin_data,
                full_tissue_boundary,
                results.boundary_section,
                analysis_region_bin_centers,
                os.path.join(plot_folder, "transcripts.png"),
                experiment.pixel_size_um,
                f"{experiment.technology}",
                dpi=600,
            )

        if plot_settings.bins_and_closest_points:
            plot_bins_and_closest_points_on_contour(
                results.filtered,
                results.boundary_section,
                os.path.join(plot_folder, "closest_points.png"),
                f"{experiment.technology}\nBins and closest point on boundary\n",
                line_width=0.5,
                dpi=600,
            )

        if plot_settings.distances_linearized:
            plot_distances_linearized(
                results.filtered,
                os.path.join(plot_folder, "distances_linearized.png"),
                f"{experiment.technology}",
                experiment.pixel_size_um,
                dpi=600,
            )

    # Spatial inspection plots remain per section; distribution plots pool all sections.
    if plot_settings.histogram or plot_settings.boxplot:
        crop_specs = [
            (experiment.technology, experiment.dataset_id, experiment.pixel_size_um, f"crop{i + 1}")
            for i in range(len(experiment.boundary_point_pairs))
        ]
        pooled = read_experiment_distance_data(crop_specs, plots_folder)[0]
        pooled_folder = os.path.join(
            plots_folder, experiment.technology, experiment.dataset_id, "analysis", "pooled",
        )
        os.makedirs(pooled_folder, exist_ok=True)

    if plot_settings.histogram:
        bin_widths_um = plot_settings.histogram_bin_widths_um
        if len(bin_widths_um) == 0:
            raise ValueError("histogram_bin_widths_um must contain at least one bin width.")

        for bin_width_um in bin_widths_um:
            bin_width_str = f"{bin_width_um:g}".replace(".", "p")
            histogram_output_path = os.path.join(
                pooled_folder,
                f"histogram_binwidth_{bin_width_str}um.png",
            )

            plot_distance_histogram(
                pooled,
                histogram_output_path,
                gene_name=results.raw_gene_bins.gene_name,
                title=f"{experiment.technology}",
                bin_width_um=bin_width_um,
                max_distance_um=plot_settings.max_plot_distance_um,
            )

    if plot_settings.boxplot:
        plot_distance_boxplot(
            pooled,
            os.path.join(pooled_folder, "boxplot_um.png"),
            title=f"{experiment.technology}"
        )
