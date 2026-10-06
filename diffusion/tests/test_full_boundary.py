"""Full exterior boundary geometry, signed distances, filtering, and output regression checks."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import Mock, patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import LineString, MultiPolygon, Polygon

from diffusion import Experiment, LevelInfo
import diffusion_plots as plots
import diffusion_utils as utils


class FullBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.ring = LineString([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
        self.addCleanup(plt.close, "all")

    def section(self, start, end=None, contour=None):
        return utils.build_boundary_section_geometry(
            self.ring if contour is None else contour,
            start, start if end is None else end,
            max_distance_px=2.0, crop_margin_px=1.0,
        )

    def sdata(self, geometries):
        return SimpleNamespace(shapes={"mask": gpd.GeoDataFrame(geometry=geometries)})

    def test_only_exteriors_are_selected_even_when_a_hole_has_a_longer_perimeter(self):
        # A narrow winding hole has a greater perimeter than the containing square.
        hole = LineString([(1 + 8 * (i % 2), 1 + i * 0.4) for i in range(20)]).buffer(0.1)
        polygon = Polygon(self.ring.coords, [hole.exterior.coords])
        self.assertTrue(polygon.is_valid)
        self.assertGreater(polygon.interiors[0].length, polygon.exterior.length)
        small = Polygon([(20, 0), (20, 2), (22, 2), (22, 0)])
        for geometries in [[polygon, small], [MultiPolygon([small, polygon])]]:
            sdata = self.sdata(geometries)
            selected = utils.longest_shape_contour(sdata, "mask")
            assert_array_equal(selected.coords, self.ring.coords)
            full = utils.full_shape_boundary_geometry(sdata, "mask")
            self.assertEqual(len(full.geoms), 2)
            self.assertAlmostEqual(full.length, 48.0)
            with_holes = utils.full_shape_boundary_geometry(sdata, "mask", include_holes=True)
            self.assertEqual(len(with_holes.geoms), 3)
            self.assertAlmostEqual(with_holes.length, 48.0 + polygon.interiors[0].length)

    def test_full_ring_rotates_to_snapped_start_preserving_all_segments_and_orientation(self):
        section = self.section((-1, 5))
        assert_array_equal(section.boundary_section.coords,
                           [(0, 5), (0, 10), (10, 10), (10, 0), (0, 0), (0, 5)])
        assert_array_equal(section.snapped_start_xy, [0, 5])
        assert_array_equal(section.snapped_end_xy, [0, 5])
        self.assertEqual(section.boundary_section.length, self.ring.length)
        self.assertEqual(utils._signed_ring_area(section.boundary_section.coords),
                         utils._signed_ring_area(self.ring.coords))
        assert_array_equal(section.crop_rect, [-3, -3, 13, 13])
        assert_array_equal(self.section((0, 0)).boundary_section.coords, self.ring.coords)

    def test_signed_distances_are_independent_of_start_and_ring_orientation(self):
        points = np.array([[-1, -1], [0, -1], [-1, 5], [1, 5], [5, 5], [0, 0], [0, 5], [11, 11]])
        expected = np.array([-np.sqrt(2), -1, -1, 1, 5, 0, 0, -np.sqrt(2)])
        for contour in [self.ring, LineString(list(self.ring.coords)[::-1])]:
            for start in [(0, 0), (0, 5), (10, 10), (4.3, 0)]:
                with self.subTest(start=start, signed_area=utils._signed_ring_area(contour.coords)):
                    section = self.section(start, contour=contour).boundary_section
                    self.assertTrue(section.is_closed)
                    self.assertAlmostEqual(section.length, contour.length)
                    self.assertAlmostEqual(utils._signed_ring_area(section.coords),
                                           utils._signed_ring_area(contour.coords))
                    data = utils.closest_points_and_distances(points, section, use_signed_distances=True)
                    assert_allclose(data.closest_distances, expected, atol=1e-12)
                    self.assertTrue(np.all((data.closest_along_boundary >= 0) & (data.closest_along_boundary <= 1)))
        unsigned = utils.closest_points_and_distances(points, self.ring, use_signed_distances=False)
        assert_allclose(unsigned.closest_distances, np.abs(expected))

    def test_full_boundary_keeps_start_point_bins_and_still_applies_distance_limit(self):
        points = np.array([[-1, -1], [0, -1], [-1, 0], [0, 0], [1, 1], [5, 5]])
        distances = utils.closest_points_and_distances(points, self.ring, use_signed_distances=True)
        data = utils.DiffusionBinData(utils.GeneBinData("Ptgds", points, np.ones(len(points))), distances)
        filtered = utils.filter_bins_for_diffusion_analysis(data, self.ring, max_distance_px=2.0)
        assert_array_equal(filtered.gene_bins.bin_centers, points[:-1])
        assert_allclose(filtered.distance_data.closest_distances, [-np.sqrt(2), -1, -1, 0, 1])

    def test_distinct_near_or_identically_snapped_points_still_raise(self):
        for start, end in [((-1, 5), (-2, 5)), ((0, 5), (0, 5 + 1e-7))]:
            with self.subTest(start=start, end=end):
                with self.assertRaisesRegex(ValueError, "too close"):
                    self.section(start, end)
        with self.assertRaisesRegex(ValueError, "closed contour"):
            self.section((0, 5), contour=LineString([(0, 0), (0, 10)]))
        # Large image coordinates must not turn a short open section into a closed ring.
        self.assertFalse(utils._is_closed_linestring(LineString([(10000, 10000), (10000, 10000.01)])))

    def test_open_sections_keep_their_order_endpoint_filter_and_distance_convention(self):
        section = self.section((0, 2), (10, 8)).boundary_section
        assert_array_equal(section.coords, [(0, 2), (0, 10), (10, 10), (10, 8)])
        wrapped = self.section((10, 5), (0, 0)).boundary_section
        assert_array_equal(wrapped.coords, [(10, 5), (10, 0), (0, 0)])
        wrapped = self.section((10, 5), (0, 5)).boundary_section
        assert_array_equal(wrapped.coords, [(10, 5), (10, 0), (0, 0), (0, 5)])
        line = LineString([(0, 0), (0, 10)])
        points = np.array([[1, -1], [-1, 5], [1, 5], [1, 11]])
        distances = utils.closest_points_and_distances(points, line, use_signed_distances=True)
        assert_allclose(distances.closest_distances, [1, -1, 1, 1])
        assert_array_equal(utils.mask_excluding_contour_endpoint_closest_points(distances.closest_xy, line),
                           [False, True, True, False])

    def test_compute_save_load_pool_and_plots_support_the_full_exterior(self):
        polygon = Polygon(self.ring.coords, [[(4, 4), (6, 4), (6, 6), (4, 6)]])
        sdata = self.sdata([polygon])
        sdata.query = SimpleNamespace(bounding_box=Mock(return_value=sdata))
        # The middle bin lies inside the hole but its distance to the exterior is +5.
        points = np.array([[-1, -1], [0, -1], [-1, 5], [1, 5], [5, 5], [0, 0], [10, 9]])
        counts = np.arange(1., len(points) + 1)
        gene_bins = utils.GeneBinData("Ptgds", points, counts)
        level = LevelInfo("level", "table", "labels")
        experiment = Experiment("Nova-ST", "brain", "sample", "unused.zarr", None, 0.5,
                                "mask", [level], level, boundary_point_pairs=[((0, 0), (0, 0))])
        with tempfile.TemporaryDirectory() as root:
            with patch.object(utils.sd, "read_zarr", return_value=sdata), \
                 patch.object(utils, "get_gene_bin_data", return_value=gene_bins):
                utils.compute_diffusion("unused", experiment, "table", "Ptgds", root,
                                        max_distance_px=6.0, crop_margin_px=1.0, use_signed_distances=True)
            analysis = Path(root) / "Nova-ST" / "sample" / "analysis"
            results = utils.load_diffusion_results(analysis / "crop1" / "diffusion.npz")
            self.assertTrue(results.boundary_section.is_closed)
            self.assertEqual(results.boundary_section.length, 40.0)
            expected = np.array([-np.sqrt(2), -1, -1, 1, 5, 0, 0])
            assert_allclose(results.filtered.distance_data.closest_distances, expected)
            assert_array_equal(results.filtered.gene_bins.transcript_counts, counts)
            specs = [("Nova-ST", "sample", 0.5, "crop1")]
            pooled, = plots.read_experiment_distance_data(specs, root)
            assert_allclose(pooled.distances_um, expected * 0.5)
            settings = plots.DiffusionPlotSettings(
                tissue_outline=False, boundary=True, bin_gene_counts=True, full_tissue_bin_gene_counts=True,
                bins_and_closest_points=True, distances_linearized=True,
                histogram=True, histogram_bin_widths_um=[1.0], boxplot=True,
            )
            # Render actual figures at modest test resolution, without changing plot logic.
            save = plots.save_or_show_figure
            def save_small(fig, path, dpi=300):
                fig.set_dpi(50)
                save(fig, path, dpi=50)
            with patch.object(utils.sd, "read_zarr", return_value=sdata), \
                 patch.object(utils, "get_gene_bin_data", return_value=gene_bins), \
                 patch.object(plots, "save_or_show_figure", side_effect=save_small):
                utils.generate_diffusion_plots("unused", experiment, "table", root, settings)
                plots.plot_distance_cdf(specs, str(analysis / "cdf.png"), root, origin="peak")
            self.assertTrue((analysis / "transcripts.png").is_file())
            for name in ["boundary.png", "transcripts_closeup.png",
                         "closest_points.png", "distances_linearized.png"]:
                self.assertTrue((analysis / "crop1" / name).is_file(), name)
            for name in ["histogram_binwidth_1um.png", "boxplot_um.png"]:
                self.assertTrue((analysis / "pooled" / name).is_file(), name)
            self.assertTrue((analysis / "cdf.png").is_file())


if __name__ == "__main__":
    unittest.main()
