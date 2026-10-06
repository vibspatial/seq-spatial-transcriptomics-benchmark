"""The preview and analysis must apply the same endpoint-aware spatial selection."""
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_array_equal
import shapely
from shapely.geometry import LineString, Point

import diffusion_utils as utils


class AnalysisRegionPreviewTests(unittest.TestCase):
    def reference_mask(self, points, section, radius):
        """Original analysis rule, evaluating nearest points everywhere."""
        distances = utils.closest_points_and_distances(points, section, use_signed_distances=True)
        inside = shapely.covers(section.buffer(radius, cap_style="round", join_style="round"),
                               shapely.points(points))
        return (inside & (np.abs(distances.closest_distances) <= radius + 1e-6)
                & utils.mask_excluding_contour_endpoint_closest_points(distances.closest_xy, section))

    def test_straight_section_excludes_caps_and_the_endpoint_cross_sections(self):
        section = LineString([(0, 0), (0, 10)])
        points = np.array([[0, -1], [1, -1], [0, 0], [1, 0],
                           [1, 0.1], [-1, 5], [1, 9.9], [1, 10], [0, 11], [3, 5]])
        expected = [False, False, False, False, True, True, True, False, False, False]
        assert_array_equal(utils.diffusion_analysis_mask(points, section, 2), expected)
        assert_array_equal(self.reference_mask(points, section, 2), expected)

    def test_flat_buffer_is_not_equivalent_for_a_folded_section(self):
        section = LineString([(0, 0), (10, 0), (10, 10), (-2, 10), (-2, -2)])
        # Another part of the line's flat buffer covers this point, but its nearest
        # point is the start endpoint, so the actual analysis must exclude it.
        self.assertTrue(section.buffer(3, cap_style="flat").covers(Point(0, -1)))
        points = np.array([[0, -1], [-1.5, -1], [5, 1]])
        expected = [False, True, True]
        assert_array_equal(utils.diffusion_analysis_mask(points, section, 3), expected)
        assert_array_equal(self.reference_mask(points, section, 3), expected)

    def test_sampled_preview_matches_original_rule_and_actual_bin_filter(self):
        sections = [
            LineString([(0, 0), (0, 10)]),
            LineString([(0, 0), (10, 0), (10, 10), (-2, 10), (-2, -2)]),
            LineString([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)]),
        ]
        for section in sections:
            with self.subTest(closed=section.is_closed, coords=list(section.coords)):
                region = utils.build_analysis_region_preview(section, 3, max_grid_size=96)
                xx, yy = np.meshgrid(region.x_px, region.y_px)
                points = np.column_stack((xx.ravel(), yy.ravel()))
                expected = self.reference_mask(points, section, 3)
                assert_array_equal(region.eligible.ravel(), expected)
                self.assertFalse(region.eligible[0].any())
                self.assertFalse(region.eligible[-1].any())
                distances = utils.closest_points_and_distances(points, section, use_signed_distances=True)
                data = utils.DiffusionBinData(utils.GeneBinData("Ptgds", points, np.ones(len(points))), distances)
                selected = utils.filter_bins_for_diffusion_analysis(data, section, 3)
                assert_array_equal(selected.gene_bins.bin_centers, points[expected])

    def test_endpoint_neighborhood_optimization_preserves_coordinate_tolerance(self):
        section = LineString([(10000, 10000), (10000, 10010)])
        points = np.array([[10001, 10000.05], [10001, 10005], [9999, 10009.95]])
        # Existing np.isclose endpoint tolerances also exclude these near-endpoint projections.
        expected = [False, True, False]
        assert_array_equal(utils.diffusion_analysis_mask(points, section, 1), expected)
        assert_array_equal(self.reference_mask(points, section, 1), expected)

    def test_only_points_near_open_endpoints_need_projection_for_the_preview(self):
        section = LineString([(0, 0), (0, 1000)])
        points = np.array([[1, -1], [1, 500], [1, 1001], [5, 500]])
        nearest = shapely.ops.nearest_points
        with patch.object(shapely.ops, "nearest_points", wraps=nearest) as query:
            assert_array_equal(utils.diffusion_analysis_mask(points, section, 2), [False, True, False, False])
        query.assert_called_once()
        self.assertEqual(len(query.call_args.args[0]), 2)
        ring = LineString([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
        with patch.object(shapely.ops, "nearest_points", side_effect=AssertionError("Closed ring needs no projection")):
            assert_array_equal(utils.diffusion_analysis_mask(np.array([[-1, -1], [5, 5]]), ring, 2), [True, False])


if __name__ == "__main__":
    unittest.main()
