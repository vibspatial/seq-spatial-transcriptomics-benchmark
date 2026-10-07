"""Select and analyse sections on separate tissue exteriors in one experiment."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import Mock, patch

import geopandas as gpd
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import MultiPolygon, Polygon

from diffusion import Experiment, LevelInfo
import diffusion_utils as utils


class MultipleContourTests(unittest.TestCase):
    def setUp(self):
        self.large = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)],
                             [[(2, 2), (4, 2), (4, 4), (2, 4)]])
        self.small = Polygon([(30, 0), (30, 4), (34, 4), (34, 0)])

    def sdata(self, geometries):
        return SimpleNamespace(shapes={"mask": gpd.GeoDataFrame(geometry=geometries)})

    def test_selects_each_exterior_including_multipolygon_parts_and_ignores_holes(self):
        for geometries in [[self.large, self.small], [MultiPolygon([self.large, self.small])]]:
            with self.subTest(multipart=len(geometries) == 1):
                sdata = self.sdata(geometries)
                for start, end, polygon in [((2, 3), (0, 8), self.large),
                                             ((29, 1), (29, 3), self.small)]:
                    contour = utils.shape_contour_for_points(sdata, "mask", start, end)
                    assert_array_equal(contour.coords, polygon.exterior.coords)
                contour = utils.shape_contour_for_points(sdata, "mask", (29, 2), (29, 2))
                section = utils.build_boundary_section_geometry(contour, (29, 2), (29, 2), 1, 0)
                self.assertTrue(section.boundary_section.is_closed)
                self.assertEqual(section.boundary_section.length, self.small.length)
                assert_array_equal(section.snapped_start_xy, [30, 2])

    def test_rejects_endpoints_on_different_contours(self):
        with self.assertRaisesRegex(ValueError, "selects different exterior contours"):
            utils.shape_contour_for_points(self.sdata([self.large, self.small]),
                                           "mask", (0, 2), (30, 2))

    def test_empty_mask_has_clear_error(self):
        for geometries in [[], [Polygon()]]:
            with self.assertRaisesRegex(ValueError, "No exterior tissue contour"):
                utils.shape_contour_for_points(self.sdata(geometries), "mask", (0, 0), (0, 0))

    def test_exact_distance_ties_use_shape_order(self):
        for geometries in [[self.large, self.small], [self.small, self.large]]:
            contour = utils.shape_contour_for_points(self.sdata(geometries), "mask", (20, 2), (20, 2))
            assert_array_equal(contour.coords, geometries[0].exterior.coords)

    def test_compute_saves_sections_from_both_geometries_and_pools_results(self):
        sdata = self.sdata([self.large, self.small])
        sdata.query = SimpleNamespace(bounding_box=Mock(return_value=sdata))
        points = np.array([[-1., 5.], [1., 5.], [29., 2.], [31., 2.]])
        gene_bins = utils.GeneBinData("Ptgds", points, np.array([1., 2., 3., 4.]))
        level = LevelInfo("level", "table", "labels")
        experiment = Experiment(
            "Nova-ST", "brain", "sample", "unused.zarr", None, 1., "mask", [level], level,
            boundary_point_pairs=[((0, 2), (0, 8)), ((30, 1), (30, 3))],
        )
        with tempfile.TemporaryDirectory() as root:
            with patch.object(utils.sd, "read_zarr", return_value=sdata), \
                 patch.object(utils, "get_gene_bin_data", return_value=gene_bins):
                utils.compute_diffusion("unused", experiment, "table", "Ptgds", root,
                                        max_distance_px=2., use_signed_distances=True)
            analysis = Path(root) / "Nova-ST" / "sample" / "analysis"
            for i, expected_section in enumerate([[(0, 2), (0, 8)], [(30, 1), (30, 3)]]):
                result = utils.load_diffusion_results(analysis / f"crop{i + 1}" / "diffusion.npz")
                assert_array_equal(result.boundary_section.coords, expected_section)
                assert_array_equal(result.filtered.gene_bins.bin_centers, points[2*i:2*i+2])
                assert_allclose(result.filtered.distance_data.closest_distances, [-1., 1.])
            specs = [("Nova-ST", "sample", 1., f"crop{i}") for i in (1, 2)]
            pooled, = utils.read_experiment_distance_data(specs, root)
            assert_allclose(pooled.distances_um, [-1., 1., -1., 1.])
            assert_array_equal(pooled.weights, gene_bins.transcript_counts)


if __name__ == "__main__":
    unittest.main()
