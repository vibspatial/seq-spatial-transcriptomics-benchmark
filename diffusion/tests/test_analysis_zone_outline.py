"""Analysis plots must retain the inner edge of a closed distance zone."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import LineString, MultiLineString, Point, Polygon

from diffusion import Experiment, LevelInfo
import diffusion_plots as plots
import diffusion_utils as utils


class AnalysisZoneOutlineTests(unittest.TestCase):
    def setUp(self):
        self.ring = LineString([(0, 0), (0, 20), (20, 20), (20, 0), (0, 0)])
        self.addCleanup(plt.close, "all")

    def render_zone(self, section, radius):
        x, y = np.meshgrid(np.arange(-12, 33), np.arange(-12, 33))
        points = np.column_stack((x.ravel(), y.ravel()))
        selected = utils.diffusion_analysis_mask(points, section, radius)
        bins = utils.GeneBinData("Ptgds", points, np.ones(len(points)))
        result = SimpleNamespace(
            crop_rect=np.array([-12, -12, 32, 32]),
            boundary_section=section,
            max_distance_px=radius,
            analysis_region_bin_centers=points[selected],
            raw_gene_bins=bins,
            approx_start_xy=np.asarray(section.coords[0]),
            approx_end_xy=np.asarray(section.coords[-1]),
            snapped_start_xy=np.asarray(section.coords[0]),
            snapped_end_xy=np.asarray(section.coords[-1]),
        )
        level = LevelInfo("level", "table", "labels")
        experiment = Experiment(
            "test", "brain", "sample", "unused.zarr", "image", 1., "mask",
            [level], level, boundary_point_pairs=[(section.coords[0], section.coords[-1])],
        )
        sdata = SimpleNamespace(shapes={"mask": gpd.GeoDataFrame(geometry=[Polygon(self.ring)])})
        sdata.query = SimpleNamespace(bounding_box=lambda **kwargs: sdata)
        settings = plots.DiffusionPlotSettings(
            tissue_outline=True, boundary=False, bin_gene_counts=False,
            full_tissue_bin_gene_counts=True, bins_and_closest_points=False,
            distances_linearized=False, histogram=False, boxplot=False,
        )
        root = self.enterContext(tempfile.TemporaryDirectory())
        saved_result = Path(root) / "test/sample/analysis/crop1/diffusion.npz"
        saved_result.parent.mkdir(parents=True)
        saved_result.touch()
        with patch.object(utils.sd, "read_zarr", return_value=sdata), \
             patch.object(utils, "load_diffusion_results", return_value=result), \
             patch.object(utils, "get_gene_bin_data", return_value=bins), \
             patch.object(utils, "plot_tissue_outline_on_HE_image") as he_plot, \
             patch.object(plots, "save_or_show_figure") as save:
            utils.generate_diffusion_plots("unused", experiment, "table", root, settings)
        ax = save.call_args.args[0].axes[0]
        segments = [segment for collection in ax.collections
                    if isinstance(collection, LineCollection) and collection.get_zorder() == 2
                    for segment in collection.get_segments()]
        return MultiLineString(segments), he_plot.call_args.kwargs

    def test_closed_zone_draws_inside_and_outside_edges_at_the_distance_limit(self):
        outline, he_options = self.render_zone(self.ring, 2.)
        self.assertEqual(len(outline.geoms), 2)
        self.assertTrue(outline.equals(self.ring.buffer(2.).boundary))
        self.assertAlmostEqual(outline.distance(Point(10, 10)), 8.)
        self.assertIsNone(he_options["overlay_bin_centers"])
        self.assertTrue(he_options["overlay_region"].boundary.equals(outline))
        # Each ring bounds the eligible zone, rather than tracing the tissue edge.
        for edge in outline.geoms:
            distances = shapely.distance(shapely.points(np.asarray(edge.coords)), self.ring)
            np.testing.assert_allclose(distances, 2., atol=1e-10)

    def test_closed_zone_has_no_inner_edge_when_distance_reaches_the_whole_interior(self):
        outline, _ = self.render_zone(self.ring, 11.)
        self.assertEqual(len(outline.geoms), 1)
        self.assertTrue(outline.equals(self.ring.buffer(11.).boundary))

    def test_open_section_keeps_bin_based_outline(self):
        section = LineString([(0, 0), (0, 20)])
        outline, he_options = self.render_zone(section, 2.)
        self.assertEqual(len(outline.geoms), 1)
        self.assertIsNone(he_options["overlay_region"])
        self.assertGreater(len(he_options["overlay_bin_centers"]), 0)
        xmin, ymin, xmax, ymax = outline.bounds
        self.assertLess(xmin, -2.)
        self.assertGreater(xmax, 2.)
        self.assertGreater(ymin, 0.)
        self.assertLess(ymax, 20.)


if __name__ == "__main__":
    unittest.main()
