"""Regression checks for pooling disjoint boundary sections per experiment.

Run from diffusion/: .venv/bin/python -m unittest discover -s tests -v
"""
from dataclasses import asdict
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

import diffusion_plots as plots
from diffusion import Experiment, LevelInfo
from diffusion_utils import generate_diffusion_plots


class ExperimentPoolingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(self.temp.name)
        self.distances = np.array([-6., -4., -2., -1., 1., 2., 4., 6.])
        self.weights = np.array([1., 2., 5., 9., 7., 5., 2., 1.])

    def crop(self, dataset, crop, distances, weights, pixel_size=0.5):
        folder = self.root / "Nova-ST" / dataset / "analysis" / crop
        folder.mkdir(parents=True, exist_ok=True)
        distances, weights = np.asarray(distances), np.asarray(weights)
        centers = np.column_stack((np.arange(distances.size), distances))
        np.savez_compressed(
            folder / "diffusion.npz",
            gene_name=np.asarray("Ptgds"),
            raw_bin_centers=centers,
            raw_transcript_counts=weights,
            analysis_region_bin_centers=centers,
            filtered_bin_centers=centers,
            filtered_transcript_counts=weights,
            filtered_closest_xy=np.column_stack((np.arange(distances.size), np.zeros(distances.size))),
            filtered_closest_distances=distances,
            filtered_closest_along_boundary=np.linspace(0, 1, distances.size),
            boundary_section_xy=np.array([[0., 0.], [10., 0.]]),
            crop_rect=np.array([-10., -10., 10., 10.]),
            approx_start_xy=np.array([0., 0.]),
            approx_end_xy=np.array([10., 0.]),
            snapped_start_xy=np.array([0., 0.]),
            snapped_end_xy=np.array([10., 0.]),
            max_distance_px=np.array([20.]),
        )
        return ("Nova-ST", dataset, pixel_size, crop)

    def split_and_reference(self):
        specs = [
            self.crop("sample", "crop1", self.distances, self.weights),
            self.crop("sample", "crop2", self.distances + 3, self.weights * 3),
        ]
        reference = [self.crop(
            "reference", "crop1",
            np.concatenate([self.distances, self.distances + 3]),
            np.concatenate([self.weights, self.weights * 3]),
        )]
        return specs, reference

    def test_single_section_preserves_distances_weights_and_label(self):
        spec = self.crop("sample", "crop1", self.distances, self.weights)
        result, = plots.read_experiment_distance_data([spec], self.root)
        self.assertEqual(result.label, "Nova-ST")
        self.assertEqual((result.technology, result.dataset_id), ("Nova-ST", "sample"))
        assert_array_equal(result.distances_um, self.distances * 0.5)
        assert_array_equal(result.weights, self.weights)

    def test_pool_before_normalization_peak_alignment_and_statistics(self):
        specs, reference = self.split_and_reference()
        pooled, = plots.read_experiment_distance_data(specs, self.root)
        expected, = plots.read_experiment_distance_data(reference, self.root)
        assert_array_equal(pooled.distances_um, expected.distances_um)
        assert_array_equal(pooled.weights, expected.weights)
        for align in [False, True]:
            actual_curve, = plots._distance_kde_curves_from_crop_specs(specs, self.root, align)
            reference_curve, = plots._distance_kde_curves_from_crop_specs(reference, self.root, align)
            assert_array_equal(actual_curve.x_um, reference_curve.x_um)
            assert_array_equal(actual_curve.density, reference_curve.density)
            self.assertEqual(actual_curve.alignment_shift_um, reference_curve.alignment_shift_um)
        for calculate in [
            plots.calculate_distance_histogram_hwhm,
            plots.calculate_distance_to_kde_peak_median,
            plots.summarize_distance_distributions,
        ]:
            actual, = calculate(specs, self.root)
            expected, = calculate(reference, self.root)
            self.assertEqual(asdict(actual), asdict(expected))
        summary, = plots.summarize_distance_distributions(specs, self.root)
        expected_outside = pooled.weights[pooled.distances_um <= 0].sum() / pooled.weights.sum()
        self.assertEqual(summary.side_mass_fraction["outside_fraction"], expected_outside)
        # CDF ordinates and peak-based origins also come from the pooled observations.
        for origin in ["boundary", "peak"]:
            captured = []
            def capture(fig, output_path):
                captured.append([[line.get_data() for line in ax.lines] for ax in fig.axes])
                plt.close(fig)
            with patch.object(plots, "save_or_show_figure", side_effect=capture):
                for selection in [specs, reference]:
                    plots.plot_distance_cdf(selection, str(self.root / "cdf.png"), self.root, origin=origin)
            for actual_axes, expected_axes in zip(*captured, strict=True):
                self.assertEqual(len(actual_axes), len(expected_axes))
                for actual_line, expected_line in zip(actual_axes, expected_axes, strict=True):
                    assert_allclose(actual_line[0], expected_line[0])
                    assert_allclose(actual_line[1], expected_line[1])

    def test_same_technology_experiments_remain_separate_and_render(self):
        specs, _ = self.split_and_reference()
        specs.insert(1, self.crop("second", "crop1", self.distances - 2, self.weights))
        distributions = plots.read_experiment_distance_data(specs, self.root)
        self.assertEqual([d.label for d in distributions], ["Nova-ST (sample)", "Nova-ST (second)"])
        self.assertEqual([d.weights.sum() for d in distributions], [4 * self.weights.sum(), self.weights.sum()])
        for fn in [plots.plot_distance_kde_combined, plots.plot_distance_kde_ridgeline]:
            fn(specs, str(self.root / f"{fn.__name__}.png"), self.root, "Ptgds", align_peaks=True)
        plots.plot_distance_cdf(specs, str(self.root / "cdf.png"), self.root, origin="peak")
        hwhm = plots.calculate_distance_histogram_hwhm(specs, self.root)
        medians = plots.calculate_distance_to_kde_peak_median(specs, self.root)
        plots.plot_distance_histogram_hwhm_bars(hwhm, str(self.root / "hwhm.png"))
        plots.plot_distance_to_kde_peak_median_bars(medians, str(self.root / "medians.png"))
        self.assertEqual(len(list(self.root.glob("*.png"))), 5)

    def test_only_configured_crops_are_loaded_and_empty_sections_are_allowed(self):
        first = self.crop("sample", "crop1", self.distances, self.weights)
        empty = self.crop("sample", "crop2", [], [])
        self.crop("sample", "crop3", [1000.], [999.])  # stale/unselected results
        result, = plots.read_experiment_distance_data([first, empty], self.root)
        assert_array_equal(result.distances_um, self.distances * 0.5)
        with self.assertRaisesRegex(ValueError, "No distance rows found for experiment"):
            plots.read_experiment_distance_data([empty], self.root)
        with self.assertRaisesRegex(FileNotFoundError, "compute_diffusion"):
            plots.read_experiment_distance_data([first, (*first[:3], "crop4")], self.root)

    def test_spatial_plots_stay_per_crop_and_distributions_are_pooled(self):
        self.split_and_reference()
        level = LevelInfo("level", "table", "labels")
        experiment = Experiment(
            "Nova-ST", "brain", "sample", "unused.zarr", None, 0.5, "mask", [level], level,
            boundary_point_pairs=[((0., 0.), (10., 0.)), ((100., 0.), (110., 0.))],
        )
        settings = plots.DiffusionPlotSettings(
            tissue_outline=False, boundary=True, bin_gene_counts=False,
            full_tissue_bin_gene_counts=False, bins_and_closest_points=False,
            distances_linearized=False, histogram=True, histogram_bin_widths_um=[2.], boxplot=True,
        )
        generate_diffusion_plots("unused", experiment, "table", str(self.root), settings)
        analysis = self.root / "Nova-ST" / "sample" / "analysis"
        for crop in ["crop1", "crop2"]:
            self.assertTrue((analysis / crop / "diffusion.npz").is_file())
            self.assertTrue((analysis / crop / "boundary.png").is_file())
            self.assertFalse((analysis / crop / "boxplot_um.png").exists())
        self.assertTrue((analysis / "pooled" / "histogram_binwidth_2um.png").is_file())
        self.assertTrue((analysis / "pooled" / "boxplot_um.png").is_file())


if __name__ == "__main__":
    unittest.main()
