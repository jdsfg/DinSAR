import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    from osgeo import gdal, osr
except ImportError:
    gdal = None
    osr = None


REPO_ROOT = Path(__file__).resolve().parents[1]
QA_SCRIPT = REPO_ROOT / "dinsar_qa.py"
GMT_AVAILABLE = shutil.which("gmt") is not None


@unittest.skipUnless(GMT_AVAILABLE, "GMT is required for synthetic grid tests")
class DinsarQASyntheticTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory(prefix="dinsar-qa-")
        self.root = Path(self.temp_dir.name)
        self.intf_dir = self.root / "intf"
        self.project_dir = self.root / "project"
        for directory in (
            self.intf_dir,
            self.project_dir / "debug",
            self.project_dir / "topo",
            self.project_dir / "results",
            self.project_dir / "SLC",
        ):
            directory.mkdir(parents=True, exist_ok=True)
        self.nx, self.ny = 200, 120
        y, x = np.mgrid[0:self.ny, 0:self.nx].astype(np.float32)
        phase = 0.002 * x + 0.003 * y
        self.grids = {
            "corr.grd": np.full((self.ny, self.nx), 0.8, dtype=np.float32),
            "phasefilt.grd": phase,
            "unwrap.grd": phase.copy(),
            "unwrap_detrended.grd": phase.copy(),
            "los_displacement.grd": 5.0 * np.sin(x / 20.0) * np.cos(y / 30.0),
        }
        for name, values in self.grids.items():
            self._write_grid(self.intf_dir / name, values)
        self._write_grid(self.project_dir / "debug" / "orbital_ramp.grd", 0.02 * x + 0.04 * y)
        self._write_grid(self.project_dir / "topo" / "dem.grd", np.full_like(x, 100.0))
        self.master_prm = self.project_dir / "SLC" / "master.PRM"
        self.master_prm.write_text("radar_wavelength = 0.0312284\nnear_range = 640074.4387\n", encoding="utf-8")
        (self.project_dir / "SLC" / "coreg_report.json").write_text(
            json.dumps(
                {
                    "gate": {"passed": True},
                    "verify": {
                        "gate": {"passed": True},
                        "check": {
                            "n_valid": 195,
                            "n_points": 195,
                            "frac_within_0p3": 0.98,
                            "robust_sigma_r": 0.02,
                            "robust_sigma_a": 0.03,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        self.snaphu_conf = self.root / "snaphu.conf"
        self.snaphu_conf.write_text("LAMBDA 0.0312284\nNEARRANGE 640074.4387\n", encoding="utf-8")
        self._write_geotiffs()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_grid(self, path, values):
        y, x = np.mgrid[0:self.ny, 0:self.nx]
        points = self.root / "grid.xyz"
        points.write_text(
            "".join(
                f"{xi} {yi} {float(value):.9g}\n"
                for xi, yi, value in zip(x.ravel(), y.ravel(), np.asarray(values).ravel())
            ),
            encoding="ascii",
        )
        subprocess.run(
            [
                "gmt",
                "xyz2grd",
                str(points),
                f"-R0/{self.nx - 1}/0/{self.ny - 1}",
                "-I1",
                f"-G{path}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def _write_geotiffs(self):
        for name in ("unwrap_detrended", "corr", "final_data_values"):
            path = self.project_dir / "results" / f"{name}.tif"
            if gdal is not None and osr is not None:
                dataset = gdal.GetDriverByName("GTiff").Create(
                    str(path), self.nx, self.ny, 1, gdal.GDT_Float32
                )
                dataset.SetGeoTransform((-180.0, 1.0, 0.0, 60.0, 0.0, -1.0))
                spatial_ref = osr.SpatialReference()
                spatial_ref.ImportFromEPSG(4326)
                dataset.SetProjection(spatial_ref.ExportToWkt())
                dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
                dataset.GetRasterBand(1).WriteArray(np.ones((self.ny, self.nx), dtype=np.float32))
                dataset = None
            else:
                source = self.intf_dir / "corr.grd"
                subprocess.run(
                    ["gdal_translate", "-q", "-a_srs", "EPSG:4326", str(source), str(path)],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

    def _run_qa(self, strict=False):
        out_dir = self.root / "qa"
        command = [
            sys.executable,
            str(QA_SCRIPT),
            "--intf-dir",
            str(self.intf_dir),
            "--project-dir",
            str(self.project_dir),
            "--master-prm",
            str(self.master_prm),
            "--wavelength",
            "0.0312284",
            "--downsample",
            "4",
            "--snaphu-tiles",
            "2 2",
            "--snaphu-overlap",
            "10",
            "--snaphu-conf",
            str(self.snaphu_conf),
            "--out-dir",
            str(out_dir),
        ]
        if strict:
            command.append("--strict")
        result = subprocess.run(command, capture_output=True, text=True)
        report_path = out_dir / "qa_report.json"
        report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else None
        return result, report

    @staticmethod
    def _check(report, check_id):
        return next(check for check in report["checks"] if check["id"] == check_id)

    def test_clean_smooth_phase_passes_congruence_and_discontinuity(self):
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self._check(report, "congruence")["status"], "PASS")
        self.assertEqual(self._check(report, "unwrap_discontinuity")["status"], "PASS")
        edge = self._check(report, "edge_strip")
        self.assertEqual(edge["status"], "PASS")
        for widths in edge["value"]["width_pixels"].values():
            self.assertEqual(widths, {"invalid_width": 0, "noisy_width": 0})

    def test_zero_filled_and_noisy_edge_band_is_detected(self):
        phase = self.grids["phasefilt.grd"].copy()
        rng = np.random.default_rng(13)
        phase[:, -18:-3] = rng.uniform(-float(np.pi), float(np.pi), (self.ny, 15))
        phase[:, -3:] = 0.0
        self._write_grid(self.intf_dir / "phasefilt.grd", phase)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        edge = self._check(report, "edge_strip")
        self.assertEqual(edge["status"], "WARN")
        self.assertNotIn("mask_coverage_pct", edge["value"])
        right = edge["value"]["width_pixels"]["right"]
        self.assertEqual(right["invalid_width"], 3)
        self.assertAlmostEqual(right["noisy_width"], 15, delta=2)

    def test_product_mask_reports_noisy_edge_band_coverage(self):
        phase = self.grids["phasefilt.grd"].copy()
        rng = np.random.default_rng(31)
        phase[:, -18:-3] = rng.uniform(-float(np.pi), float(np.pi), (self.ny, 15))
        phase[:, -3:] = 0.0
        self._write_grid(self.intf_dir / "phasefilt.grd", phase)

        mask = np.ones((self.ny, self.nx), dtype=np.float32)
        mask[:, -18:] = np.nan
        self._write_grid(self.intf_dir / "product_mask.grd", mask)
        summary = {
            "fractions": {
                "masked_fraction_by_overlap": 0.1,
                "masked_fraction_by_corr": 0.2,
                "masked_fraction_total": 0.3,
            }
        }
        (self.intf_dir / "product_mask_summary.json").write_text(
            json.dumps(summary), encoding="utf-8"
        )

        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        product_mask = self._check(report, "product_mask")
        self.assertEqual(product_mask["status"], "INFO")
        self.assertEqual(product_mask["value"]["masked_fraction_by_overlap"], 0.1)
        self.assertEqual(product_mask["value"]["masked_fraction_by_corr"], 0.2)
        self.assertEqual(product_mask["value"]["masked_fraction_total"], 0.09)
        self.assertAlmostEqual(
            product_mask["value"]["noisy_band_coverage_by_edge"]["right"]["masked_fraction"],
            1.0,
        )
        edge = self._check(report, "edge_strip")
        self.assertEqual(edge["status"], "PASS")
        self.assertIn("已被产品掩膜覆盖", edge["note_zh"])
        self.assertAlmostEqual(edge["value"]["mask_coverage_pct"]["right"], 100.0)

    def test_partial_product_mask_keeps_noisy_edge_warning(self):
        phase = self.grids["phasefilt.grd"].copy()
        rng = np.random.default_rng(31)
        phase[:, -18:-3] = rng.uniform(-float(np.pi), float(np.pi), (self.ny, 15))
        phase[:, -3:] = 0.0
        self._write_grid(self.intf_dir / "phasefilt.grd", phase)
        mask = np.ones((self.ny, self.nx), dtype=np.float32)
        mask[:, -9:] = np.nan
        self._write_grid(self.intf_dir / "product_mask.grd", mask)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        edge = self._check(report, "edge_strip")
        self.assertEqual(edge["status"], "WARN")
        right_coverage = edge["value"]["mask_coverage_pct"]["right"]
        self.assertGreater(right_coverage, 0.0)
        self.assertLess(right_coverage, 99.0)

    def test_interior_noisy_band_is_not_reported_as_edge_noise(self):
        phase = self.grids["phasefilt.grd"].copy()
        rng = np.random.default_rng(17)
        phase[:, 40:55] = rng.uniform(-float(np.pi), float(np.pi), (self.ny, 15))
        self._write_grid(self.intf_dir / "phasefilt.grd", phase)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        edge = self._check(report, "edge_strip")
        self.assertEqual(edge["status"], "PASS")
        for widths in edge["value"]["width_pixels"].values():
            self.assertEqual(widths["noisy_width"], 0)

    def test_noisy_band_crossing_tile_seam_matches_local_reference(self):
        unwrap = self.grids["unwrap.grd"].copy()
        corr = self.grids["corr.grd"].copy()
        rng = np.random.default_rng(23)
        seam = self.ny // 2
        start, end = seam - 55, seam + 55
        unwrap[start:end, :] = rng.uniform(-float(np.pi), float(np.pi), (end - start, self.nx))
        corr[start:end, :] = 0.35
        self._write_grid(self.intf_dir / "unwrap.grd", unwrap)
        self._write_grid(self.intf_dir / "corr.grd", corr)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        seam_check = self._check(report, "tile_seams")
        row_seam = seam_check["value"]["row_seams"][0]
        self.assertEqual(seam_check["status"], "PASS")
        self.assertGreater(row_seam["fraction"], 0.005)
        self.assertLess(row_seam["ratio"], 1.5)
        self.assertEqual(len(row_seam["local_reference_samples"]), 20)
        self.assertAlmostEqual(
            row_seam["local_reference_fraction"],
            float(np.median([sample["fraction"] for sample in row_seam["local_reference_samples"]])),
        )

    def test_clean_field_step_on_tile_seam_fails(self):
        unwrap = self.grids["unwrap.grd"].copy()
        unwrap[self.ny // 2 :, :] += np.float32(2.0 * np.pi)
        self._write_grid(self.intf_dir / "unwrap.grd", unwrap)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        seam = self._check(report, "tile_seams")
        row_seam = seam["value"]["row_seams"][0]
        self.assertEqual(seam["status"], "FAIL")
        self.assertEqual(row_seam["fraction"], 1.0)
        self.assertEqual(row_seam["local_reference_fraction"], 0.0)
        self.assertGreater(row_seam["ratio"], 3.0)

    def test_low_coherence_tile_seam_pixels_are_excluded(self):
        unwrap = self.grids["unwrap.grd"].copy()
        corr = self.grids["corr.grd"].copy()
        seam_index = self.ny // 2
        unwrap[seam_index:, :] += np.float32(2.0 * np.pi)
        corr[seam_index - 1 : seam_index + 1, :] = 0.2
        self._write_grid(self.intf_dir / "unwrap.grd", unwrap)
        self._write_grid(self.intf_dir / "corr.grd", corr)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        seam = self._check(report, "tile_seams")
        row_seam = seam["value"]["row_seams"][0]
        self.assertEqual(seam["status"], "PASS")
        self.assertEqual(row_seam["fraction"], 0.0)
        self.assertEqual(row_seam["neighbor_pairs"], 0)

    def test_noncongruent_unwrap_fails_strict_mode_with_exit_two(self):
        unwrap = self.grids["unwrap.grd"] + np.float32(0.5)
        self._write_grid(self.intf_dir / "unwrap.grd", unwrap)
        result, report = self._run_qa(strict=True)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertEqual(self._check(report, "congruence")["status"], "FAIL")
        self.assertEqual(report["overall"], "FAIL")

    def test_snaphu_wavelength_mismatch_warns(self):
        self.snaphu_conf.write_text("LAMBDA 0.0565647\nNEARRANGE 640074.4387\n", encoding="utf-8")
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self._check(report, "snaphu_config")["status"], "WARN")

    def _write_coreg(self, **check):
        base = {"n_valid": 100, "n_points": 100, "frac_within_0p3": 0.98,
                "robust_sigma_r": 0.02, "robust_sigma_a": 0.03}
        base.update(check)
        (self.project_dir / "SLC" / "coreg_report.json").write_text(
            json.dumps({"gate": {"passed": True}, "verify": {"gate": {"passed": True}, "check": base}}),
            encoding="utf-8",
        )

    def test_low_fraction_within_0p3_warns(self):
        self._write_coreg(frac_within_0p3=0.8)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self._check(report, "coreg")["status"], "WARN")

    def test_coreg_outlier_is_reported_next_to_rmse(self):
        self._write_coreg(rmse_r=0.39, max_abs_r=5.4, rmse_a=0.07, max_abs_a=0.6)
        result, report = self._run_qa()
        check = self._check(report, "coreg")
        self.assertEqual(check["status"], "PASS")
        self.assertEqual(check["value"]["max_abs_range_px"], 5.4)
        self.assertEqual(check["value"]["rmse_range_px"], 0.39)
        self.assertIn("离群点", check["note_zh"])

    def test_large_removed_ramp_warns(self):
        x = np.arange(self.nx, dtype=np.float32)[None, :] * np.ones((self.ny, 1), dtype=np.float32)
        y = np.arange(self.ny, dtype=np.float32)[:, None] * np.ones((1, self.nx), dtype=np.float32)
        self._write_grid(self.project_dir / "debug" / "orbital_ramp.grd", 0.02 * x + 0.5 * y)
        result, report = self._run_qa()
        check = self._check(report, "detrend")
        self.assertEqual(check["status"], "WARN")
        self.assertGreater(check["value"]["orbital_ramp_peak_to_peak_weeks"], 3.0)

    def test_small_removed_ramp_passes(self):
        result, report = self._run_qa()
        self.assertEqual(self._check(report, "detrend")["status"], "PASS")

    def test_orbit_inconsistency_warns_and_missing_is_skipped(self):
        result, report = self._run_qa()
        self.assertEqual(self._check(report, "orbit")["status"], "SKIP")
        (self.project_dir / "SLC" / "20300105.convert.txt").write_text(
            "led_pv_mismatch_mps = 0.058\nled_pv_drift_m = 0.20\n", encoding="utf-8")
        result, report = self._run_qa()
        self.assertEqual(self._check(report, "orbit")["status"], "WARN")
        (self.project_dir / "SLC" / "20300105.convert.txt").write_text(
            "led_pv_mismatch_mps = 0.005\nled_pv_drift_m = 0.01\n", encoding="utf-8")
        result, report = self._run_qa()
        self.assertEqual(self._check(report, "orbit")["status"], "PASS")

    def test_missing_coreg_report_is_skipped(self):
        (self.project_dir / "SLC" / "coreg_report.json").unlink()
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self._check(report, "coreg")["status"], "SKIP")


class RunDinsarQAOptionTests(unittest.TestCase):
    def _parse(self, args):
        command = (
            "source ./run_dinsar.sh; LOG_FILE=/dev/null; "
            f"parse_args {args} /tmp; printf '%s' \"$QA_MODE\""
        )
        return subprocess.run(
            ["bash", "-c", command],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )

    def test_default_mode_is_on(self):
        result = subprocess.run(
            ["bash", "-c", "source ./run_dinsar.sh; printf '%s' \"$QA_MODE\""],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "on")

    def test_valid_qa_modes_parse(self):
        for mode in ("on", "off", "strict"):
            with self.subTest(mode=mode):
                result = self._parse(f"--qa {mode}")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, mode)

    def test_invalid_qa_mode_is_rejected(self):
        result = self._parse("--qa invalid")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("无效 --qa", result.stdout + result.stderr)
        self.assertNotIn("reached", result.stdout)


if __name__ == "__main__":
    unittest.main()
