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
        right = edge["value"]["width_pixels"]["right"]
        self.assertEqual(right["invalid_width"], 3)
        self.assertAlmostEqual(right["noisy_width"], 15, delta=2)

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

    def test_injected_tile_seam_is_flagged(self):
        unwrap = self.grids["unwrap.grd"].copy()
        unwrap[self.ny // 2 :, :] += np.float32(2.0 * np.pi)
        self._write_grid(self.intf_dir / "unwrap.grd", unwrap)
        result, report = self._run_qa()
        self.assertEqual(result.returncode, 0, result.stderr)
        seam = self._check(report, "tile_seams")
        self.assertIn(seam["status"], ("WARN", "FAIL"))
        self.assertGreater(seam["value"]["row_seams"][0]["ratio"], 1.5)

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
