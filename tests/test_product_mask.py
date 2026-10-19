#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression tests for the overlap/coherence product mask."""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "product_mask.py"
RUN_DINSAR = REPO_ROOT / "run_dinsar.sh"


@unittest.skipUnless(shutil.which("gmt"), "需要 GMT")
class ProductMaskTests(unittest.TestCase):
    nx = 20
    ny = 12
    dx = 4
    dy = 4
    x_max = 80
    y_max = 48

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.grid = self.root / "corr.grd"
        self.slave_prm = self.root / "slave_shift.PRM"
        self.size_prm = self.root / "slave_size.PRM"
        self.out_mask = self.root / "mask.grd"
        self.summary_path = self.root / "summary.json"
        self.corr = np.ones((self.ny, self.nx), dtype=np.float32)
        self._write_grid(self.corr)
        self._write_prms()

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_grid(self, values):
        xyz = self.root / "corr.xyz"
        with xyz.open("w", encoding="utf-8") as stream:
            for row in range(self.ny):
                y = self.y_max - (row + 0.5) * self.dy
                for col in range(self.nx):
                    x = (col + 0.5) * self.dx
                    stream.write(f"{x} {y} {float(values[row, col])}\n")
        result = subprocess.run(
            [
                "gmt",
                "xyz2grd",
                str(xyz),
                f"-R0/{self.x_max}/0/{self.y_max}",
                f"-I{self.dx}/{self.dy}",
                "-r",
                f"-G{self.grid}",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def _write_prms(self, shift_values=None, nx_slave=100, ny_slave=100, include_subpixel=True):
        values = {
            "rshift": 0,
            "ashift": 4,
            "stretch_r": 0,
            "a_stretch_r": 0,
            "stretch_a": 0,
            "a_stretch_a": 0,
        }
        if include_subpixel:
            values["sub_int_r"] = 0
            values["sub_int_a"] = 0
        if shift_values:
            values.update(shift_values)
        self.slave_prm.write_text(
            "".join(f"{key} = {value}\n" for key, value in values.items()),
            encoding="utf-8",
        )
        self.size_prm.write_text(
            f"num_rng_bins = {nx_slave}\nnum_lines = {ny_slave}\n",
            encoding="utf-8",
        )

    def _run(self, corr_thresh=0.3, margin_px=4):
        command = [
            sys.executable,
            str(SCRIPT),
            "--grid",
            str(self.grid),
            "--slave-prm",
            str(self.slave_prm),
            "--slave-size-prm",
            str(self.size_prm),
            "--corr-thresh",
            str(corr_thresh),
            "--margin-px",
            str(margin_px),
            "--out-mask",
            str(self.out_mask),
            "--summary-json",
            str(self.summary_path),
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        summary = json.loads(self.summary_path.read_text(encoding="utf-8")) if self.summary_path.exists() else None
        if result.returncode == 0:
            info = subprocess.run(
                ["gmt", "grdinfo", "-C", str(self.out_mask)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.split()
            nx, ny = int(float(info[9])), int(float(info[10]))
            raw = subprocess.run(
                ["gmt", "grd2xyz", str(self.out_mask), "-ZTLf"],
                check=True,
                capture_output=True,
            ).stdout
            mask = np.frombuffer(raw, dtype=np.float32).reshape((ny, nx))
        else:
            mask = None
        return result, mask, summary

    def test_range_shift_masks_exact_right_columns(self):
        self._write_prms(
            shift_values={"rshift": 8, "stretch_r": 0.1},
            nx_slave=60,
            ny_slave=100,
        )
        result, mask, summary = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        valid = np.isfinite(mask)
        self.assertTrue(np.all(valid[:, :11]))
        self.assertTrue(np.all(~valid[:, 11:]))
        self.assertEqual(summary["overlap_invalid_widths_px"]["right"], 9)
        self.assertAlmostEqual(summary["fractions"]["masked_fraction_by_overlap"], 9 / 20)

    def test_azimuth_shift_masks_top_and_bottom_rows(self):
        self._write_prms(shift_values={"ashift": 55})
        result, mask, summary = self._run(margin_px=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(np.all(np.isfinite(mask[:-1, :])))
        self.assertTrue(np.all(~np.isfinite(mask[-1, :])))
        self.assertEqual(summary["overlap_invalid_widths_px"]["bottom"], 1)

        self._write_prms(shift_values={"ashift": -10})
        result, mask, summary = self._run(margin_px=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(np.all(~np.isfinite(mask[:2, :])))
        self.assertTrue(np.all(np.isfinite(mask[2:, :])))
        self.assertEqual(summary["overlap_invalid_widths_px"]["top"], 2)

    def test_coherence_threshold_masks_low_corr_pixels(self):
        self.corr[:, 5] = 0.2
        self._write_grid(self.corr)
        result, mask, summary = self._run(corr_thresh=0.3, margin_px=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(np.all(~np.isfinite(mask[:, 5])))
        self.assertTrue(np.all(np.isfinite(mask[:, :5])))
        self.assertAlmostEqual(summary["fractions"]["masked_fraction_by_corr"], 1 / 20)

    def test_zero_threshold_disables_coherence_mask(self):
        self.corr[:, 5] = 0.2
        self._write_grid(self.corr)
        result, mask, summary = self._run(corr_thresh=0.0, margin_px=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(np.all(np.isfinite(mask)))
        self.assertEqual(summary["fractions"]["masked_fraction_by_corr"], 0.0)

    def test_missing_subpixel_keys_default_to_zero_and_are_reported(self):
        self._write_prms(include_subpixel=False)
        result, mask, summary = self._run(margin_px=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(summary["subpixel_defaults"], ["sub_int_r", "sub_int_a"])
        self.assertTrue(np.all(np.isfinite(mask)))

    def test_missing_required_prm_key_exits_one_with_key_name(self):
        self.slave_prm.write_text(
            "ashift = 0\nstretch_r = 0\na_stretch_r = 0\n"
            "stretch_a = 0\na_stretch_a = 0\n",
            encoding="utf-8",
        )
        result, _, _ = self._run()
        self.assertEqual(result.returncode, 1)
        self.assertIn("rshift", result.stderr)


class RunDinsarMaskCorrOptionTests(unittest.TestCase):
    def _parse(self, args):
        command = (
            "source ./run_dinsar.sh; LOG_FILE=/dev/null; "
            f"parse_args {args} /tmp; printf '%s' \"$MASK_CORR\""
        )
        return subprocess.run(
            ["bash", "-c", command],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )

    def test_default_mask_threshold(self):
        result = subprocess.run(
            ["bash", "-c", "source ./run_dinsar.sh; printf '%s' \"$MASK_CORR\""],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "0.3")

    def test_valid_mask_thresholds_parse(self):
        for threshold in ("0", ".25", "0.3", "0.999", "1e-3"):
            with self.subTest(threshold=threshold):
                result = self._parse(f"--mask-corr {threshold}")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, threshold)

    def test_invalid_mask_thresholds_are_rejected(self):
        for threshold in ("-0.1", "1", "nan", "abc"):
            with self.subTest(threshold=threshold):
                result = self._parse(f"--mask-corr {threshold}")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("无效 --mask-corr", result.stdout + result.stderr)

    def test_help_mentions_mask_flag(self):
        result = subprocess.run(
            ["bash", str(RUN_DINSAR), "--help"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("--mask-corr T", result.stdout)


if __name__ == "__main__":
    unittest.main()
