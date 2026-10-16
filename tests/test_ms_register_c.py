#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parity checks for the native ms_register core; numerical tests need the built binary."""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_ms_register as base  # noqa: E402

ms = base.ms
REPO = os.path.dirname(HERE)
LOCAL_BIN = os.path.join(REPO, "gmtsar-master", "gmtsar", "ms_register")
C_BIN = os.environ.get("MS_REGISTER_BIN") or LOCAL_BIN
HAVE_C = os.path.isfile(C_BIN) and os.access(C_BIN, os.X_OK)


def _mask_numbers(text):
    return re.sub(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", "#", text)


class EngineSelection(unittest.TestCase):
    def test_configured_binary_has_priority(self):
        with tempfile.TemporaryDirectory() as d:
            binary = os.path.join(d, "ms_register")
            with open(binary, "w") as f:
                f.write("#!/bin/sh\nexit 0\n")
            os.chmod(binary, 0o755)
            with mock.patch.dict(os.environ, {"MS_REGISTER_BIN": binary}):
                self.assertEqual(ms._find_c_engine(), binary)

    def test_missing_binary_explains_python_fallback(self):
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(ms.shutil, "which", return_value=None), \
                mock.patch.object(ms.os.path, "isfile", return_value=False):
            with self.assertRaisesRegex(ms.InputError, "--engine python"):
                ms._find_c_engine()

    def test_register_parser_exposes_engine_choice(self):
        with self.assertRaises(SystemExit) as caught:
            ms.main(["register", "--help"])
        self.assertEqual(caught.exception.code, 0)


@unittest.skipUnless(HAVE_C, "build gmtsar-master/gmtsar/ms_register to run C parity tests")
class NativeCore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ms_register_c_")
        cls.good = os.path.join(cls.tmp, "good")
        cls.bad = os.path.join(cls.tmp, "unrelated")
        os.makedirs(cls.good)
        os.makedirs(cls.bad)
        base.make_pair(cls.good, 0.9)
        base.make_pair(cls.bad, 0.9, unrelated=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _estimate(self, engine, directory, prefix, output=None):
        report = os.path.join(directory, prefix + ".json")
        output = output or os.path.join(directory, prefix + ".PRM")
        log = os.path.join(directory, prefix + ".log")
        points = os.path.join(directory, prefix + ".dat")
        common = [os.path.join(directory, "m.PRM"), os.path.join(directory, "s.PRM"),
                  "--output", output, "--report", report, "--log", log,
                  "--points", points, "--write-on-fail", "--threads", "1"]
        if engine == "python":
            rc = ms.main(["estimate"] + common + ["-q"])
        else:
            proc = subprocess.run([C_BIN, "estimate"] + common + ["--quiet"],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  universal_newlines=True)
            rc = proc.returncode
            self.assertIn(rc, (ms.EXIT_OK, ms.EXIT_GATE), proc.stdout)
        with open(report) as f:
            return rc, json.load(f)

    def test_seed_zero_bounded_draws_match_numpy(self):
        proc = subprocess.run([C_BIN, "--rng-dump", "0", "10000"],
                              check=True, stdout=subprocess.PIPE, universal_newlines=True)
        actual = np.fromstring(proc.stdout, sep=" ", dtype=np.uint32)
        expected = np.random.default_rng(0).integers(0, 100000, size=10000, dtype=np.uint32)
        np.testing.assert_array_equal(actual, expected)

    def test_ransac_triplet_stream_matches_numpy_across_scales(self):
        bounds = (137, 73)
        proc = subprocess.run([C_BIN, "--rng-pattern", "0", str(bounds[0]), str(bounds[1])],
                              check=True, stdout=subprocess.PIPE, universal_newlines=True)
        actual = np.fromstring(proc.stdout, sep=" ", dtype=np.uint32)
        rng = np.random.default_rng(0)
        expected = np.concatenate([
            rng.integers(0, n, size=(1500, 3)).reshape(-1) for n in bounds
        ])
        np.testing.assert_array_equal(actual, expected)

    def test_prm_roundtrip_preserves_unknown_lines_and_negative_fraction(self):
        source = os.path.join(self.tmp, "roundtrip_in.PRM")
        target = os.path.join(self.tmp, "roundtrip_out.PRM")
        with open(source, "w") as f:
            f.write("num_rng_bins = 64\nopaque_key = keep-me\nrshift = 99\n"
                    "ashift = 88\nSLC_file = sample.SLC\n")
        subprocess.run([C_BIN, "--prm-roundtrip", source, target,
                        "-2.25", "0.125", "-0.00001",
                        "-0.75", "0.00002", "0.00003"], check=True)
        with open(target) as f:
            text = f.read()
        self.assertIn("opaque_key = keep-me\n", text)
        self.assertIn("SLC_file = sample.SLC\n", text)
        self.assertIn("rshift = -3\n", text)
        self.assertIn("sub_int_r = 0.750000\n", text)
        self.assertIn("ashift = -1\n", text)
        self.assertIn("sub_int_a = 0.250000\n", text)

    def test_shifted_pair_stage_offsets_counts_and_gate_match(self):
        python_rc, py = self._estimate("python", self.good, "python_good")
        c_rc, c = self._estimate("c", self.good, "c_good")
        self.assertEqual(c_rc, python_rc)
        self.assertEqual(c["gate"]["passed"], py["gate"]["passed"])
        for k in ("dr", "da"):
            self.assertLess(abs(c["global"][k] - py["global"][k]), 0.01, k)
        self.assertEqual(len(c["stages"]), len(py["stages"]))
        for cs, ps in zip(c["stages"], py["stages"]):
            self.assertEqual(cs["name"], ps["name"])
            self.assertEqual(cs["n_valid"], ps["n_valid"], cs["name"])
            self.assertEqual(cs["n_inliers"], ps["n_inliers"], cs["name"])
            for key in ("center_dr", "center_da"):
                self.assertLess(abs(cs[key] - ps[key]), 0.01, "%s %s" % (cs["name"], key))
        for key in ("range", "azimuth"):
            self.assertLess(abs(c["center_offset"][key] - py["center_offset"][key]), 0.01, key)
        self.assertEqual(c["check"]["n_valid"], py["check"]["n_valid"])
        for key in ("rms_r", "rms_a"):
            self.assertLessEqual(abs(c["check"][key] - py["check"][key]), 0.005, key)

    def test_outputs_match_python_structure_after_masking_numbers(self):
        shared_output = os.path.join(self.good, "structure_estimate.PRM")
        python_rc, py = self._estimate("python", self.good, "python_structure", shared_output)
        c_rc, c = self._estimate("c", self.good, "c_structure", shared_output)
        self.assertEqual(c_rc, python_rc)
        self.assertEqual(
            _mask_numbers(json.dumps(py, ensure_ascii=False, separators=(",", ":"))),
            _mask_numbers(json.dumps(c, ensure_ascii=False, separators=(",", ":"))),
            "report JSON structure or key order differs",
        )
        for extension in ("log", "dat"):
            with open(os.path.join(self.good, "python_structure." + extension)) as f:
                python_text = f.read()
            with open(os.path.join(self.good, "c_structure." + extension)) as f:
                c_text = f.read()
            self.assertEqual(
                _mask_numbers(python_text), _mask_numbers(c_text),
                "%s structure differs" % extension,
            )

    def test_unrelated_pair_failure_status_matches(self):
        python_rc, py = self._estimate("python", self.bad, "python_bad")
        c_rc, c = self._estimate("c", self.bad, "c_bad")
        self.assertEqual(c_rc, python_rc)
        self.assertEqual(c["gate"]["passed"], py["gate"]["passed"])
        self.assertEqual(c["status"], py["status"])


if __name__ == "__main__":
    unittest.main()
