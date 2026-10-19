#!/usr/bin/env python3
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
GENERATOR = REPO / "snaphu_conf_gen.py"
RUN_DINSAR = REPO / "run_dinsar.sh"


class SnaphuConfGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="snaphu_conf_gen_")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base = self.root / "base.conf"
        self.prm = self.root / "master.PRM"
        self.out = self.root / "auto.conf"
        self.base.write_text(
            "INFILEFORMAT FLOAT_DATA\n"
            "LAMBDA 0.9\n"
            "NEARRANGE 1\n"
            "EARTHRADIUS 2\n"
            "NCORRLOOKS 3\n"
            "NCORRLOOKSRANGE 4\n"
            "NCORRLOOKSAZ 5\n"
            "KEEP THIS\n",
            encoding="utf-8",
        )
        self.prm.write_text(
            "radar_wavelength = 0.0312284\n"
            "near_range = 640074.4387\n"
            "earth_radius = 6368727.157867\n"
            "SC_height = 530982.849128\n"
            "SC_vel = 7352.266526\n"
            "rng_samp_rate = 200000000\n"
            "PRF = 4240.882103\n",
            encoding="utf-8",
        )

    def run_generator(self, prm=None):
        return subprocess.run(
            [
                sys.executable,
                str(GENERATOR),
                "--base",
                str(self.base),
                "--prm",
                str(prm or self.prm),
                "--downsample",
                "4",
                "--corr-filter-width",
                "20",
                "--out",
                str(self.out),
            ],
            cwd=REPO,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )

    @staticmethod
    def active_config_values(text):
        values = {}
        for line in text.splitlines():
            match = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s+([^\s#]+)", line)
            if match and not line.lstrip().startswith("#"):
                values[match.group(1).upper()] = match.group(2)
        return values

    def test_derives_geometry_and_corr_looks_and_replaces_or_appends(self):
        result = self.run_generator()
        self.assertEqual(result.returncode, 0, result.stdout)
        config = self.out.read_text(encoding="utf-8")
        values = self.active_config_values(config)

        c = 299792458.0
        earth = 6368727.157867
        height = 530982.849128
        velocity = 7352.266526
        prf = 4240.882103
        range_res = c / (2.0 * 200000000.0)
        az_res = velocity * earth / (earth + height) / prf
        self.assertAlmostEqual(float(values["DR"]), range_res * 4, places=10)
        self.assertAlmostEqual(float(values["DA"]), az_res * 4, places=10)
        self.assertAlmostEqual(float(values["NCORRLOOKS"]), 4 * math.pi * (20 / 6) ** 2, places=9)
        self.assertAlmostEqual(float(values["LAMBDA"]), 0.0312284, places=10)
        self.assertAlmostEqual(float(values["NEARRANGE"]), 640074.4387, places=7)
        self.assertAlmostEqual(float(values["ORBITRADIUS"]), earth + height, delta=1e-4)
        self.assertEqual(values["NLOOKSRANGE"], "4")
        self.assertEqual(values["NLOOKSAZ"], "4")
        self.assertEqual(values["NLOOKSOTHER"], "1")
        self.assertEqual(values["RANGERES"], format(range_res, ".12g"))
        self.assertEqual(values["AZRES"], format(az_res, ".12g"))
        self.assertEqual(values["KEEP"], "THIS")
        self.assertEqual(config.count("LAMBDA\t"), 1)
        self.assertIn("# NCORRLOOKSRANGE 4", config)
        self.assertIn("# NCORRLOOKSAZ 5", config)
        self.assertNotIn("\nNCORRLOOKSRANGE ", config)
        self.assertNotIn("\nNCORRLOOKSAZ ", config)
        self.assertIn("# Source PRM: %s" % self.prm, config)

    def test_missing_prm_key_fails_clearly(self):
        incomplete = self.root / "incomplete.PRM"
        incomplete.write_text("near_range = 10\n", encoding="utf-8")
        result = self.run_generator(incomplete)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing required PRM key", result.stdout)
        self.assertIn("radar_wavelength", result.stdout)
        self.assertFalse(self.out.exists())


class RunDinsarSnaphuOptionTests(unittest.TestCase):
    def run_parse(self, args):
        return subprocess.run(
            ["bash", "-c", 'source ./run_dinsar.sh; LOG_FILE=/dev/null; parse_args %s; echo "mode=$SNAPHU_CONF_MODE path=$SNAPHU_CONF_PATH overlap=$SNAPHU_OVERLAP user=$SNAPHU_OVERLAP_USER_SET"; apply_tier_defaults; echo "after_tier=$SNAPHU_OVERLAP"' % args],
            cwd=REPO,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )

    def test_conf_and_overlap_modes_parse_and_override_mode_default(self):
        with tempfile.TemporaryDirectory(prefix="snaphu_cli_") as directory:
            conf = Path(directory) / "custom.conf"
            conf.write_text("LAMBDA 0.03\n", encoding="utf-8")
            r = self.run_parse('--mode quality --snaphu-conf auto --snaphu-overlap 321 "%s"' % directory)
            self.assertEqual(r.returncode, 0, r.stdout)
            self.assertIn("mode=auto path= overlap=321 user=true", r.stdout)
            self.assertIn("after_tier=321", r.stdout)

            r = self.run_parse('--snaphu-conf base "%s"' % directory)
            self.assertEqual(r.returncode, 0, r.stdout)
            self.assertIn("mode=base path= overlap=200 user=false", r.stdout)

            r = self.run_parse('--snaphu-conf "%s" --snaphu-overlap 77 "%s"' % (conf, directory))
            self.assertEqual(r.returncode, 0, r.stdout)
            self.assertIn("mode=path path=%s overlap=77 user=true" % conf, r.stdout)

    def test_invalid_overlap_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="snaphu_cli_") as directory:
            for value in ("0", "-2", "abc"):
                with self.subTest(value=value):
                    r = self.run_parse('--snaphu-overlap %s "%s"' % (value, directory))
                    self.assertNotEqual(r.returncode, 0, r.stdout)
                    self.assertIn("无效 --snaphu-overlap", r.stdout)
                    self.assertNotIn("mode=", r.stdout)

    def test_help_documents_new_options(self):
        r = subprocess.run(
            ["bash", str(RUN_DINSAR), "--help"],
            cwd=REPO,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertIn("--snaphu-conf auto|base|<path>", r.stdout)
        self.assertIn("--snaphu-overlap N", r.stdout)


if __name__ == "__main__":
    unittest.main()
