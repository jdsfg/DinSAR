import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from osgeo import gdal


REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "gmtsar-master/gmtsar/csh"


@unittest.skipUnless(all(shutil.which(c) for c in ("bash", "csh", "gmt")), "需要 bash、csh、GMT")
class AtmosphericWorkflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="atm_workflow_")
        self.addCleanup(self.tmp.cleanup)
        self.project = Path(self.tmp.name)
        self.intf = self.project / "intf/20260904_20260912"
        for path in (self.intf, self.project / "topo", self.project / "debug", self.project / "weather", self.project / "SLC"):
            path.mkdir(parents=True, exist_ok=True)
        self.prm = self.intf / "20260904.PRM"
        self.prm.write_text("led_file = 20260904.LED\nrng_samp_rate = 60000000\nPRF = 2000\n"
                            "SC_vel = 7500\nradar_wavelength = 0.031\nSC_height = 600000\nSC_height = 600001\n")
        self.led = self.project / "SLC/20260904.LED"
        self.led.write_text("orbit fixture\n")
        self.log = self.project / "processing.log"
        self.stub = self.project / "correction.csh"
        self.gmt("grdmath", "-R100/140/200/240", "-I4", "X", "Y", "ADD", "=", "unwrap.grd")
        self.gmt("grdmath", "-R100/140/200/240", "-I4", "1500", "X", "ADD", "=", "dem_ra.grd")
        shutil.copyfile(self.intf / "dem_ra.grd", self.project / "topo/dem.grd")
        (self.intf / "trans.dat").write_bytes(b"fixture")
        for date in ("20260904", "20260912"):
            for ext in ("ztd", "ztd.rsc"):
                (self.project / "weather" / f"{date}.{ext}").write_bytes(b"fixture")
        self.original = (self.intf / "unwrap.grd").read_bytes()
        self.write_stub('gmt grdmath "$1" 1 SUB = "$3"\nexit $status\n')

    def gmt(self, *args):
        return subprocess.run(["gmt", *args], cwd=self.intf, check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)

    def write_stub(self, body):
        self.stub.write_text("#!/bin/csh -f\n" + body)

    def run_atm(self, mode="EMPIRICAL", script=None, ref=None, extra=""):
        env = dict(os.environ, TEST_PROJECT=str(self.project), TEST_SCRIPT=str(script or self.stub),
                   TEST_MODE=mode, PYTHONDONTWRITEBYTECODE="1")
        snippet = '''
source ./run_dinsar.sh
PROJECT_DIR="$TEST_PROJECT"
INTF_DIR="$PROJECT_DIR/intf/20260904_20260912"
LOG_FILE="$PROJECT_DIR/processing.log"
MASTER_DATE=20260904
SLAVE_DATE=20260912
ATM_MODE="$TEST_MODE"
GACOS_PATH="$PROJECT_DIR/weather"
find_atm_script() { printf '%s\\n' "$TEST_SCRIPT"; }
'''
        if ref:
            snippet += f"REF_PX={ref[0]}; REF_PY={ref[1]}\n"
        snippet += extra + "\nrun_atm_correction\n"
        return subprocess.run(["bash", "-c", snippet], cwd=REPO, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)

    def assert_failed_without_replacing_input(self, result):
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), self.original)
        self.assertFalse((self.intf / "unwrap_pre_atm.grd").exists())
        self.assertFalse((self.intf / "atm_correction.state").exists())
        self.assertEqual(list((self.project / "debug").glob(".atm-work.*")), [])

    def values(self, name):
        return gdal.Open(str(self.intf / name)).ReadAsArray()

    def test_none_does_not_touch_files_or_run_scripts(self):
        self.write_stub("exit 99\n")
        result = self.run_atm("NONE")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), self.original)
        self.assertFalse((self.intf / "unwrap_atm.grd").exists())

    def test_valid_result_and_matching_checkpoint(self):
        original_values = self.values("unwrap.grd")
        first = self.run_atm()
        self.assertEqual(first.returncode, 0, first.stdout + self.log.read_text())
        np.testing.assert_allclose(self.values("unwrap.grd"), original_values - 1)
        self.assertEqual((self.intf / "unwrap_pre_atm.grd").read_bytes(), self.original)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        second = self.run_atm()
        self.assertEqual(second.returncode, 0, second.stdout)
        self.assertIn("Checkpoint", second.stdout)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)
        self.assertEqual(list((self.project / "debug").glob(".atm-work.*")), [])

    def test_nonzero_exit_rejects_new_and_stale_outputs(self):
        stale = b"previous failed correction"
        (self.intf / "unwrap_atm.grd").write_bytes(stale)
        self.write_stub('gmt grdmath "$1" 1 SUB = "$3"\nexit 8\n')
        self.assert_failed_without_replacing_input(self.run_atm())
        self.assertEqual((self.intf / "unwrap_atm.grd").read_bytes(), stale)

    def test_zero_exit_without_output_does_not_accept_stale_file(self):
        shutil.copyfile(self.intf / "unwrap.grd", self.intf / "unwrap_atm.grd")
        self.write_stub("exit 0\n")
        self.assert_failed_without_replacing_input(self.run_atm())

    def test_wrong_dimensions_rejected(self):
        self.write_stub('gmt grdmath -R0/3/0/3 -I1 7 = "$3"\n')
        self.assert_failed_without_replacing_input(self.run_atm())

    def test_shifted_grid_rejected(self):
        self.write_stub('gmt grdmath -R101/141/200/240 -I4 7 = "$3"\n')
        self.assert_failed_without_replacing_input(self.run_atm())

    def test_invalid_footprint_rejected(self):
        self.write_stub('gmt grdmath "$1" 0 DIV = "$3"\n')
        self.assert_failed_without_replacing_input(self.run_atm())

    def test_downstream_checkpoint_blocks_new_correction(self):
        shutil.copyfile(self.intf / "unwrap.grd", self.intf / "unwrap_detrended.grd")
        self.assert_failed_without_replacing_input(self.run_atm())

    def test_unverifiable_old_checkpoint_is_preserved(self):
        shutil.copyfile(self.intf / "unwrap.grd", self.intf / "unwrap_pre_atm.grd")
        result = self.run_atm()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.intf / "unwrap_pre_atm.grd").read_bytes(), self.original)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), self.original)

    def test_changed_mode_and_none_do_not_reuse_correction(self):
        self.assertEqual(self.run_atm().returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        for mode in ("ST_FILTER", "NONE"):
            with self.subTest(mode=mode):
                result = self.run_atm(mode)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_modified_corrected_phase_rejected(self):
        self.assertEqual(self.run_atm().returncode, 0)
        self.gmt("grdmath", "unwrap.grd", "2", "ADD", "=", "unwrap.grd")
        changed = (self.intf / "unwrap.grd").read_bytes()
        self.assertNotEqual(self.run_atm().returncode, 0)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), changed)

    def test_modified_model_input_rejected(self):
        self.assertEqual(self.run_atm().returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        self.gmt("grdmath", "dem_ra.grd", "10", "ADD", "=", "dem_ra.grd")
        self.assertNotEqual(self.run_atm().returncode, 0)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_modified_script_rejected(self):
        self.assertEqual(self.run_atm().returncode, 0)
        self.write_stub('gmt grdmath "$1" 2 SUB = "$3"\n')
        self.assertNotEqual(self.run_atm().returncode, 0)

    def test_incomplete_checkpoint_with_missing_phase_is_rejected(self):
        self.assertEqual(self.run_atm().returncode, 0)
        (self.intf / "unwrap.grd").unlink()
        self.assertNotEqual(self.run_atm().returncode, 0)
        self.assertEqual((self.intf / "unwrap_pre_atm.grd").read_bytes(), self.original)

    def test_checkpoint_with_missing_dem_is_rejected(self):
        self.assertEqual(self.run_atm().returncode, 0)
        (self.project / "topo/dem.grd").unlink()
        self.assertNotEqual(self.run_atm().returncode, 0)

    def test_orphan_state_rejects_correction_and_none(self):
        self.assertEqual(self.run_atm().returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        state = (self.intf / "atm_correction.state").read_bytes()
        (self.intf / "unwrap_pre_atm.grd").unlink()
        for mode in ("EMPIRICAL", "ST_FILTER", "GACOS", "NONE"):
            with self.subTest(mode=mode):
                self.assertNotEqual(self.run_atm(mode).returncode, 0)
                self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)
                self.assertEqual((self.intf / "atm_correction.state").read_bytes(), state)
                self.assertFalse((self.intf / "unwrap_pre_atm.grd").exists())

    def test_dangling_checkpoint_is_not_treated_as_new_run(self):
        (self.intf / "atm_correction.state").symlink_to(self.intf / "missing_state")
        for mode in ("EMPIRICAL", "NONE"):
            self.assertNotEqual(self.run_atm(mode).returncode, 0)
            self.assertEqual((self.intf / "unwrap.grd").read_bytes(), self.original)
            self.assertTrue((self.intf / "atm_correction.state").is_symlink())
            self.assertFalse((self.intf / "unwrap_pre_atm.grd").exists())

    def test_empirical_checkpoint_ignores_unconsumed_prm_fields(self):
        self.assertEqual(self.run_atm().returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        self.prm.write_text(self.prm.read_text().replace("SC_height = 600000\n", "") + "rshift = 0\n")
        result = self.run_atm()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Checkpoint", result.stdout)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_spatial_checkpoint_uses_effective_scale_not_raw_prm(self):
        script = SCRIPTS / "ST_Filter.csh"
        self.assertEqual(self.run_atm("ST_FILTER", script=script).returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        self.prm.write_text(self.prm.read_text().replace("SC_height = 600000\n", "")
                            .replace("60000000", "6.000000000e7") + "rshift = 0\n")
        result = self.run_atm("ST_FILTER", script=script)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Checkpoint", result.stdout)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_spatial_checkpoint_rejects_changed_pixel_scale(self):
        script = SCRIPTS / "ST_Filter.csh"
        self.assertEqual(self.run_atm("ST_FILTER", script=script).returncode, 0)
        corrected = (self.intf / "unwrap.grd").read_bytes()
        self.prm.write_text(self.prm.read_text().replace("PRF = 2000", "PRF = 4000"))
        self.assertNotEqual(self.run_atm("ST_FILTER", script=script).returncode, 0)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_spatial_checkpoint_rejects_changed_fallback_pixel_scale(self):
        self.prm.unlink()
        script = SCRIPTS / "ST_Filter.csh"
        self.assertEqual(self.run_atm("ST_FILTER", script=script, extra="export ST_PIXEL_M=3").returncode, 0)
        self.assertNotEqual(self.run_atm("ST_FILTER", script=script, extra="export ST_PIXEL_M=4").returncode, 0)

    def assert_gacos_reads_resolved_orbit(self):
        original_prm = self.prm.read_bytes()
        self.write_stub('set orbit = `awk \'$1 == "led_file" {print $3; exit}\' 20260904.PRM`\n'
                        'if ( ! -f "$orbit" ) exit 7\n'
                        'if ( "`cat $orbit`" != "orbit fixture" ) exit 8\n'
                        'gmt grdmath unwrap.grd 2 SUB = unwrap_gacos_corrected.grd\n')
        result = self.run_atm("GACOS")
        self.assertEqual(result.returncode, 0, result.stdout + self.log.read_text())
        self.assertEqual(self.prm.read_bytes(), original_prm)
        self.assertEqual(self.run_atm("GACOS").returncode, 0)

    def test_gacos_resolves_orbit_from_slc_without_local_link(self):
        self.assertFalse((self.intf / self.led.name).exists())
        self.assert_gacos_reads_resolved_orbit()

    def test_gacos_resolves_orbit_from_raw(self):
        (self.project / "raw").mkdir()
        self.led.rename(self.project / "raw" / self.led.name)
        self.assert_gacos_reads_resolved_orbit()

    def test_gacos_preserves_relative_orbit_reference(self):
        self.prm.write_text(self.prm.read_text().replace("led_file = 20260904.LED", "led_file = ../../SLC/20260904.LED"))
        self.assert_gacos_reads_resolved_orbit()

    def test_gacos_preserves_absolute_orbit_reference(self):
        self.prm.write_text(self.prm.read_text().replace("led_file = 20260904.LED", f"led_file = {self.led}"))
        self.assert_gacos_reads_resolved_orbit()

    def test_gacos_checkpoint_rejects_changed_orbit(self):
        self.assert_gacos_reads_resolved_orbit()
        corrected = (self.intf / "unwrap.grd").read_bytes()
        self.led.write_text("changed orbit\n")
        self.assertNotEqual(self.run_atm("GACOS").returncode, 0)
        self.assertEqual((self.intf / "unwrap.grd").read_bytes(), corrected)

    def test_gacos_missing_orbit_rejected_before_script(self):
        self.led.unlink()
        self.write_stub('gmt grdmath unwrap.grd 2 SUB = unwrap_gacos_corrected.grd\n')
        self.assert_failed_without_replacing_input(self.run_atm("GACOS"))

    def test_gacos_reference_uses_grid_coordinates(self):
        self.write_stub('echo "REFERENCE $3 $4"\n'
                        'gmt grdmath unwrap.grd 2 SUB = unwrap_gacos_corrected.grd\n')
        result = self.run_atm("GACOS")
        self.assertEqual(result.returncode, 0, result.stdout + self.log.read_text())
        self.assertIn("REFERENCE 120 220", self.log.read_text())
        self.assertTrue((self.intf / "unwrap_gacos_corrected.grd").exists())
        self.assertEqual(self.run_atm("GACOS").returncode, 0)
        self.assertNotEqual(self.run_atm("GACOS", ref=(124, 220)).returncode, 0)

    def test_changed_gacos_data_rejected(self):
        self.write_stub('gmt grdmath unwrap.grd 2 SUB = unwrap_gacos_corrected.grd\n')
        self.assertEqual(self.run_atm("GACOS").returncode, 0)
        (self.project / "weather/20260912.ztd").write_bytes(b"new weather")
        self.assertNotEqual(self.run_atm("GACOS").returncode, 0)

    def test_gacos_failure_does_not_accept_old_product(self):
        shutil.copyfile(self.intf / "unwrap.grd", self.intf / "unwrap_gacos_corrected.grd")
        self.write_stub("exit 9\n")
        self.assert_failed_without_replacing_input(self.run_atm("GACOS"))

    def test_real_empirical_script_preserves_grid_and_known_height_model(self):
        result = self.run_atm(script=SCRIPTS / "Re_atm.csh")
        self.assertEqual(result.returncode, 0, result.stdout + self.log.read_text())
        corrected = self.values("unwrap.grd")
        np.testing.assert_allclose(corrected[:, 0], corrected[:, -1], atol=1e-4)
        self.assertGreater(float(np.ptp(corrected[:, 0])), 30)

    def test_real_spatial_filter_checkpoint_is_repeatable(self):
        result = self.run_atm("ST_FILTER", script=SCRIPTS / "ST_Filter.csh")
        self.assertEqual(result.returncode, 0, result.stdout + self.log.read_text())
        again = self.run_atm("ST_FILTER", script=SCRIPTS / "ST_Filter.csh")
        self.assertEqual(again.returncode, 0, again.stdout)

    def test_raytrace_mode_valid_result_and_matching_checkpoint(self):
        self.write_stub('gmt grdmath "$1" 1 SUB = "$2"\nexit 0\n')
        weather_dir = self.project / "weather_nc"
        weather_dir.mkdir(parents=True, exist_ok=True)
        (weather_dir / "master.nc").write_bytes(b"mock master nc")
        (weather_dir / "slave.nc").write_bytes(b"mock slave nc")

        extra = f'WEATHER_DIR="{weather_dir}"'
        first = self.run_atm("RAYTRACE", extra=extra)
        self.assertEqual(first.returncode, 0, first.stdout + self.log.read_text())
        np.testing.assert_allclose(self.values("unwrap.grd"), self.values("unwrap_pre_atm.grd") - 1)

        second = self.run_atm("RAYTRACE", extra=extra)
        self.assertEqual(second.returncode, 0, second.stdout)
        self.assertIn("Checkpoint", second.stdout)

    def test_raytrace_missing_dem_rejected(self):
        self.write_stub('gmt grdmath "$1" 1 SUB = "$2"\nexit 0\n')
        (self.project / "topo/dem.grd").unlink()
        result = self.run_atm("RAYTRACE")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("RAYTRACE 模式需要有效的 topo/dem.grd", result.stdout + self.log.read_text())

    def test_raytrace_missing_trans_dat_rejected(self):
        self.write_stub('gmt grdmath "$1" 1 SUB = "$2"\nexit 0\n')
        (self.intf / "trans.dat").unlink()
        result = self.run_atm("RAYTRACE")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("RAYTRACE 模式需要有效的 trans.dat", result.stdout + self.log.read_text())

    def test_raytrace_changed_weather_dir_invalidates_checkpoint(self):
        self.write_stub('gmt grdmath "$1" 1 SUB = "$2"\nexit 0\n')
        weather_dir = self.project / "weather_nc"
        weather_dir.mkdir(parents=True, exist_ok=True)
        (weather_dir / "master.nc").write_bytes(b"mock master nc")
        first = self.run_atm("RAYTRACE", extra=f'WEATHER_DIR="{weather_dir}"')
        self.assertEqual(first.returncode, 0, first.stdout + self.log.read_text())

        # Modify weather file contents
        (weather_dir / "master.nc").write_bytes(b"updated master nc")
        second = self.run_atm("RAYTRACE", extra=f'WEATHER_DIR="{weather_dir}"')
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("大气校正检查点缺少来源记录或输入/模式/结果已变", second.stdout + self.log.read_text())


if __name__ == "__main__":
    unittest.main()
