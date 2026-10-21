from datetime import datetime, timezone
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from osgeo import gdal
from scipy.io import netcdf_file


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "gmtsar-master/gmtsar/csh"
sys.path.insert(0, str(SCRIPTS))
import era5_weather
import raytrace_atm_workflow as workflow


class Era5Weather(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build_dir = tempfile.TemporaryDirectory(prefix="era5 raytrace build ")
        cls.addClassCleanup(cls.build_dir.cleanup)
        cls.binary = Path(cls.build_dir.name) / "atm_raytrace"
        flags = shlex.split(subprocess.check_output(["nc-config", "--cflags"], text=True))
        libs = shlex.split(subprocess.check_output(["nc-config", "--libs"], text=True))
        source = ROOT / "gmtsar-master/gmtsar/atm_raytrace.cpp"
        subprocess.run(
            ["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
             *flags, str(source), "-o", str(cls.binary), *libs],
            check=True,
            timeout=90,
        )

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="era5 weather fixture ")
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.epoch0 = np.datetime64("2024-02-29T12:00:00", "ns")
        self.epoch1 = np.datetime64("2024-02-29T13:00:00", "ns")
        self.raw = self.synthetic_raw()
        self.geoid = self.folder / "zero geoid.grd"
        subprocess.run(
            ["gmt", "grdmath", "-R107/110/39/42", "-I0.25", "0", "=", str(self.geoid)],
            check=True,
            capture_output=True,
            timeout=30,
        )

    def synthetic_raw(self):
        times = np.array([self.epoch0, self.epoch1])
        levels = np.array([1000, 900, 800, 700, 500, 300, 100, 1], dtype=float)
        latitude = np.array([41.0, 40.5, 40.0, 39.5, 39.0])
        longitude = np.array([108.0, 108.25, 108.5, 108.75, 109.0])
        heights = np.array([0, 1000, 2000, 3000, 5000, 9000, 16000, 32000], dtype=float)
        shape = (len(times), len(levels), len(latitude), len(longitude))
        temperature = np.empty(shape)
        humidity = np.empty(shape)
        geopotential = np.empty(shape)
        for t in range(2):
            temperature[t] = (285 + 1.5 * t - heights[:, None, None] * 0.004
                              + latitude[None, :, None] * 0.01)
            humidity[t] = np.maximum(
                1e-5,
                np.array([0.010, 0.008, 0.006, 0.004, 0.002, 0.0005, 0.00005, 0.00001])[:, None, None]
                * (1 + 0.1 * t),
            )
            geopotential[t] = (heights[:, None, None] + 0.2 * t) * 9.80665
        surface = np.zeros((2, len(latitude), len(longitude)))
        return {
            "time": times,
            "level": levels,
            "latitude": latitude,
            "longitude": longitude,
            "geopotential": geopotential,
            "temperature": temperature,
            "specific_humidity": humidity,
            "surface_pressure": np.full((2, len(latitude), len(longitude)), 101000.0)
                                + np.arange(2)[:, None, None] * 100,
            "geopotential_at_surface": surface,
            "2m_temperature": np.full((2, len(latitude), len(longitude)), 285.0)
                              + np.arange(2)[:, None, None] * 1.5,
            "2m_dewpoint_temperature": np.full((2, len(latitude), len(longitude)), 280.0)
                                       + np.arange(2)[:, None, None] * 0.5,
        }

    def build(self, raw, epoch, name):
        output = self.folder / name
        era5_weather.build_cube(raw, epoch, (108.0, 109.0, 39.5, 40.5), self.geoid, output)
        return output

    def read_variable(self, path, variable):
        with netcdf_file(str(path), "r", mmap=False) as nc:
            return np.array(nc.variables[variable][:]), getattr(nc, "epoch_utc", None)

    def test_cube_is_accepted_by_atm_raytrace(self):
        raw_path = self.folder / "synthetic_era5_raw.npz"
        np.savez_compressed(raw_path, **self.raw)
        with np.load(raw_path, allow_pickle=False) as archive:
            raw = {key: archive[key] for key in archive.files}
        cube = self.build(raw, datetime(2024, 2, 29, 12, 30, tzinfo=timezone.utc), "weather.nc")
        rays = self.folder / "rays.txt"
        delays = self.folder / "delays.txt"
        rays.write_text("108.5 40 0 0 0 1\n")
        result = subprocess.run(
            [str(self.binary), str(cube), str(rays), str(delays), "100"],
            capture_output=True,
            text=True,
            timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        row = [float(value) for value in delays.read_text().splitlines()[-1].split()]
        self.assertGreater(row[6], 0.1)
        self.assertLess(row[6], 10.0)

    def test_half_hour_interpolates_cube_temperature_linearly(self):
        raw0 = {name: values[:1] if values.ndim >= 3 and values.shape[0] == 2 else values
                for name, values in self.raw.items()}
        raw1 = dict(raw0)
        for name, values in self.raw.items():
            if values.ndim >= 3 and values.shape[0] == 2:
                raw1[name] = values[1:2]
        raw0["time"] = self.raw["time"][:1]
        raw1["time"] = self.raw["time"][1:]
        cube0 = self.build(raw0, datetime(2024, 2, 29, 12, tzinfo=timezone.utc), "t0.nc")
        cube1 = self.build(raw1, datetime(2024, 2, 29, 13, tzinfo=timezone.utc), "t1.nc")
        middle = self.build(
            self.raw, datetime(2024, 2, 29, 12, 30, tzinfo=timezone.utc), "mid.nc"
        )
        raw_static_surface = dict(self.raw)
        raw_static_surface["geopotential_at_surface"] = self.raw["geopotential_at_surface"][0]
        static_surface = self.build(
            raw_static_surface, datetime(2024, 2, 29, 12, 30, tzinfo=timezone.utc), "static.nc"
        )
        raw_seconds = dict(self.raw)
        raw_seconds["time"] = self.raw["time"].astype("datetime64[s]").astype(np.int64)
        integer_times = self.build(
            raw_seconds, datetime(2024, 2, 29, 12, 30, tzinfo=timezone.utc), "seconds.nc"
        )
        t0, _ = self.read_variable(cube0, "temperature")
        t1, _ = self.read_variable(cube1, "temperature")
        tm, _ = self.read_variable(middle, "temperature")
        self.assertAlmostEqual(tm[40, 1, 1], (t0[40, 1, 1] + t1[40, 1, 1]) / 2, delta=1e-6)
        ts, _ = self.read_variable(static_surface, "temperature")
        ti, _ = self.read_variable(integer_times, "temperature")
        np.testing.assert_allclose(ts, tm, atol=1e-6)
        np.testing.assert_allclose(ti, tm, atol=1e-6)

    def test_epoch_and_bbox_must_be_covered_by_raw(self):
        with self.assertRaisesRegex(ValueError, "do not bracket"):
            self.build(self.raw, datetime(2024, 2, 29, 14, tzinfo=timezone.utc), "late.nc")
        with self.assertRaisesRegex(ValueError, "outside raw ERA5 extent"):
            era5_weather.build_cube(
                self.raw, datetime(2024, 2, 29, 12, tzinfo=timezone.utc),
                (108.0, 110.0, 39.5, 40.5), self.geoid, self.folder / "outside.nc",
            )

    def test_prm_clock_start_parses_fractional_day_as_utc(self):
        prm = self.folder / "clock.PRM"
        prm.write_text("SC_clock_start = 2030005.1831530258\n")
        parsed = workflow.prm_epoch_utc(prm)
        expected = datetime(2030, 1, 5, 4, 23, 44, 400000, tzinfo=timezone.utc)
        self.assertLessEqual(abs((parsed - expected).total_seconds()), 1.0)

    def write_epoch_cube(self, path, epoch):
        with netcdf_file(str(path), "w") as nc:
            nc.epoch_utc = epoch

    def test_weather_discovery_has_no_fallback_and_checks_epochs(self):
        project = self.folder / "project"
        stage = project / "stage"
        project.mkdir()
        stage.mkdir()
        self.assertEqual(workflow.locate_weather_cubes(stage, project, ""), (None, None))
        master = stage / "master_weather.nc"
        slave = stage / "slave_weather.nc"
        self.write_epoch_cube(master, "2024-02-29T12:00:00Z")
        self.write_epoch_cube(slave, "2024-02-29T13:00:00Z")
        self.assertTrue(workflow.validate_cube_epoch(
            master, datetime(2024, 2, 29, 12, 10, tzinfo=timezone.utc)
        ))
        with self.assertRaisesRegex(ValueError, "epoch mismatch"):
            workflow.validate_cube_epoch(
                slave, datetime(2024, 2, 29, 12, tzinfo=timezone.utc)
            )

    def write_grid(self, path, values, geotransform, nodata=-9999.0):
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(str(path), values.shape[1], values.shape[0], 1, gdal.GDT_Float32)
        ds.SetGeoTransform(geotransform)
        band = ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        band.WriteArray(values)
        ds = None

    def test_nearest_fill_radius_and_target_geotransform(self):
        geotransform = (100.0, 1.0, 0.0, 200.0, 0.0, -1.0)
        values = np.full((7, 7), -9999.0, dtype=np.float32)
        values[2:5, 2:5] = 7.0
        values[3, 3] = -9999.0
        temp_grid = self.folder / "projected.tif"
        source_grid = self.folder / "source.tif"
        output_grid = self.folder / "resampled.nc"
        self.write_grid(temp_grid, values, geotransform)
        self.write_grid(source_grid, np.ones((7, 7), dtype=np.float32), geotransform)
        result = workflow.fill_and_resample(temp_grid, source_grid, output_grid, radius=2)
        self.assertEqual(result[3, 3], 7.0)
        self.assertTrue(np.isnan(result[0, 0]))
        output = gdal.Open(str(output_grid))
        np.testing.assert_allclose(output.GetGeoTransform(), geotransform)

    def test_correction_subtracts_scaled_plane_removed_atmosphere_and_preserves_mask(self):
        rows, cols = np.mgrid[:12, :12]
        unw = (2 + 0.1 * cols - 0.2 * rows + 0.1 * np.sin(cols + rows)).astype(np.float32)
        atm = (5 + 0.4 * cols + 0.2 * rows + 0.3 * np.sin(cols * 0.7 + rows)).astype(np.float32)
        valid = np.ones(unw.shape, dtype=bool)
        valid[0, 0] = False
        unw[0, 0] = -9999
        fit_mask = valid.copy()
        result, stratified = workflow.correct_phase(unw, atm, valid, fit_mask, scale=0.75)
        raw_corrected = unw - 0.75 * stratified
        reference = np.median(raw_corrected[fit_mask]) - np.median(unw[fit_mask])
        expected = unw.copy()
        expected[valid] = (raw_corrected - reference)[valid]
        np.testing.assert_allclose(result, expected, atol=1e-6)
        self.assertEqual(result[0, 0], -9999)
        unchanged, _ = workflow.correct_phase(unw, atm, valid, fit_mask, scale=0)
        np.testing.assert_array_equal(unchanged, unw)


if __name__ == "__main__":
    unittest.main()
