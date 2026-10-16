#!/usr/bin/env python3

import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "gmtsar-master/gmtsar/atm_raytrace.cpp"
WRAPPER = ROOT / "gmtsar-master/gmtsar/csh/raytrace_atm.csh"
EPOCH = "2024-02-29T12:34:56Z"
AXES = ("height", "latitude", "longitude")
UNITS = {"height": "m", "latitude": "degrees_north", "longitude": "degrees_east",
         "pressure": "Pa", "temperature": "K", "specific_humidity": "kg kg-1"}
A = 6378137.0
F = 1.0 / 298.257223563


def refractivity(p, t, q):
    e_pa = q * p / (0.622 + 0.378 * q)
    dry_hpa, wet_hpa = (p - e_pa) / 100.0, e_pa / 100.0
    return 77.6 * dry_hpa / t + 70.4 * wet_hpa / t + 373900.0 * wet_hpa / t**2


class AtmRaytrace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = tempfile.TemporaryDirectory(prefix="atm build with spaces ")
        cls.addClassCleanup(cls.build.cleanup)
        cls.binary = Path(cls.build.name) / "atm_raytrace"
        flags = shlex.split(subprocess.check_output(["nc-config", "--cflags"], text=True))
        libs = shlex.split(subprocess.check_output(["nc-config", "--libs"], text=True))
        subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
                        *flags, str(SOURCE), "-o", str(cls.binary), *libs], check=True, timeout=90)

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(prefix="atm small cube ")
        self.addCleanup(self.work.cleanup)
        self.folder = Path(self.work.name)
        self.weather = self.folder / "weather cube.nc"
        self.rays = self.folder / "rays input.txt"
        self.output = self.folder / "delays output.txt"
        self.cube()

    def cube(self, heights=(0, 1000, 2000), latitudes=(-1, 0, 1), longitudes=(-1, 0, 1),
             pressure=100000.0, temperature=288.0, humidity=0.0, dims=None, units=None):
        coords = dict(zip(AXES, (heights, latitudes, longitudes)))
        with netcdf_file(self.weather, "w") as nc:
            nc.height_reference = "WGS84_ellipsoid"
            nc.epoch_utc = EPOCH
            for name, values in coords.items():
                nc.createDimension(name, len(values))
                var = nc.createVariable(name, "d", (name,))
                var.units = (units or {}).get(name, UNITS[name])
                var[:] = values
            for name, values in (("pressure", pressure), ("temperature", temperature),
                                 ("specific_humidity", humidity)):
                order = (dims or {}).get(name, AXES)
                var = nc.createVariable(name, "d", order)
                var.units = (units or {}).get(name, UNITS[name])
                shape = tuple(len(coords[d]) for d in order)
                var[:] = np.broadcast_to(values, shape)

    def invoke(self, rays="0 0 0 0 0 1\n", step=None, command=None, output=None, env=None):
        self.rays.write_text(rays)
        args = [str(self.binary)] if command is None else [str(p) for p in command]
        args += [str(self.weather), str(self.rays), str(output or self.output)]
        if step is not None:
            args.append(str(step))
        return subprocess.run(args, capture_output=True, text=True, timeout=20, env=env)

    def success(self, **kwargs):
        result = self.invoke(**kwargs)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = self.output.read_text().splitlines()
        return np.array([[float(v) for v in line.split()] for line in lines if not line.startswith("#")])

    def failure(self, message=None, **kwargs):
        result = self.invoke(**kwargs)
        self.assertNotEqual(result.returncode, 0)
        if message:
            self.assertIn(message, result.stderr)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.folder.glob("*.tmp.*")), [])
        return result

    def test_help_and_usage_contract(self):
        result = subprocess.run([str(self.binary), "--help"], text=True, capture_output=True, timeout=10)
        self.assertEqual(result.returncode, 0)
        for text in ("normalized", "NOT a raw WRF/ERA5", "WGS84_ellipsoid", "epoch_utc",
                     "no automatic main-pipeline", "NO ray bending", "128 MiB", "1000000"):
            self.assertIn(text, result.stdout)
        result = subprocess.run([str(self.binary)], capture_output=True, timeout=10)
        self.assertEqual(result.returncode, 2)

    def test_vertical_constant_refractivity_and_default_header(self):
        row = self.success()[0]
        self.assertEqual(len(row), 8)
        self.assertAlmostEqual(row[7], 2000.0, delta=2e-7)
        self.assertAlmostEqual(row[6], refractivity(100000, 288, 0) * 0.002, delta=1e-10)
        text = self.output.read_text()
        for marker in (EPOCH, "step_m=100", "truncated straight-ray", "no bending or above-top delay",
                       "truncated_los_delay_m path_length_m"):
            self.assertIn(marker, text)

    def test_dry_pressure_units(self):
        self.cube(pressure=80000, temperature=300)
        row = self.success()[0]
        self.assertAlmostEqual(row[6], 77.6 * 800 / 300 * 0.002, delta=1e-10)

    def test_wet_refractivity_units(self):
        self.cube(pressure=95000, temperature=280, humidity=0.012)
        row = self.success()[0]
        expected = refractivity(95000, 280, 0.012) * 0.002
        self.assertAlmostEqual(row[6], expected, delta=1e-10)
        self.assertGreater(expected, refractivity(95000, 280, 0) * 0.002)

    def test_equatorial_slant_independent_circle_intersection(self):
        h, top, up = 350.0, 2000.0, 0.6
        row = self.success(rays=f"0 0 {h} 0.8 0 {up}\n")[0]
        b = (A + h) * up
        c = (top - h) * (2 * A + top + h)
        expected = c / (math.sqrt(b * b + c) + b)
        self.assertAlmostEqual(row[7], expected, delta=2e-7)
        self.assertAlmostEqual(row[6], refractivity(100000, 288, 0) * expected * 1e-6, delta=1e-10)
        self.assertGreater(abs(expected - (top - h) / up), 0.5)

    def test_slant_independent_wgs84_ellipsoid_quadratic(self):
        self.cube(heights=(-1500, -500, 0), latitudes=(44, 45, 46), longitudes=(10, 11, 12))
        lon, lat, h = math.radians(11), math.radians(45), -1000.0
        east, north, up = 0.4, 0.5, math.sqrt(0.59)
        row = self.success(rays=f"11 45 {h} {east} {north} {up:.17g}\n")[0]
        e2 = F * (2 - F)
        radius = A / math.sqrt(1 - e2 * math.sin(lat)**2)
        origin = np.array([(radius + h) * math.cos(lat) * math.cos(lon),
                           (radius + h) * math.cos(lat) * math.sin(lon),
                           (radius * (1 - e2) + h) * math.sin(lat)])
        enu_basis = np.array([[-math.sin(lon), -math.sin(lat) * math.cos(lon), math.cos(lat) * math.cos(lon)],
                              [math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat) * math.sin(lon)],
                              [0, math.cos(lat), math.sin(lat)]])
        direction = enu_basis @ np.array([east, north, up])
        axes2 = np.array([A**2, A**2, (A * (1 - F))**2])
        qa = np.sum(direction**2 / axes2)
        qb = 2 * np.sum(origin * direction / axes2)
        qc = np.sum(origin**2 / axes2) - 1
        expected = -2 * qc / (qb + math.sqrt(qb * qb - 4 * qa * qc))
        self.assertAlmostEqual(row[7], expected, delta=1e-6)
        self.assertAlmostEqual(row[6], refractivity(100000, 288, 0) * expected * 1e-6, delta=1e-9)

    def test_start_heights_order_and_comments(self):
        rows = self.success(rays="# heading\n\n0 0 0 0 0 1 # first\n0 0 750 0 0 1\n0 0 1500 0 0 1\n")
        np.testing.assert_allclose(rows[:, 2], [0, 750, 1500], atol=0)
        np.testing.assert_allclose(rows[:, 7], [2000, 1250, 500], atol=2e-7, rtol=0)
        np.testing.assert_allclose(rows[:, 6], refractivity(100000, 288, 0) * np.array([.002, .00125, .0005]), atol=1e-10)

    def test_vertical_at_high_latitude(self):
        self.cube(latitudes=(66, 67, 68), longitudes=(20, 21, 22))
        row = self.success(rays="21 67 575 0 0 1\n")[0]
        self.assertAlmostEqual(row[7], 1425, delta=2e-7)

    def test_linear_height_profile_analytic_integral(self):
        heights = np.array([0.0, 450, 1300, 2000])
        self.cube(heights=heights, pressure=(100000 - 10 * heights)[:, None, None])
        row = self.success(rays="0 0 325 0 0 1\n", step=137)[0]
        expected = (refractivity(96750, 288, 0) + refractivity(80000, 288, 0)) / 2 * 1675e-6
        self.assertAlmostEqual(row[6], expected, delta=1e-10)

    def test_lateral_longitude_gradient_changes_slant(self):
        self.cube(pressure=np.array([70000, 100000, 130000])[None, None, :])
        rows = self.success(rays="0 0 0 0.8 0 0.6\n0 0 0 -0.8 0 0.6\n")
        self.assertAlmostEqual(rows[0, 7], rows[1, 7], delta=1e-7)
        self.assertGreater(rows[0, 6], rows[1, 6] + 0.001)
        self.assertAlmostEqual(np.mean(rows[:, 6]), refractivity(100000, 288, 0) * rows[0, 7] * 1e-6, delta=1e-10)

    def test_lateral_latitude_gradient_changes_slant(self):
        self.cube(pressure=np.array([70000, 100000, 130000])[None, :, None])
        rows = self.success(rays="0 0 0 0 0.8 0.6\n0 0 0 0 -0.8 0.6\n")
        self.assertGreater(rows[0, 6], rows[1, 6] + 0.001)

    def test_step_refinement_converges_to_analytic_profile(self):
        self.cube(temperature=np.array([250, 275, 300])[:, None, None])
        exact = 77.6 * 1000 / 0.025 * math.log(300 / 250) * 1e-6
        errors = []
        for step in (1000, 200, 50):
            row = self.success(step=step)[0]
            errors.append(abs(row[6] - exact))
            self.output.unlink()
        self.assertLess(errors[1], errors[0] / 20)
        self.assertLess(errors[2], errors[1] / 10)
        self.assertLess(errors[2], 3e-6)

    def test_top_layer_zero_integral(self):
        rows = self.success(rays="0 0 2000 0 0 1\n1 1 2000 0.8 0 0.6\n", step="1e-12")
        np.testing.assert_array_equal(rows[:, 6:], np.zeros((2, 2)))

    def test_near_unit_los_is_normalized(self):
        row = self.success(rays="0 0 0 0 0 1.0000005\n")[0]
        self.assertEqual(row[5], 1.0000005)
        self.assertAlmostEqual(row[7], 2000, delta=2e-7)

    def test_later_horizontal_exit_removes_partial_output(self):
        self.failure("longitude", rays="0 0 0 0 0 1\n0.999 0 0 0.8 0 0.6\n")

    def test_start_outside_model_rejected(self):
        for ray in ("2 0 0 0 0 1", "0 2 0 0 0 1", "0 0 -1 0 0 1", "0 0 2001 0 0 1"):
            with self.subTest(ray=ray):
                self.failure("outside", rays=ray)

    def test_longitude_seam_crossing_rejected(self):
        self.cube(longitudes=(-180, 0, 180))
        self.failure("seam", rays="179.999 0 0 0.8 0 0.6\n")

    def test_nonunit_horizontal_and_downward_los_rejected(self):
        for vector in ("1 0 1", "0 0 0.9", "0 0 -1", "1 0 0", "0 0 0"):
            with self.subTest(vector=vector):
                self.failure("LOS ENU", rays="0 0 0 " + vector)

    def test_nonfinite_ray_values_rejected(self):
        for position in range(6):
            for invalid in ("NaN", "inf", "-inf"):
                with self.subTest(position=position, invalid=invalid):
                    values = ["0", "0", "0", "0", "0", "1"]
                    values[position] = invalid
                    self.failure("finite", rays=" ".join(values))

    def test_invalid_geographic_rays_and_columns_rejected(self):
        for ray in ("181 0 0 0 0 1", "0 91 0 0 0 1", "0 -91 0 0 0 1", "-181 0 0 0 0 1",
                    "0 0 0 0 1", "0 0 0 0 0 1 extra", "0 0 0 text 0 1", "# no rays\n", ""):
            with self.subTest(ray=ray):
                self.failure(rays=ray)

    def test_invalid_step_and_resource_limit(self):
        for step in ("0", "-100", "nan", "inf", "1e999", "100m", "", "1e-12"):
            with self.subTest(step=step):
                self.failure(step=step)

    def test_wrong_dimension_order_or_rank(self):
        for name in ("pressure", "temperature", "specific_humidity"):
            for dims in (("latitude", "height", "longitude"), ("height", "latitude")):
                with self.subTest(name=name, dims=dims):
                    self.cube(dims={name: dims})
                    self.failure("dimension")

    def test_wrong_or_missing_units(self):
        for name in UNITS:
            with self.subTest(name=name):
                self.cube(units={name: "wrong"})
                self.failure("units")
                self.cube()
                with netcdf_file(self.weather, "a") as nc:
                    del nc.variables[name]._attributes["units"]
                self.failure("units")

    def test_missing_or_wrong_height_reference(self):
        for value in (None, "geoid", "WGS84"):
            with self.subTest(value=value):
                self.cube()
                with netcdf_file(self.weather, "a") as nc:
                    if value is None:
                        del nc._attributes["height_reference"]
                    else:
                        nc.height_reference = value
                self.failure("height_reference")

    def test_epoch_format_calendar_and_single_time(self):
        bad = (None, "", "2023-02-29T00:00:00Z", "1900-02-29T00:00:00Z", "2024-04-31T00:00:00Z",
               "2024-13-01T00:00:00Z", "2024-00-01T00:00:00Z", "0000-01-01T00:00:00Z",
               "2024-01-00T00:00:00Z", "2024-01-01T24:00:00Z", "2024-01-01T00:60:00Z",
               "2024-01-01T00:00:60Z", "2024-01-01 00:00:00Z", "2024-01-01T00:00:00+00:00",
               "2024-01-01T00:00:00Z,2024-01-02T00:00:00Z")
        for epoch in bad:
            with self.subTest(epoch=epoch):
                self.cube()
                with netcdf_file(self.weather, "a") as nc:
                    if epoch is None:
                        del nc._attributes["epoch_utc"]
                    else:
                        nc.epoch_utc = epoch
                self.failure("epoch_utc")
        with netcdf_file(self.weather, "a") as nc:
            nc.epoch_utc = "2000-02-29T00:00:00Z"
        self.success()

    def test_coordinate_order_size_and_range(self):
        for kwargs in ({"heights": (2000, 1000, 0)}, {"heights": (0, 0, 2000)}, {"heights": (0,)},
                       {"latitudes": (1, 0, -1)}, {"latitudes": (-91, 0, 1)},
                       {"longitudes": (-1, 0, 181)}, {"longitudes": (179, -179)},
                       {"heights": (0, float("nan"), 2000)}):
            with self.subTest(kwargs=kwargs):
                self.cube(**kwargs)
                self.failure()

    def test_nonfinite_and_unphysical_cube_values(self):
        for name, invalid in (("pressure", 0), ("pressure", -1), ("pressure", np.inf),
                              ("temperature", 0), ("temperature", -273), ("temperature", np.nan),
                              ("specific_humidity", -0.001), ("specific_humidity", 1),
                              ("specific_humidity", np.nan)):
            with self.subTest(name=name, invalid=invalid):
                self.cube()
                with netcdf_file(self.weather, "a") as nc:
                    nc.variables[name][0, 0, 0] = invalid
                self.failure()

    def test_missing_and_default_fill_values(self):
        for name, attr, fill in (("pressure", "_FillValue", -9999.0),
                                 ("specific_humidity", "missing_value", 0.5),
                                 ("temperature", None, 9.969209968386869e36),
                                 ("height", "_FillValue", 1000.0)):
            with self.subTest(name=name, attr=attr):
                self.cube()
                with netcdf_file(self.weather, "a") as nc:
                    var = nc.variables[name]
                    if attr:
                        setattr(var, attr, np.array([fill, -9998]) if attr == "missing_value" else fill)
                    var[(1,) if name == "height" else (0, 0, 0)] = fill
                self.failure("missing/fill")

    def test_unused_fill_metadata_is_allowed(self):
        with netcdf_file(self.weather, "a") as nc:
            nc.variables["pressure"]._FillValue = -9999.0
        self.success()

    def test_explicit_fill_overrides_type_default(self):
        with netcdf_file(self.weather, "a") as nc:
            del nc.variables["longitude"]
            var = nc.createVariable("longitude", "b", ("longitude",))
            var.units = "degrees_east"
            var._FillValue = np.int8(-128)
            var[:] = [-127, -126, -125]
        row = self.success(rays="-126 0 0 0 0 1\n")[0]
        self.assertAlmostEqual(row[7], 2000, delta=2e-7)

    def test_packed_fields_and_coordinates_rejected(self):
        for name in UNITS:
            for attr, value in (("scale_factor", 1.0), ("add_offset", 0.0)):
                with self.subTest(name=name, attr=attr):
                    self.cube()
                    with netcdf_file(self.weather, "a") as nc:
                        setattr(nc.variables[name], attr, value)
                    self.failure("packed input")

    def test_oversized_cube_rejected_from_header_without_allocating(self):
        with netcdf_file(self.weather, "w") as nc:
            nc.height_reference = "WGS84_ellipsoid"
            nc.epoch_utc = EPOCH
            for name in AXES:
                nc.createDimension(name, 200)
        self.assertLess(self.weather.stat().st_size, 1024)
        self.failure("128 MiB")

    def test_missing_variable_and_invalid_netcdf(self):
        with netcdf_file(self.weather, "a") as nc:
            del nc.variables["specific_humidity"]
        self.failure("specific_humidity")
        self.weather.write_text("not a weather cube")
        self.failure("open weather cube")

    def test_existing_output_protected(self):
        self.output.write_text("do not replace\n")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.output.read_text(), "do not replace\n")
        self.assertEqual(list(self.folder.glob("*.tmp.*")), [])

    def test_input_output_aliases_protected(self):
        original_cube = self.weather.read_bytes()
        result = self.invoke(output=self.weather)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.weather.read_bytes(), original_cube)
        result = self.invoke(output=self.rays)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.rays.read_text(), "0 0 0 0 0 1\n")
        self.output.symlink_to(self.weather)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.weather.read_bytes(), original_cube)
        self.output.unlink()
        os.link(self.weather, self.output)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.weather.read_bytes(), original_cube)
        self.assertEqual(list(self.folder.glob("*.tmp.*")), [])

    def test_dangling_output_symlink_not_overwritten(self):
        target = self.folder / "absent destination"
        self.output.symlink_to(target)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(self.output.is_symlink())
        self.assertFalse(target.exists())

    def wrapper_layout(self, name):
        directory = self.folder / name / "csh scripts"
        directory.mkdir(parents=True)
        wrapper = directory / "raytrace_atm.csh"
        shutil.copy2(WRAPPER, wrapper)
        return wrapper

    @unittest.skipUnless(shutil.which("csh"), "csh required for wrapper tests")
    def test_wrapper_same_directory_parent_and_path_with_spaces(self):
        self.assertTrue(os.access(WRAPPER, os.X_OK))
        for location in ("same", "parent", "path"):
            with self.subTest(location=location):
                wrapper = self.wrapper_layout(location)
                env = os.environ.copy()
                if location == "same":
                    (wrapper.parent / "atm_raytrace").symlink_to(self.binary)
                elif location == "parent":
                    (wrapper.parent.parent / "atm_raytrace").symlink_to(self.binary)
                else:
                    env["PATH"] = str(self.binary.parent) + os.pathsep + env.get("PATH", "")
                row = self.success(command=[wrapper], env=env, step=125)[0]
                self.assertAlmostEqual(row[7], 2000, delta=2e-7)
                self.assertIn("step_m=125", self.output.read_text())
                self.output.unlink()

    @unittest.skipUnless(shutil.which("csh"), "csh required for wrapper tests")
    def test_wrapper_preserves_binary_exit_status_and_lookup_precedence(self):
        wrapper = self.wrapper_layout("exit status")
        local = wrapper.parent / "atm_raytrace"
        parent = wrapper.parent.parent / "atm_raytrace"
        local.write_text("#!/bin/sh\nexit 37\n")
        parent.write_text("#!/bin/sh\nexit 42\n")
        local.chmod(0o755)
        parent.chmod(0o755)
        result = self.invoke(command=[wrapper])
        self.assertEqual(result.returncode, 37)
        local.unlink()
        result = self.invoke(command=[wrapper])
        self.assertEqual(result.returncode, 42)
        self.assertFalse(self.output.exists())

    @unittest.skipUnless(shutil.which("csh"), "csh required for wrapper tests")
    def test_wrapper_reports_real_validation_error_and_no_auto_build(self):
        wrapper = self.wrapper_layout("validation")
        (wrapper.parent / "atm_raytrace").symlink_to(self.binary)
        direct = self.failure(step=-1)
        wrapped = self.failure(step=-1, command=[wrapper])
        self.assertEqual(wrapped.returncode, direct.returncode)
        self.assertEqual(wrapped.stderr, direct.stderr)
        (wrapper.parent / "atm_raytrace").unlink()
        env = os.environ.copy()
        env["PATH"] = str(wrapper.parent)
        result = self.invoke(command=[wrapper], env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sorted(p.name for p in wrapper.parent.iterdir()), ["raytrace_atm.csh"])


if __name__ == "__main__":
    unittest.main()
