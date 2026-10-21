#!/usr/bin/env python3
"""
raytrace_atm_workflow.py:
End-to-end ray-tracing physical atmospheric correction workflow for GMTSAR / run_dinsar.sh.
Called by Raytrace_atm.csh in the pipeline.

Workflow:
1. Parse master and slave PRM files (wavelength, acquisition UTC timestamps).
2. Locate or download/generate WGS84 NetCDF weather cubes (NOAA GFS).
3. Sample DEM and compute 3D look vectors via SAT_look.
4. Run atm_raytrace numerical integration along LOS for master and slave.
5. Compute differential phase, grid to geographic domain, and project to radar coordinates.
6. Isolate topographic stratification (remove planar trend).
7. Apply referenced correction to unwrapped phase and export unwrap_atm.grd.
"""

import sys
import os
import math
import shutil
import argparse
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
from osgeo import gdal
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

def parse_args():
    parser = argparse.ArgumentParser(description="GMTSAR Physical Atmospheric Ray-Tracing Workflow")
    parser.add_argument("--source-phase", required=True, help="Input unwrapped phase grid (unwrap.grd)")
    parser.add_argument("--output-phase", required=True, help="Output corrected unwrapped phase grid (unwrap_atm.grd)")
    parser.add_argument("--project-dir", required=True, help="Root project directory")
    parser.add_argument("--stage-dir", required=True, help="Stage execution directory")
    parser.add_argument("--weather-dir", default="", help="Optional weather directory containing .nc cubes")
    parser.add_argument("--scale", type=float, default=1.0, help="Atmospheric phase correction scale")
    return parser.parse_args()

def locate_raytrace_bin():
    candidates = [
        "/tmp/atm_raytrace",
        Path(__file__).resolve().parents[1] / "atm_raytrace",
    ]
    for c in candidates:
        if os.path.isfile(str(c)) and os.access(str(c), os.X_OK):
            return str(c)

    # Try which atm_raytrace
    which_bin = shutil.which("atm_raytrace")
    if which_bin:
        return which_bin

    # Try compiling to /tmp/atm_raytrace
    cpp_src = Path(__file__).resolve().parents[1] / "atm_raytrace.cpp"
    if cpp_src.is_file():
        try:
            cmd = f"g++ -std=c++17 -O2 -Wall -Wextra $(nc-config --cflags) '{cpp_src}' -o /tmp/atm_raytrace $(nc-config --libs)"
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if os.path.isfile("/tmp/atm_raytrace") and os.access("/tmp/atm_raytrace", os.X_OK):
                return "/tmp/atm_raytrace"
        except Exception:
            pass

    raise FileNotFoundError("Could not find or compile executable atm_raytrace")

def parse_prms(stage_dir):
    prm_files = sorted(list(Path(stage_dir).glob("*.PRM")))
    if not prm_files:
        raise FileNotFoundError(f"No PRM files found in {stage_dir}")

    master_prm = prm_files[0]
    slave_prm = prm_files[1] if len(prm_files) > 1 else prm_files[0]

    wavelength = 0.0312284
    with open(master_prm) as f:
        for line in f:
            if "radar_wavelength" in line:
                parts = line.split("=")
                if len(parts) > 1:
                    wavelength = float(parts[1].strip().split()[0])
                break

    return master_prm, slave_prm, wavelength

def locate_weather_cubes(stage_dir, project_dir, weather_dir):
    candidates_dirs = []
    if weather_dir and os.path.isdir(weather_dir):
        candidates_dirs.append(Path(weather_dir))
    candidates_dirs.extend([
        Path(stage_dir) / "weather",
        Path(stage_dir),
        Path(project_dir) / "weather",
    ])

    m_nc = None
    s_nc = None
    for d in candidates_dirs:
        if not d.is_dir():
            continue
        m_c = list(d.glob("*master*weather*.nc")) or list(d.glob("master_weather.nc"))
        s_c = list(d.glob("*slave*weather*.nc")) or list(d.glob("slave_weather.nc"))
        if m_c and s_c:
            m_nc = m_c[0]
            s_nc = s_c[0]
            break

    return m_nc, s_nc


def prm_epoch_utc(prm_path):
    clock = None
    with open(prm_path) as prm:
        for line in prm:
            if line.strip().startswith("SC_clock_start"):
                clock = line.split("=", 1)[1].strip().split()[0]
                break
    if clock is None:
        raise ValueError(f"SC_clock_start not found in PRM: {prm_path}")
    try:
        year_text, day_text = clock.split(".", 1)
        year, day = int(year_text[:4]), int(year_text[4:])
        fraction = float("0." + day_text)
        start = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=day - 1)
        return start + timedelta(seconds=fraction * 86400.0)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"Invalid SC_clock_start in {prm_path}: {clock}") from error


def validate_cube_epoch(cube_path, expected_epoch, tolerance_seconds=1800):
    with netcdf_file(str(cube_path), "r", mmap=False) as nc:
        value = getattr(nc, "epoch_utc", None)
    if isinstance(value, bytes):
        value = value.decode("ascii", errors="replace")
    if value is None:
        raise ValueError(f"Weather cube has no epoch_utc attribute: {cube_path}")
    try:
        actual = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"Invalid epoch_utc in weather cube {cube_path}: {value}") from error
    if actual.tzinfo is None:
        actual = actual.replace(tzinfo=timezone.utc)
    if expected_epoch.tzinfo is None:
        expected_epoch = expected_epoch.replace(tzinfo=timezone.utc)
    difference = abs((actual - expected_epoch).total_seconds())
    if difference > tolerance_seconds:
        raise ValueError(
            f"Weather cube epoch mismatch for {cube_path}: cube={actual.isoformat()}, "
            f"PRM={expected_epoch.isoformat()}, difference={difference:.1f}s "
            f"(tolerance {tolerance_seconds}s)"
        )
    return True


def _weather_dirs(stage_dir, project_dir, weather_dir):
    candidates = []
    if weather_dir and os.path.isdir(weather_dir):
        candidates.append(Path(weather_dir))
    candidates.extend((Path(stage_dir) / "weather", Path(stage_dir), Path(project_dir) / "weather"))
    return candidates


def locate_era5_raw(stage_dir, project_dir, weather_dir):
    for directory in _weather_dirs(stage_dir, project_dir, weather_dir):
        if directory.is_dir():
            matches = sorted(directory.glob("*era5*raw*.npz"))
            if matches:
                return matches[0]
    return None


def dem_bbox(dem_path, margin=0.75):
    ds = gdal.Open(str(dem_path))
    if ds is None:
        raise FileNotFoundError(f"Could not open geographic DEM: {dem_path}")
    gt = ds.GetGeoTransform()
    if gt[2] != 0 or gt[4] != 0:
        raise ValueError(f"Rotated geographic DEM geotransform is unsupported: {dem_path}")
    x0, x1 = gt[0], gt[0] + ds.RasterXSize * gt[1]
    y0, y1 = gt[3], gt[3] + ds.RasterYSize * gt[5]
    return min(x0, x1) - margin, max(x0, x1) + margin, min(y0, y1) - margin, max(y0, y1) + margin


def ensure_weather_cubes(stage_dir, project_dir, weather_dir, master_prm, slave_prm, dem_path):
    m_nc, s_nc = locate_weather_cubes(stage_dir, project_dir, weather_dir)
    master_epoch = prm_epoch_utc(master_prm)
    slave_epoch = prm_epoch_utc(slave_prm)
    if not m_nc or not s_nc:
        raw_path = locate_era5_raw(stage_dir, project_dir, weather_dir)
        if raw_path is not None:
            from era5_weather import build_cube

            with np.load(raw_path, allow_pickle=False) as archive:
                raw = {name: archive[name] for name in archive.files}
            bbox = dem_bbox(dem_path)
            weather_out = Path(stage_dir) / "weather"
            weather_out.mkdir(parents=True, exist_ok=True)
            m_nc = build_cube(raw, master_epoch, bbox, None, weather_out / "master_weather.nc")
            s_nc = build_cube(raw, slave_epoch, bbox, None, weather_out / "slave_weather.nc")
    if not m_nc or not s_nc:
        raise FileNotFoundError(
            "Weather cubes (master_weather.nc / slave_weather.nc) not found; "
            "provide matching cubes or an *era5*raw*.npz file."
        )
    validate_cube_epoch(m_nc, master_epoch)
    validate_cube_epoch(s_nc, slave_epoch)
    return m_nc, s_nc, master_epoch, slave_epoch


def _grid_axes(ds):
    gt = ds.GetGeoTransform()
    if gt[2] != 0 or gt[4] != 0:
        raise ValueError("Rotated radar grid geotransforms are unsupported for atmospheric resampling")
    x = gt[0] + (np.arange(ds.RasterXSize) + 0.5) * gt[1]
    y = gt[3] + (np.arange(ds.RasterYSize) + 0.5) * gt[5]
    return x, y


def fill_and_resample(temp_grid, source_grid, output_grid, radius=300):
    temp = gdal.Open(str(temp_grid))
    if temp is None:
        raise FileNotFoundError(f"Could not open projected atmospheric grid: {temp_grid}")
    source = gdal.Open(str(source_grid))
    if source is None:
        raise FileNotFoundError(f"Could not open source phase grid: {source_grid}")
    values = temp.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nodata = temp.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        values[values == nodata] = np.nan
    bad = ~np.isfinite(values)
    if np.all(bad):
        filled = np.full_like(values, np.nan)
    else:
        distance, (rows, cols) = ndimage.distance_transform_edt(bad, return_indices=True)
        filled = values[rows, cols]
        filled[distance > radius] = np.nan
    x, y = _grid_axes(temp)
    if x[0] > x[-1]:
        x, filled = x[::-1], filled[:, ::-1]
    if y[0] > y[-1]:
        y, filled = y[::-1], filled[::-1, :]
    target_x, target_y = _grid_axes(source)
    xx, yy = np.meshgrid(target_x, target_y)
    resampled = RegularGridInterpolator(
        (y, x), filled, bounds_error=False, fill_value=np.nan
    )(np.column_stack((yy.ravel(), xx.ravel()))).reshape(xx.shape)
    output_grid = Path(output_grid)
    output_grid.parent.mkdir(parents=True, exist_ok=True)
    if output_grid.exists():
        output_grid.unlink()
    driver = gdal.GetDriverByName("netCDF")
    out = driver.Create(
        str(output_grid), source.RasterXSize, source.RasterYSize, 1, gdal.GDT_Float32
    )
    if out is None:
        raise RuntimeError(f"Could not create resampled atmospheric grid: {output_grid}")
    out.SetGeoTransform(source.GetGeoTransform())
    out.SetProjection(source.GetProjection())
    band = out.GetRasterBand(1)
    band.SetNoDataValue(float("nan"))
    band.WriteArray(resampled.astype(np.float32))
    out.FlushCache()
    out = None
    return resampled.astype(np.float32)


def _remove_plane(values, mask, step=4):
    nrows, ncols = values.shape
    sampled = mask[::step, ::step] & np.isfinite(values[::step, ::step])
    if sampled.sum() < 3:
        rows, cols = np.nonzero(mask & np.isfinite(values))
        zs = values[rows, cols].astype(np.float64)
        xs, ys = cols / ncols, rows / nrows
    else:
        rows, cols = np.nonzero(sampled)
        rows, cols = rows * step, cols * step
        zs = values[::step, ::step][sampled].astype(np.float64)
        xs, ys = cols / ncols, rows / nrows
    design = np.column_stack((np.ones_like(xs), xs, ys))
    coefficients, _, _, _ = np.linalg.lstsq(design, zs, rcond=None)
    residual = np.array(values, dtype=np.float32, copy=True)
    residual -= coefficients[0]
    residual -= (coefficients[1] * np.arange(ncols, dtype=np.float32) / ncols)[None, :]
    residual -= (coefficients[2] * np.arange(nrows, dtype=np.float32) / nrows)[:, None]
    return residual, coefficients


def correct_phase(unw_arr, atm_arr, valid_unw, fit_mask, scale=1.0):
    usable = valid_unw & np.isfinite(atm_arr)
    fit_mask = fit_mask & usable
    if np.count_nonzero(fit_mask) < 3:
        raise ValueError("Insufficient valid pixels to fit atmospheric plane")
    atm_strat, _ = _remove_plane(atm_arr, fit_mask)
    corrected = unw_arr.astype(np.float32) - float(scale) * atm_strat.astype(np.float32)
    ref_val = float(np.median(corrected[fit_mask]) - np.median(unw_arr[fit_mask]))
    corrected -= ref_val
    output = unw_arr.copy()
    output[usable] = corrected[usable]
    return output, atm_strat


def _plane_removed_std(values, mask):
    residual, _ = _remove_plane(values.astype(np.float64), mask)
    return float(np.std(residual[mask]))


def main():
    args = parse_args()
    stage = Path(args.stage_dir)
    proj = Path(args.project_dir)

    source_phase = Path(args.source_phase)
    output_phase = Path(args.output_phase)

    print(f"[Raytrace] Starting physical atmospheric ray-tracing workflow...")
    print(f"[Raytrace] Source: {source_phase} -> Output: {output_phase}")

    raytrace_bin = locate_raytrace_bin()
    master_prm, slave_prm, wavelength = parse_prms(stage)
    print(f"[Raytrace] PRMs: Master={master_prm.name}, Wavelength={wavelength:.5f} m")

    dem_geo = proj / "topo/dem.grd"
    if not dem_geo.is_file():
        raise FileNotFoundError(f"Geographic DEM not found: {dem_geo}")
    m_nc, s_nc, master_epoch, slave_epoch = ensure_weather_cubes(
        stage, proj, args.weather_dir, master_prm, slave_prm, dem_geo
    )
    print(f"[Raytrace] Weather cubes: Master={m_nc}, Slave={s_nc}")

    trans_dat = stage / "trans.dat"
    if not trans_dat.is_file():
        raise FileNotFoundError(f"trans.dat not found in stage: {trans_dat}")

    # 1. Generate rays: resample DEM to 15s for fast & smooth ray generation
    rays_file = stage / "rays.txt"
    sample_dem = stage / "tmp_dem_sample.grd"
    print(f"[Raytrace] Sampling DEM to 15s and computing look vectors via SAT_look...")
    subprocess.run(["gmt", "grdsample", str(dem_geo), "-I15s", f"-G{sample_dem}"], cwd=str(stage), check=True)

    xyz_proc = subprocess.run(["gmt", "grd2xyz", str(sample_dem), "-s"], cwd=str(stage),
                              capture_output=True, text=True, check=True)
    orig_lines = xyz_proc.stdout.strip().splitlines()

    sat_proc = subprocess.run(["SAT_look", str(master_prm)], cwd=str(stage),
                              input=xyz_proc.stdout, capture_output=True, text=True, check=True)
    sat_lines = sat_proc.stdout.strip().splitlines()

    with open(rays_file, "w") as f_out:
        for orig_l, sat_l in zip(orig_lines, sat_lines):
            o_parts = orig_l.split()
            s_parts = sat_l.split()
            if len(o_parts) < 3 or len(s_parts) < 6:
                continue
            lon, lat, h = float(o_parts[0]), float(o_parts[1]), float(o_parts[2])
            e, n, u = float(s_parts[3]), float(s_parts[4]), float(s_parts[5])
            norm = math.hypot(e, n, u)
            if norm == 0 or u <= 0 or not math.isfinite(norm):
                continue
            e /= norm; n /= norm; u /= norm
            f_out.write(f"{lon:.8f} {lat:.8f} {h:.3f} {e:.8f} {n:.8f} {u:.8f}\n")

    # 2. Run atm_raytrace
    delays_m = stage / "delays_master.txt"
    delays_s = stage / "delays_slave.txt"
    print(f"[Raytrace] Executing 3D ray-tracing integration...")
    subprocess.run([raytrace_bin, str(m_nc), str(rays_file), str(delays_m), "100.0"], check=True)
    subprocess.run([raytrace_bin, str(s_nc), str(rays_file), str(delays_s), "100.0"], check=True)

    # 3. Compute differential delay and phase
    m_delays, s_delays, lons, lats = [], [], [], []
    with open(delays_m) as fm, open(delays_s) as fs:
        for lm, ls in zip(fm, fs):
            if lm.startswith("#") or not lm.strip() or ls.startswith("#") or not ls.strip():
                continue
            pm = lm.split()
            ps = ls.split()
            lons.append(float(pm[0]))
            lats.append(float(pm[1]))
            m_delays.append(float(pm[6]))
            s_delays.append(float(ps[6]))

    dL = np.array(m_delays, dtype=np.float64) - np.array(s_delays, dtype=np.float64)
    # Phase delay: phi_atm = -4*pi/lambda * dL
    dPhi = - (4.0 * math.pi / wavelength) * dL

    atm_xyz = stage / "atm_phase_ll.xyz"
    with open(atm_xyz, "w") as f:
        for lon, lat, phi in zip(lons, lats, dPhi):
            f.write(f"{lon:.8f} {lat:.8f} {phi:.6f}\n")

    # Grid with surface
    ds_s = gdal.Open(str(sample_dem))
    gt_s = ds_s.GetGeoTransform()
    nx_s, ny_s = ds_s.RasterXSize, ds_s.RasterYSize
    x_min = gt_s[0]
    x_max = gt_s[0] + nx_s * gt_s[1]
    y_min = gt_s[3] + ny_s * gt_s[5]
    y_max = gt_s[3]
    ds_s = None

    atm_ll_grd = stage / "atm_phase_ll.grd"
    region = f"-R{x_min:.8f}/{x_max:.8f}/{y_min:.8f}/{y_max:.8f}"
    inc = f"-I{gt_s[1]:.8f}/{abs(gt_s[5]):.8f}"
    subprocess.run(["gmt", "surface", str(atm_xyz), region, inc, "-T0.35", f"-G{atm_ll_grd}"],
                   cwd=str(stage), check=True)

    # Project to radar coordinates
    temp_ra = stage / "temp_atm_ra.grd"
    atm_ra = stage / "atm_phase_ra.grd"

    subprocess.run(["proj_ll2ra.csh", str(trans_dat), str(atm_ll_grd), str(temp_ra), "-I4/4"],
                   cwd=str(stage), check=True)
    atm_arr = fill_and_resample(temp_ra, source_phase, atm_ra, radius=300)

    # 4. Load radar grids and perform correction
    ds_unw = gdal.Open(str(source_phase))
    unw_arr = ds_unw.GetRasterBand(1).ReadAsArray().astype(np.float32)
    gt_unw = ds_unw.GetGeoTransform()
    proj_unw = ds_unw.GetProjection()
    nrows, ncols = unw_arr.shape

    corr_file = stage / "corr.grd"
    if corr_file.is_file():
        ds_c = gdal.Open(str(corr_file))
        corr_arr = ds_c.GetRasterBand(1).ReadAsArray().astype(np.float32)
        ds_c = None
    else:
        corr_arr = np.ones_like(unw_arr)

    nodata_val = ds_unw.GetRasterBand(1).GetNoDataValue()
    if nodata_val is not None:
        valid_unw = np.isfinite(unw_arr) & (unw_arr != nodata_val)
    else:
        valid_unw = np.isfinite(unw_arr) & (unw_arr != 0)
    fit_mask = valid_unw & np.isfinite(atm_arr) & (corr_arr > 0.6)
    if np.sum(fit_mask) < 1000:
        fit_mask = valid_unw & np.isfinite(atm_arr)

    # phasediff uses SLC1*conj(SLC2), so the delay correction has this sign.
    out_arr, atm_strat = correct_phase(unw_arr, atm_arr, valid_unw, fit_mask, args.scale)
    before_std = _plane_removed_std(unw_arr, fit_mask)
    after_std = _plane_removed_std(out_arr, fit_mask)
    print(f"[Raytrace] Plane-removed unwrapped-phase std on fit_mask: before={before_std:.6f}, after={after_std:.6f} rad")
    ds_unw = None

    # Export output NetCDF matching input exact structure
    drv = gdal.GetDriverByName("netCDF")
    ds_out = drv.Create(str(output_phase), ncols, nrows, 1, gdal.GDT_Float32)
    ds_out.SetGeoTransform(gt_unw)
    ds_out.SetProjection(proj_unw)
    b_out = ds_out.GetRasterBand(1)
    if nodata_val is not None:
        b_out.SetNoDataValue(nodata_val)
    b_out.WriteArray(out_arr)
    ds_out.FlushCache()
    ds_out = None

    print(f"[Raytrace] Atmospheric correction successfully written to: {output_phase}")
    return 0

if __name__ == '__main__':
    sys.exit(main())
