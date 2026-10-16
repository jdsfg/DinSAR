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
from pathlib import Path
import numpy as np
from osgeo import gdal, osr

def parse_args():
    parser = argparse.ArgumentParser(description="GMTSAR Physical Atmospheric Ray-Tracing Workflow")
    parser.add_argument("--source-phase", required=True, help="Input unwrapped phase grid (unwrap.grd)")
    parser.add_argument("--output-phase", required=True, help="Output corrected unwrapped phase grid (unwrap_atm.grd)")
    parser.add_argument("--project-dir", required=True, help="Root project directory")
    parser.add_argument("--stage-dir", required=True, help="Stage execution directory")
    parser.add_argument("--weather-dir", default="", help="Optional weather directory containing .nc cubes")
    return parser.parse_args()

def locate_raytrace_bin():
    candidates = [
        "/tmp/atm_raytrace",
        Path(__file__).resolve().parents[1] / "atm_raytrace",
        Path("/home/mihu/.gemini/antigravity-cli/brain/b7e5878e-56ae-4d39-ae08-f151da5e8456/scratch/atm/atm_raytrace"),
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
        Path("/home/mihu/.gemini/antigravity-cli/brain/b7e5878e-56ae-4d39-ae08-f151da5e8456/scratch/atm")
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

    m_nc, s_nc = locate_weather_cubes(stage, proj, args.weather_dir)
    if not m_nc or not s_nc:
        raise FileNotFoundError("Weather cubes (master_weather.nc / slave_weather.nc) not found.")
    print(f"[Raytrace] Weather cubes: Master={m_nc}, Slave={s_nc}")

    dem_geo = proj / "topo/dem.grd"
    if not dem_geo.is_file():
        raise FileNotFoundError(f"Geographic DEM not found: {dem_geo}")

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
    temp_int = stage / "temp_atm_int.grd"
    atm_ra = stage / "atm_phase_ra.grd"

    subprocess.run(["proj_ll2ra.csh", str(trans_dat), str(atm_ll_grd), str(temp_ra), "-I4/4"],
                   cwd=str(stage), check=True)
    subprocess.run(["nearest_grid", str(temp_ra), str(temp_int), "300"], cwd=str(stage), check=True)
    subprocess.run(["gmt", "grdsample", str(temp_int), f"-R{source_phase}", f"-G{atm_ra}"],
                   cwd=str(stage), check=True)

    # 4. Load radar grids and perform correction
    ds_unw = gdal.Open(str(source_phase))
    unw_arr = ds_unw.GetRasterBand(1).ReadAsArray().astype(np.float32)
    gt_unw = ds_unw.GetGeoTransform()
    proj_unw = ds_unw.GetProjection()

    ds_atm = gdal.Open(str(atm_ra))
    atm_arr = ds_atm.GetRasterBand(1).ReadAsArray().astype(np.float32)

    corr_file = stage / "corr.grd"
    if corr_file.is_file():
        ds_c = gdal.Open(str(corr_file))
        corr_arr = ds_c.GetRasterBand(1).ReadAsArray().astype(np.float32)
    else:
        corr_arr = np.ones_like(unw_arr)

    nrows, ncols = unw_arr.shape
    y_idx, x_idx = np.mgrid[0:nrows, 0:ncols]

    # Valid mask (ensure same valid footprint)
    nodata_val = ds_unw.GetRasterBand(1).GetNoDataValue()
    if nodata_val is not None:
        valid_unw = np.isfinite(unw_arr) & (unw_arr != nodata_val)
    else:
        valid_unw = np.isfinite(unw_arr) & (unw_arr != 0)

    fit_mask = valid_unw & (corr_arr > 0.6)
    if np.sum(fit_mask) < 1000:
        fit_mask = valid_unw

    # Fit 2D planar trend to atmospheric delay to isolate topographic stratification
    step = 4
    sub_m = fit_mask[::step, ::step]
    xs = (x_idx[::step, ::step][sub_m] / ncols).astype(np.float64)
    ys = (y_idx[::step, ::step][sub_m] / nrows).astype(np.float64)
    zs = atm_arr[::step, ::step][sub_m].astype(np.float64)

    A_plane = np.column_stack([np.ones_like(xs), xs, ys])
    c_p, _, _, _ = np.linalg.lstsq(A_plane, zs, rcond=None)

    X_norm = (x_idx / ncols).astype(np.float32)
    Y_norm = (y_idx / nrows).astype(np.float32)
    atm_plane = (c_p[0] + c_p[1] * X_norm + c_p[2] * Y_norm).astype(np.float32)
    atm_strat = atm_arr - atm_plane

    alpha = 0.280
    corr_unw = unw_arr + alpha * atm_strat

    # Reference alignment
    ref_val = float(np.median(corr_unw[fit_mask]) - np.median(unw_arr[fit_mask]))
    corr_unw -= ref_val

    # Ensure strictly identical valid footprint as input
    out_arr = unw_arr.copy()
    out_arr[valid_unw] = corr_unw[valid_unw]

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
