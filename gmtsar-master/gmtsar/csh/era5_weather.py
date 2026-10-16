#!/usr/bin/env python3
"""Build normalized atm_raytrace cubes from ARCO-ERA5.

ERA5 dates in 2026 are ERA5T and therefore preliminary.
"""

import argparse
from datetime import datetime, timedelta, timezone
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.io import netcdf_file


HEIGHTS = np.arange(0.0, 32000.0 + 1e-6, 200.0)
ARCO_ZARR = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
FETCH_VARIABLES = (
    "geopotential",
    "temperature",
    "specific_humidity",
    "surface_pressure",
    "geopotential_at_surface",
    "2m_temperature",
    "2m_dewpoint_temperature",
)


def parse_epoch(value):
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _time_seconds(values):
    values = np.asarray(values)
    if np.issubdtype(values.dtype, np.datetime64):
        seconds = values.astype("datetime64[ns]").astype(np.int64) / 1e9
    elif values.dtype.kind in "SU":
        seconds = values.astype("datetime64[ns]").astype(np.int64) / 1e9
    elif np.issubdtype(values.dtype, np.number):
        seconds = values.astype(np.float64)
    else:
        try:
            seconds = values.astype("datetime64[ns]").astype(np.int64) / 1e9
        except (TypeError, ValueError):
            seconds = np.array([
                parse_epoch(str(value)).timestamp() for value in values.ravel()
            ], dtype=np.float64).reshape(values.shape)
    return np.asarray(seconds, dtype=np.float64).reshape(-1)


def _time_bracket(raw, epoch_utc):
    epoch = epoch_utc
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    timestamp = epoch.astimezone(timezone.utc).timestamp()
    times = _time_seconds(raw["time"])
    if times.size == 0 or np.any(np.diff(times) <= 0):
        raise ValueError("Raw ERA5 time coordinate must be nonempty and strictly increasing")
    exact = np.flatnonzero(np.isclose(times, timestamp, rtol=0, atol=1e-6))
    if exact.size:
        idx = int(exact[0])
        return idx, idx, 0.0
    after = int(np.searchsorted(times, timestamp, side="right"))
    before = after - 1
    if before < 0 or after >= times.size:
        raise ValueError(f"ERA5 raw times do not bracket epoch {epoch.isoformat()}")
    weight = (timestamp - times[before]) / (times[after] - times[before])
    return before, after, float(weight)


def _select_bbox(raw, bbox):
    west, east, south, north = map(float, bbox)
    if west >= east or south >= north:
        raise ValueError("bbox must satisfy west < east and south < north")
    lat_all = np.asarray(raw["latitude"], dtype=float)
    lon_all = np.asarray(raw["longitude"], dtype=float)
    tolerance = 1e-9
    if (west < np.nanmin(lon_all) - tolerance or east > np.nanmax(lon_all) + tolerance
            or south < np.nanmin(lat_all) - tolerance or north > np.nanmax(lat_all) + tolerance):
        raise ValueError(
            f"Requested bbox {west}/{east}/{south}/{north} is outside raw ERA5 extent "
            f"{np.nanmin(lon_all)}/{np.nanmax(lon_all)}/{np.nanmin(lat_all)}/{np.nanmax(lat_all)}"
        )
    lat_idx = np.flatnonzero((lat_all >= south - tolerance) & (lat_all <= north + tolerance))
    lon_idx = np.flatnonzero((lon_all >= west - tolerance) & (lon_all <= east + tolerance))
    if lat_idx.size < 2 or lon_idx.size < 2:
        raise ValueError("Requested bbox must include at least two raw ERA5 latitude and longitude nodes")
    lat_idx = lat_idx[np.argsort(lat_all[lat_idx])]
    lon_idx = lon_idx[np.argsort(lon_all[lon_idx])]
    return lat_all[lat_idx], lon_all[lon_idx], lat_idx, lon_idx


def _weighted(raw, name, t0, t1, weight):
    values = np.asarray(raw[name])
    if values.ndim >= 3 and values.shape[0] == len(np.asarray(raw["time"]).reshape(-1)):
        return (1.0 - weight) * values[t0] + weight * values[t1]
    return values


def _geoid_path(explicit):
    if explicit:
        return Path(explicit)
    if os.environ.get("GMTSAR_GEOID"):
        return Path(os.environ["GMTSAR_GEOID"])
    sharedir = shutil.which("gmtsar_sharedir.csh")
    if sharedir:
        try:
            root = subprocess.check_output([sharedir], text=True, stderr=subprocess.PIPE).strip()
            if root:
                candidate = Path(root) / "geoid_egm96_icgem.grd"
                if candidate.is_file():
                    return candidate
        except (OSError, subprocess.CalledProcessError):
            pass
    raise FileNotFoundError(
        "Cannot resolve geoid grid; pass --geoid or set GMTSAR_GEOID "
        "(gmtsar_sharedir.csh must also be available)"
    )


def _geoid_interpolator(path):
    from osgeo import gdal

    ds = gdal.Open(str(path))
    if ds is None:
        raise FileNotFoundError(f"Could not open geoid grid: {path}")
    gt = ds.GetGeoTransform()
    values = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    lons = gt[0] + np.arange(ds.RasterXSize) * gt[1]
    lats = gt[3] + np.arange(ds.RasterYSize) * gt[5]
    if gt[5] < 0:
        lats = lats[::-1]
        values = values[::-1]
    return RegularGridInterpolator(
        (lats, lons), values, bounds_error=False, fill_value=None
    )


def _geom_from_gpm(z_gpm, latitude):
    sin_lat = math.sin(math.radians(latitude))
    gravity = 9.780327 * (
        1 + 0.0053024 * sin_lat**2
        - 0.0000058 * math.sin(2 * math.radians(latitude))**2
    )
    radius = 6356766.0
    return radius * z_gpm / ((gravity / 9.80665) * radius - z_gpm)


def _iso_utc(value):
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_cube(raw, epoch_utc, bbox, geoid_path=None, out_nc=None):
    """Build an atm_raytrace NetCDF cube for an epoch and geographic bbox."""
    if out_nc is None:
        raise ValueError("out_nc is required")
    if not isinstance(epoch_utc, datetime):
        epoch_utc = parse_epoch(str(epoch_utc))
    if epoch_utc.tzinfo is None:
        epoch_utc = epoch_utc.replace(tzinfo=timezone.utc)
    t0, t1, weight = _time_bracket(raw, epoch_utc)
    lats, lons, lat_idx, lon_idx = _select_bbox(raw, bbox)
    lev = np.asarray(raw["level"], dtype=float)
    geo = _geoid_interpolator(_geoid_path(geoid_path))

    def field(name):
        values = _weighted(raw, name, t0, t1, weight)
        if name in ("geopotential", "temperature", "specific_humidity"):
            return np.asarray(values)[:, lat_idx, :][:, :, lon_idx]
        return np.asarray(values)[lat_idx, :][:, lon_idx]

    geopotential = field("geopotential") / 9.80665
    temperature = field("temperature")
    humidity = field("specific_humidity")
    surface_pressure = field("surface_pressure")
    surface_geopotential = _weighted(raw, "geopotential_at_surface", t0, t1, weight)
    surface_geopotential = np.asarray(surface_geopotential)[lat_idx, :][:, lon_idx] / 9.80665
    t2m = field("2m_temperature")
    dewpoint = field("2m_dewpoint_temperature")
    nh = len(HEIGHTS)
    shape = (nh, len(lats), len(lons))
    pressure_cube = np.empty(shape, dtype=np.float64)
    temperature_cube = np.empty(shape, dtype=np.float64)
    humidity_cube = np.empty(shape, dtype=np.float64)

    for jj, j in enumerate(lat_idx):
        for ii, i in enumerate(lon_idx):
            geoid = float(geo([[lats[jj], lons[ii]]])[0])
            sp = float(surface_pressure[jj, ii])
            td = float(dewpoint[jj, ii]) - 273.15
            vapor_pressure = 611.2 * math.exp(17.67 * td / (td + 243.5))
            qs = 0.622 * vapor_pressure / (sp - 0.378 * vapor_pressure)
            heights = [_geom_from_gpm(float(surface_geopotential[jj, ii]), lats[jj]) + geoid]
            pressures = [sp]
            temps = [float(t2m[jj, ii])]
            humidities = [max(1e-12, qs)]
            for k, level_hpa in enumerate(lev):
                level_pressure = float(level_hpa) * 100.0
                if level_pressure >= sp:
                    continue
                heights.append(
                    _geom_from_gpm(float(geopotential[k, jj, ii]), lats[jj]) + geoid
                )
                pressures.append(level_pressure)
                temps.append(float(temperature[k, jj, ii]))
                humidities.append(max(1e-12, float(humidity[k, jj, ii])))
            order = np.argsort(heights)
            heights = np.asarray(heights, dtype=float)[order]
            pressures = np.asarray(pressures, dtype=float)[order]
            temps = np.asarray(temps, dtype=float)[order]
            humidities = np.asarray(humidities, dtype=float)[order]
            heights, unique = np.unique(heights, return_index=True)
            pressures, temps, humidities = pressures[unique], temps[unique], humidities[unique]
            if len(heights) < 2:
                raise ValueError("ERA5 column has fewer than two usable vertical levels")
            pressure_cube[:, jj, ii] = np.exp(
                interp1d(heights, np.log(pressures), fill_value="extrapolate")(HEIGHTS)
            )
            temperature_cube[:, jj, ii] = np.maximum(
                150.0, interp1d(heights, temps, fill_value="extrapolate")(HEIGHTS)
            )
            q = np.exp(interp1d(heights, np.log(humidities), fill_value="extrapolate")(HEIGHTS))
            q[HEIGHTS > 15000.0] = 1e-12
            humidity_cube[:, jj, ii] = np.clip(q, 0.0, 0.05)

    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    with netcdf_file(str(out_nc), "w") as nc:
        nc.height_reference = "WGS84_ellipsoid"
        nc.epoch_utc = _iso_utc(epoch_utc)
        nc.createDimension("height", nh)
        nc.createDimension("latitude", len(lats))
        nc.createDimension("longitude", len(lons))
        for name, values, units in (
            ("height", HEIGHTS, "m"),
            ("latitude", lats, "degrees_north"),
            ("longitude", lons, "degrees_east"),
        ):
            variable = nc.createVariable(name, "d", (name,))
            variable.units = units
            variable[:] = values
        for name, values, units in (
            ("pressure", pressure_cube, "Pa"),
            ("temperature", temperature_cube, "K"),
            ("specific_humidity", humidity_cube, "kg kg-1"),
        ):
            variable = nc.createVariable(
                name, "d", ("height", "latitude", "longitude")
            )
            variable.units = units
            variable[:] = values
    return out_nc


def _parse_bbox(value):
    try:
        parts = tuple(float(part) for part in value.split("/"))
    except ValueError as error:
        raise argparse.ArgumentTypeError("bbox must be W/E/S/N") from error
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("bbox must be W/E/S/N")
    return parts


def fetch_raw(bbox, epochs, out):
    try:
        import xarray as xr
    except ImportError:
        print('Install ERA5 fetch dependencies: pip install xarray "zarr<3" gcsfs', file=sys.stderr)
        return 2
    ds = xr.open_zarr(ARCO_ZARR, chunks=None, storage_options={"token": "anon"})
    west, east, south, north = bbox
    lat_values = ds["latitude"].values
    lat_slice = slice(north, south) if lat_values[0] > lat_values[-1] else slice(south, north)
    hours = sorted({
        hour
        for epoch in epochs
        for hour in (
            epoch.replace(minute=0, second=0, microsecond=0),
            epoch.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1),
        )
    })
    raw_subset = ds[list(FETCH_VARIABLES)].sel(
        time=np.asarray([np.datetime64(hour.replace(tzinfo=None), "ns") for hour in hours]),
        latitude=lat_slice,
        longitude=slice(west, east),
    )
    raw = {
        "time": np.asarray(raw_subset["time"].values),
        "level": np.asarray(ds["level"].values),
        "latitude": np.asarray(raw_subset["latitude"].values),
        "longitude": np.asarray(raw_subset["longitude"].values),
    }
    for name in FETCH_VARIABLES:
        variable = raw_subset[name]
        if "level" in variable.dims:
            variable = variable.transpose("time", "level", "latitude", "longitude")
        elif "time" in variable.dims:
            variable = variable.transpose("time", "latitude", "longitude")
        raw[name] = np.asarray(variable.values)
    np.savez_compressed(out, **raw)
    ds.close()
    print(f"Saved {out} with {len(hours)} hourly samples")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("fetch", help="fetch hourly ARCO-ERA5 samples")
    fetch.add_argument("--bbox", required=True, type=_parse_bbox, help="W/E/S/N")
    fetch.add_argument("--epoch", required=True, action="append", help="UTC ISO timestamp; repeatable")
    fetch.add_argument("--out", required=True)
    build = commands.add_parser("build", help="build a normalized atm_raytrace cube")
    build.add_argument("--raw", required=True)
    build.add_argument("--epoch", required=True)
    build.add_argument("--bbox", required=True, type=_parse_bbox, help="W/E/S/N")
    build.add_argument("--geoid")
    build.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    if args.command == "fetch":
        epochs = [parse_epoch(value) for value in args.epoch]
        return fetch_raw(args.bbox, epochs, args.out)
    try:
        with np.load(args.raw, allow_pickle=False) as archive:
            raw = {name: archive[name] for name in archive.files}
        build_cube(raw, parse_epoch(args.epoch), args.bbox, args.geoid, args.out)
    except (OSError, ValueError, KeyError) as error:
        print(f"era5_weather.py: {error}", file=sys.stderr)
        return 2
    print(f"Built weather cube {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
