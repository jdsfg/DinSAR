#!/usr/bin/env python3
import sys
import os
import numpy as np
from osgeo import gdal, osr
from scipy.interpolate import griddata

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 geocode_gdal.py <trans.dat> <input_ra.grd> <output_ll.tif> [res_deg]")
        sys.exit(1)

    trans_file = sys.argv[1]
    input_grd = sys.argv[2]
    output_tif = sys.argv[3]
    res = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0001

    print(f"[GeoCode] Reading {input_grd} and mapping via {trans_file} (Res: {res})")
    
    # 1. Read input grid
    ds_in = gdal.Open(input_grd)
    if not ds_in:
        print(f"Error: Cannot open {input_grd}")
        sys.exit(1)
        
    band_in = ds_in.GetRasterBand(1)
    val_arr = band_in.ReadAsArray()
    nrows, ncols = val_arr.shape
    ds_in = None

    # 2. Read trans.dat
    trans = np.fromfile(trans_file, dtype=np.float64).reshape(-1, 5)
    lon = trans[:, 3]
    lat = trans[:, 4]
    
    # Get original full-resolution dimensions from trans.dat max values
    full_ncols = int(np.ceil(trans[:, 0].max())) + 1
    full_nrows = int(np.ceil(trans[:, 1].max())) + 1
    
    # If the input grid is downsampled, scale coordinates accordingly
    if ncols < full_ncols:
        scale_x = ncols / full_ncols
        scale_y = nrows / full_nrows
        print(f"[GeoCode] Detected downsampled input grid ({ncols}x{nrows}). Scaling trans.dat coordinates ({full_ncols}x{full_nrows}) by x:{scale_x:.4f}, y:{scale_y:.4f}")
        rng_f = trans[:, 0] * scale_x
        azi_f = trans[:, 1] * scale_y
    else:
        rng_f = trans[:, 0]
        azi_f = trans[:, 1]
        
    rng = rng_f.astype(int)
    azi = azi_f.astype(int)

    # Valid points mask
    mask = (rng >= 0) & (rng < ncols) & (azi >= 0) & (azi < nrows)
    lon, lat = lon[mask], lat[mask]
    rng, azi = rng[mask], azi[mask]

    vals = val_arr[azi, rng]
    
    # Remove NaNs from values
    valid = ~np.isnan(vals)
    lon, lat, vals = lon[valid], lat[valid], vals[valid]

    if len(vals) == 0:
        print("Error: No valid data points found.")
        sys.exit(1)

    # 3. Create regular LL Grid
    lon_min, lon_max = lon.min(), lon.max()
    lat_min, lat_max = lat.min(), lat.max()

    out_cols = int(np.ceil((lon_max - lon_min) / res))
    out_rows = int(np.ceil((lat_max - lat_min) / res))
    print(f"[GeoCode] Target Grid Size: {out_cols} x {out_rows} (Lon: {lon_min:.4f}~{lon_max:.4f}, Lat: {lat_min:.4f}~{lat_max:.4f})")

    grid_lon, grid_lat = np.meshgrid(
        np.linspace(lon_min, lon_max, out_cols),
        np.linspace(lat_max, lat_min, out_rows)
    )

    # 4. Interpolate using Scipy (Nearest avoids interpolation artifacts on wrapped phase, safe for displacement)
    print("\tInterpolating data into grid using Nearest-Neighbor KDTree...")
    points = np.column_stack((lon, lat))
    out_grid = griddata(points, vals, (grid_lon, grid_lat), method='nearest', fill_value=np.nan)

    # 5. Export to GeoTIFF
    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(output_tif, out_cols, out_rows, 1, gdal.GDT_Float32)
    ds_out.SetGeoTransform((lon_min, res, 0, lat_max, 0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds_out.SetProjection(srs.ExportToWkt())
    
    band_out = ds_out.GetRasterBand(1)
    band_out.SetNoDataValue(-9999.0)
    out_grid[np.isnan(out_grid)] = -9999.0
    band_out.WriteArray(out_grid)
    ds_out.FlushCache()
    ds_out = None
    
    print(f"[GeoCode] Success: {output_tif} created.")

if __name__ == '__main__':
    main()
