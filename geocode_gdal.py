#!/usr/bin/env python3
"""雷达坐标网格 -> 经纬度 GeoTIFF (EPSG:4326), 用 trans.dat 做查找表, 最近邻取值。

用法: geocode_gdal.py trans.dat input_ra.grd output_ll.tif [res_deg] [full_ncols full_nrows]

  trans.dat       SAT_llt2rat -bod 的输出, 每条记录 5 个 float64: 距离向像元, 方位向像元, 高程, 经度, 纬度
                  (只含落在影像范围内的 DEM 点)
  full_ncols/rows 整景 SLC 的列数/行数 (PRM 的 num_rng_bins / num_lines)。输入网格比整景小 (降采样过) 时, 靠它把
                  trans.dat 的坐标缩放到输入网格上。不给就按 trans.dat 坐标的最大值推断: DEM 没盖满整景时
                  推断值偏小, 缩放就有误差 (降采样 16 倍时远端可偏几个像元)。
"""
import sys

import numpy as np
from osgeo import gdal, osr
from scipy.spatial import cKDTree


def geocode(trans, val_arr, res, full_dims=None, max_dist_factor=3.0, chunk_rows=500):
    """返回 (grid, lon_min, lat_max); grid 是 float32 [rows, cols], 无数据处为 NaN。

    距离最近的有效 DEM 点超过 max_dist_factor 倍点间距的像元置 NaN: 否则最近邻会把影像边缘的值
    一直外推到经纬度外包框的四个角, 在影像范围之外画出一片假数据。
    """
    nrows, ncols = val_arr.shape
    lon = trans[:, 3]
    lat = trans[:, 4]

    if full_dims is not None:
        full_ncols, full_nrows = full_dims
    else:
        full_ncols = int(np.ceil(trans[:, 0].max())) + 1
        full_nrows = int(np.ceil(trans[:, 1].max())) + 1

    scale_x = ncols / full_ncols
    scale_y = nrows / full_nrows
    if (ncols, nrows) != (full_ncols, full_nrows):
        print(f"[GeoCode] Input grid ({ncols}x{nrows}) differs from full size ({full_ncols}x{full_nrows}). "
              f"Scaling trans.dat coordinates by x:{scale_x:.6f}, y:{scale_y:.6f}")
    # SAT_llt2rat 的坐标是像元编号 (第 k 个像元的中心在 k), 所以包含该点的像元是 floor(坐标 + 0.5)
    rng = np.floor((trans[:, 0] + 0.5) * scale_x).astype(int)
    azi = np.floor((trans[:, 1] + 0.5) * scale_y).astype(int)

    mask = (rng >= 0) & (rng < ncols) & (azi >= 0) & (azi < nrows)
    lon, lat = lon[mask], lat[mask]
    rng, azi = rng[mask], azi[mask]

    vals = val_arr[azi, rng]
    valid = ~np.isnan(vals)
    lon, lat, vals = lon[valid], lat[valid], vals[valid]
    if len(vals) == 0:
        raise ValueError("No valid data points found.")

    lon_min, lon_max = lon.min(), lon.max()
    lat_min, lat_max = lat.min(), lat.max()
    out_cols = int(np.ceil((lon_max - lon_min) / res)) + 1
    out_rows = int(np.ceil((lat_max - lat_min) / res)) + 1
    print(f"[GeoCode] Target Grid Size: {out_cols} x {out_rows} "
          f"(Lon: {lon_min:.4f}~{lon_max:.4f}, Lat: {lat_min:.4f}~{lat_max:.4f})")

    points = np.column_stack((lon, lat))
    tree = cKDTree(points)
    sample = points[:: max(1, len(points) // 20000)]
    spacing = float(np.median(tree.query(sample, k=2)[0][:, 1])) if len(points) > 1 else res
    max_dist = max_dist_factor * max(spacing, res)

    # 像元中心: 与 GeoTransform (原点在左上角, 像元大小 res) 严格对应
    gx = lon_min + (np.arange(out_cols) + 0.5) * res
    gy = lat_max - (np.arange(out_rows) + 0.5) * res

    print(f"\tNearest-neighbour fill (DEM point spacing {spacing:.6f} deg, max distance {max_dist:.6f} deg)...")
    out = np.full((out_rows, out_cols), np.nan, dtype=np.float32)
    for r0 in range(0, out_rows, chunk_rows):
        r1 = min(out_rows, r0 + chunk_rows)
        gxx, gyy = np.meshgrid(gx, gy[r0:r1])
        dist, idx = tree.query(np.column_stack((gxx.ravel(), gyy.ravel())), distance_upper_bound=max_dist)
        ok = np.isfinite(dist)
        block = np.full(dist.shape, np.nan, dtype=np.float32)
        block[ok] = vals[idx[ok]]
        out[r0:r1] = block.reshape(r1 - r0, out_cols)
    return out, lon_min, lat_max


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 geocode_gdal.py <trans.dat> <input_ra.grd> <output_ll.tif> [res_deg] [full_ncols full_nrows]")
        sys.exit(1)

    trans_file = sys.argv[1]
    input_grd = sys.argv[2]
    output_tif = sys.argv[3]
    res = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0001
    full_dims = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else None

    print(f"[GeoCode] Reading {input_grd} and mapping via {trans_file} (Res: {res})")

    ds_in = gdal.Open(input_grd)
    if not ds_in:
        print(f"Error: Cannot open {input_grd}")
        sys.exit(1)
    val_arr = ds_in.GetRasterBand(1).ReadAsArray()
    ds_in = None

    trans = np.fromfile(trans_file, dtype=np.float64).reshape(-1, 5)

    try:
        out_grid, lon_min, lat_max = geocode(trans, val_arr, res, full_dims)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(output_tif, out_grid.shape[1], out_grid.shape[0], 1, gdal.GDT_Float32)
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
