#!/usr/bin/env python3
"""把 SAT_llt2rat 给出的 DEM 点 (雷达坐标 + 高程) 内插成 phasediff 用的 topo_ra 网格。

用法: topo_ra_grid.py dem.xyz num_rng num_az dec out.f32

  dem.xyz           SAT_llt2rat -bod 的输出, 每条记录 5 个 float64: 距离向像元, 方位向像元, 高程, 经度, 纬度。
                    SAT_llt2rat 只输出落在影像外扩 10/20 像元以内的点, 所以影像边缘处也有点可用。
  num_rng, num_az   整景 SLC 的列数/行数 (PRM 的 num_rng_bins / num_lines), 都要能被 dec 整除
  dec               输出网格相对 SLC 的降采样倍数。输出的第 (j, i) 个值对应 SLC 的
                    第 j*dec ~ (j+1)*dec-1 行、第 i*dec ~ (i+1)*dec-1 列, 取该块中心处的内插值。
                    phasediff 读 topo 时就是这样用的: 第 k 个像元取第 k/dec (整除) 个节点, 不再内插。
  out.f32           float32 原始二进制, (num_az/dec) 行 x (num_rng/dec) 列, 第 j 行就是第 j 段方位线, 没有上下翻转。
                    phasediff 把 GMT 网格内存里的第 j 行当作第 j 条方位线来读, 所以转 grd 时用
                        gmt xyz2grd out.f32 -R0/num_rng/0/num_az -I<dec> -r -ZTLf -Gdem_ra.grd
                    (-ZTL: 文件第一行放在网格最上面一行, 正好是 GMT 内存里的第 0 行)。

为什么不用 gmt triangulate: 这里的 GMT 用的是 Watson 算法, 电建一号 SM 整景的 105 万个 DEM 点跑了 15 分钟还没完。
scipy 的 Qhull Delaunay 处理同样的点只要约 9 秒, 再内插到 2270 万个节点约 9 秒, 峰值内存约 0.8 GB。
两者都是对同一批点做 Delaunay 三角网上的线性内插, 结果一致 (见 tests/test_topo_ra_grid.py)。
点云是整体送进 Qhull 的, 每百万点约占 0.8 GB 内存。
"""
import sys

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree

RECORD_DOUBLES = 5


def load_points(path):
    """读 SAT_llt2rat -bod 的输出, 返回 [n, 3] (距离向, 方位向, 高程), 丢掉含 NaN/inf 的点。"""
    raw = np.fromfile(path, dtype="<f8")
    if raw.size % RECORD_DOUBLES:
        raise ValueError(f"{path} 的大小不是 {RECORD_DOUBLES} 个 float64 的整数倍 "
                         f"(SAT_llt2rat -bod 每条记录 {RECORD_DOUBLES} 个 float64), 文件可能不完整")
    pts = raw.reshape(-1, RECORD_DOUBLES)[:, :3]
    return pts[np.isfinite(pts).all(axis=1)]


def block_mean(pts, size):
    """每个 size x size 的方格只留一个点 (方格内各点的距离向/方位向/高程的平均)。

    只合并同一方格里的点: 点比方格稀疏的地方不受影响; DEM 比输出网格密得多时 (如 10 m DEM) 靠它把点数压到不多于节点数。
    """
    ix = np.floor(pts[:, 0] / size).astype(np.int64)
    iy = np.floor(pts[:, 1] / size).astype(np.int64)
    ix -= ix.min()
    iy -= iy.min()
    key = iy * (ix.max() + 1) + ix
    order = np.argsort(key, kind="stable")
    key = key[order]
    first = np.flatnonzero(np.concatenate(([True], key[1:] != key[:-1])))
    if len(first) == len(pts):
        return pts
    counts = np.diff(np.append(first, len(pts)))
    return np.add.reduceat(pts[order], first, axis=0) / counts[:, None]


def grid_topo(pts, num_rng, num_az, dec, chunk_nodes=2_000_000):
    """把 [n, 3] 的 (距离向, 方位向, 高程) 点内插到 (num_az/dec) x (num_rng/dec) 个节点上, 返回 float32 [行=方位向, 列=距离向]。

    节点取在 dec x dec 像元块的中心 (像元 k 的中心在距离向坐标 k, 所以块中心是 i*dec + (dec-1)/2)。
    Delaunay 三角网的凸包以外 (DEM 没盖到的地方) 用最近的 DEM 点填充, 并返回填充的节点数。
    """
    if num_rng % dec or num_az % dec:
        raise ValueError(f"SLC 尺寸 {num_rng}x{num_az} 不能被 dec={dec} 整除")

    pts = block_mean(pts, float(dec))
    if len(pts) < 3:
        raise ValueError(f"落在影像范围内的 DEM 点只有 {len(pts)} 个, 无法内插 (DEM 与影像没有重叠?)")

    try:
        tri = Delaunay(pts[:, :2])
    except Exception as e:
        raise ValueError(f"DEM 点无法做三角剖分: {e}") from e
    interp = LinearNDInterpolator(tri, pts[:, 2], fill_value=np.nan)

    nx, ny = num_rng // dec, num_az // dec
    xs = np.arange(nx) * float(dec) + (dec - 1) / 2.0
    ys = np.arange(ny) * float(dec) + (dec - 1) / 2.0
    out = np.empty((ny, nx), dtype=np.float32)
    rows = max(1, chunk_nodes // nx)
    for j0 in range(0, ny, rows):
        j1 = min(ny, j0 + rows)
        xx, yy = np.meshgrid(xs, ys[j0:j1])
        out[j0:j1] = interp(np.column_stack((xx.ravel(), yy.ravel()))).reshape(j1 - j0, nx)

    bad = np.isnan(out)
    n_filled = int(bad.sum())
    if n_filled:
        jj, ii = np.nonzero(bad)
        _, idx = cKDTree(pts[:, :2]).query(np.column_stack((xs[ii], ys[jj])))
        out[bad] = pts[idx, 2]
    return out, n_filled


def main(argv):
    if len(argv) != 6:
        print(__doc__)
        return 1
    dem_xyz, num_rng, num_az, dec, out_path = argv[1], int(argv[2]), int(argv[3]), int(argv[4]), argv[5]
    if dec < 1:
        print(f"Error: dec 必须 >= 1, 收到 {dec}", file=sys.stderr)
        return 1

    try:
        pts = load_points(dem_xyz)
        print(f"[topo_ra_grid] {len(pts)} 个 DEM 点, 输出 {num_rng // dec} x {num_az // dec} 节点 (dec={dec})")
        grid, n_filled = grid_topo(pts, num_rng, num_az, dec)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2

    if n_filled:
        frac = n_filled / grid.size
        print(f"[topo_ra_grid] {n_filled} 个节点 ({frac:.2%}) 在 DEM 覆盖范围以外, 已用最近的 DEM 点填充"
              + (" -- DEM 没有盖满整景, 请检查范围" if frac > 0.01 else ""))
    grid.tofile(out_path)
    print(f"[topo_ra_grid] 高程范围 {grid.min():.1f} ~ {grid.max():.1f} m, 已写入 {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
