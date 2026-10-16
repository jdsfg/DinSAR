#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DEM -> 雷达坐标 topo_ra 网格 (topo_ra_grid.py / fast_dem2topo_ra.sh) 与地理编码 (geocode_gdal.py) 的尺寸、朝向约定。

  1. topo_ra_grid.py: 平面被精确复现; 行号 = 方位向块号 (不翻转); 点比节点密/稀、DEM 没盖满整景时都出得了结果;
     非法输入给明确的错误; 与 gmt triangulate 在同一批点上结果一致 (没有 gmt 则跳过);
  2. fast_dem2topo_ra.sh (假的 SAT_llt2rat + 真实 gmt): 按 SAT_llt2rat 的 5 列记录读点, 产出 phasediff 要的网格
     (列数/行数 = SLC 尺寸/dec, GMT 内存里的第 j 行 = 第 j 段方位线), dec=auto 的取值, 不能整除时退回 1,
     失败时返回非 0 且不留下 dem_ra.grd;
  3. geocode_gdal.py: 给了整景尺寸就按它缩放, 不再用 DEM 点坐标的最大值推断; 像元编号按"中心在整数"取整。

运行: python3 -m unittest tests.test_topo_ra_grid -v    (需要 numpy、scipy; 第 2 部分还需要 bash、gmt, 第 3 部分需要 osgeo)
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)

try:
    import topo_ra_grid as tg
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
try:
    import geocode_gdal as gg
    HAVE_GDAL = True
except ImportError:
    HAVE_GDAL = False

HAVE_GMT = shutil.which("gmt") is not None
HELPER = os.path.join(REPO, "fast_dem2topo_ra.sh")
GRID_SCRIPT = os.path.join(REPO, "topo_ra_grid.py")


def jittered_lattice(nr, na, step, jitter, margin=12.0, seed=0):
    """覆盖 [-margin, nr+margin] x [-margin, na+margin] 的、带随机抖动的格点, 返回 (r, a)。"""
    rng = np.random.default_rng(seed)
    r = np.arange(-margin, nr + margin + step, step, dtype=float)
    a = np.arange(-margin, na + margin + step, step, dtype=float)
    rr, aa = np.meshgrid(r, a)
    rr = rr + rng.uniform(-jitter, jitter, rr.shape)
    aa = aa + rng.uniform(-jitter, jitter, aa.shape)
    return rr.ravel(), aa.ravel()


def block_centers(n, dec):
    return np.arange(n // dec) * dec + (dec - 1) / 2.0


@unittest.skipUnless(HAVE_SCIPY, "需要 scipy")
class GridTopo(unittest.TestCase):
    NR, NA, DEC = 64, 48, 4

    def _plane(self, r, a):
        return 100.0 + 0.5 * r - 0.25 * a

    def test_plane_is_reproduced_at_block_centres(self):
        r, a = jittered_lattice(self.NR, self.NA, step=6.0, jitter=1.5)
        out, n_filled = tg.grid_topo(np.column_stack((r, a, self._plane(r, a))), self.NR, self.NA, self.DEC)
        self.assertEqual(out.shape, (self.NA // self.DEC, self.NR // self.DEC))
        self.assertEqual(out.dtype, np.float32)
        self.assertEqual(n_filled, 0)
        rc, ac = np.meshgrid(block_centers(self.NR, self.DEC), block_centers(self.NA, self.DEC))
        np.testing.assert_allclose(out, self._plane(rc, ac), atol=2e-3)

    def test_rows_are_azimuth_and_columns_are_range(self):
        """第 j 行 = 第 j 段方位线 (不翻转), 第 i 列 = 第 i 段距离向: 用只随方位向/只随距离向变化的面分别验证。"""
        r, a = jittered_lattice(self.NR, self.NA, step=6.0, jitter=1.5)
        out_a, _ = tg.grid_topo(np.column_stack((r, a, 10.0 * a)), self.NR, self.NA, self.DEC)
        np.testing.assert_allclose(out_a[:, 0], 10.0 * block_centers(self.NA, self.DEC), atol=2e-3)
        np.testing.assert_allclose(out_a[:, -1], out_a[:, 0], atol=2e-3)
        out_r, _ = tg.grid_topo(np.column_stack((r, a, 10.0 * r)), self.NR, self.NA, self.DEC)
        np.testing.assert_allclose(out_r[0, :], 10.0 * block_centers(self.NR, self.DEC), atol=2e-3)
        np.testing.assert_allclose(out_r[-1, :], out_r[0, :], atol=2e-3)

    def test_dem_denser_than_the_grid(self):
        """DEM 比节点密得多时, 同一方格里的点被合并 (点数不超过节点数), 平面照样复现。"""
        r, a = jittered_lattice(self.NR, self.NA, step=0.5, jitter=0.2)
        pts = np.column_stack((r, a, self._plane(r, a)))
        self.assertGreater(len(pts), 4 * (self.NR // self.DEC) * (self.NA // self.DEC))
        merged = tg.block_mean(pts, float(self.DEC))
        self.assertLess(len(merged), len(pts) / 10)
        cells = {(int(np.floor(x / self.DEC)), int(np.floor(y / self.DEC))) for x, y in merged[:, :2]}
        self.assertEqual(len(cells), len(merged), "合并后每个方格只应剩一个点")
        out, _ = tg.grid_topo(pts, self.NR, self.NA, self.DEC)
        rc, ac = np.meshgrid(block_centers(self.NR, self.DEC), block_centers(self.NA, self.DEC))
        np.testing.assert_allclose(out, self._plane(rc, ac), atol=2e-3)

    def test_block_mean_only_merges_points_in_the_same_cell(self):
        pts = np.array([[0.5, 0.5, 10.0], [1.5, 2.5, 20.0], [3.9, 3.9, 30.0],     # 同在 [0,4) x [0,4)
                        [4.5, 0.5, 99.0]])                                          # 另一个方格
        merged = tg.block_mean(pts, 4.0)
        self.assertEqual(len(merged), 2)
        np.testing.assert_allclose(merged[np.argmin(merged[:, 0])], pts[:3].mean(axis=0))
        np.testing.assert_allclose(merged[np.argmax(merged[:, 0])], pts[3])
        apart = pts[[0, 3]]                                   # 没有两个点在同一个方格里: 原样返回
        self.assertIs(tg.block_mean(apart, 4.0), apart)

    def test_uncovered_nodes_take_the_nearest_dem_value(self):
        """DEM 只盖住左半幅: 其余节点不留 NaN, 用最近的 DEM 点填充, 并报告填了多少个。"""
        r, a = jittered_lattice(self.NR, self.NA, step=6.0, jitter=1.5)
        keep = r < 30.0
        pts = np.column_stack((r[keep], a[keep], self._plane(r[keep], a[keep])))
        out, n_filled = tg.grid_topo(pts, self.NR, self.NA, self.DEC)
        self.assertTrue(np.isfinite(out).all())
        self.assertGreater(n_filled, 0)
        covered = out[:, : 24 // self.DEC]
        rc, ac = np.meshgrid(block_centers(self.NR, self.DEC), block_centers(self.NA, self.DEC))
        np.testing.assert_allclose(covered, self._plane(rc, ac)[:, : 24 // self.DEC], atol=2e-3)

    def test_bad_input_is_reported(self):
        pts = np.array([[1.0, 1.0, 5.0], [2.0, 9.0, 6.0], [8.0, 3.0, 7.0]])
        with self.assertRaisesRegex(ValueError, "整除"):
            tg.grid_topo(pts, 64, 50, 4)
        with self.assertRaisesRegex(ValueError, "无法内插"):
            tg.grid_topo(pts[:2], 64, 48, 4)
        with self.assertRaisesRegex(ValueError, "三角剖分"):
            tg.grid_topo(np.array([[0.0, 0.0, 1.0], [5.0, 5.0, 2.0], [9.0, 9.0, 3.0]]), 64, 48, 4)   # 三点共线

    def test_load_points_keeps_three_columns_and_drops_nonfinite(self):
        rec = np.array([[1.0, 2.0, 3.0, 109.1, 40.1],
                        [4.0, np.nan, 6.0, 109.2, 40.2],
                        [7.0, 8.0, 9.0, 109.3, 40.3]])
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "dem.xyz")
            rec.astype("<f8").tofile(p)
            np.testing.assert_array_equal(tg.load_points(p), rec[[0, 2], :3])
            with open(p, "ab") as f:
                f.write(b"\0" * 8)
            with self.assertRaisesRegex(ValueError, "不完整"):
                tg.load_points(p)

    @unittest.skipUnless(HAVE_GMT, "需要 gmt")
    def test_matches_gmt_triangulate(self):
        """与 gmt triangulate 在同一批点上逐节点一致 (两者都是 Delaunay 三角网上的线性内插)。

        GMT 的 -r 网格把像元 k 的中心放在 k+0.5, 本脚本放在 k (SAT_llt2rat 的像元编号); 所以给 GMT 的点整体平移 +0.5。
        点间距比方格大, 不会触发方格合并, 两边的点集完全相同。"""
        nr, na, dec = 96, 64, 4
        r, a = jittered_lattice(nr, na, step=8.0, jitter=1.5, seed=3)
        z = 50.0 + 20.0 * np.sin(r / 13.0) + 15.0 * np.cos(a / 9.0)
        mine, _ = tg.grid_topo(np.column_stack((r, a, z)), nr, na, dec)
        with tempfile.TemporaryDirectory() as d:
            pts = os.path.join(d, "pts.bin")
            np.column_stack((r + 0.5, a + 0.5, z)).astype("<f8").tofile(pts)
            grd = os.path.join(d, "tri.grd")
            subprocess.run(["gmt", "triangulate", pts, "-R0/%d/0/%d" % (nr, na), "-I%d/%d" % (dec, dec),
                            "-bi3d", "-G" + grd, "-r"], cwd=d, check=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
            raw = subprocess.run(["gmt", "grd2xyz", grd, "-ZTLf"], cwd=d, check=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).stdout
        # GMT 内存顺序是第一行 = y 最大; 翻成 y 升序后与 mine (第 j 行 = 方位向块 j) 对齐
        gmt_grid = np.frombuffer(raw, dtype="<f4").reshape(na // dec, nr // dec)[::-1]
        both = np.isfinite(gmt_grid) & np.isfinite(mine)
        self.assertGreater(both.mean(), 0.95)
        np.testing.assert_allclose(mine[both], gmt_grid[both], atol=1e-3)

    def test_command_line(self):
        r, a = jittered_lattice(self.NR, self.NA, step=6.0, jitter=1.5)
        rec = np.column_stack((r, a, self._plane(r, a), 109.0 + r / 1000.0, 40.0 + a / 1000.0))
        with tempfile.TemporaryDirectory() as d:
            src, dst = os.path.join(d, "dem.xyz"), os.path.join(d, "topo.f32")
            rec.astype("<f8").tofile(src)
            run = lambda *args: subprocess.run([sys.executable, GRID_SCRIPT] + [str(x) for x in args],
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            ok = run(src, self.NR, self.NA, self.DEC, dst)
            self.assertEqual(ok.returncode, 0, ok.stderr)
            grid = np.fromfile(dst, dtype="<f4").reshape(self.NA // self.DEC, self.NR // self.DEC)
            rc, ac = np.meshgrid(block_centers(self.NR, self.DEC), block_centers(self.NA, self.DEC))
            np.testing.assert_allclose(grid, self._plane(rc, ac), atol=2e-3)

            bad = run(src, self.NR, 50, self.DEC, dst + ".2")
            self.assertEqual(bad.returncode, 2)
            self.assertIn("整除", bad.stderr)
            self.assertFalse(os.path.exists(dst + ".2"))
            self.assertEqual(run(src, self.NR).returncode, 1)


FAKE_SAT_LLT2RAT = r'''#!/usr/bin/env python3
"""假的 SAT_llt2rat master.PRM 1 -bod: 文本 "经度 纬度 高程" 进, 5 个 float64 (r, a, 高程, 经度, 纬度) 出,
只输出落在影像外扩 10/20 像元以内的点 (与真实程序一致)。映射: r = (经度-109)*1000, a = (纬度-40)*1000。"""
import re
import struct
import sys

prm = open(sys.argv[1]).read()
nr = int(re.search(r"num_rng_bins\s*=\s*(\d+)", prm).group(1))
na = int(re.search(r"num_lines\s*=\s*(\d+)", prm).group(1))
out = sys.stdout.buffer
for line in sys.stdin:
    f = line.split()
    if len(f) < 3:
        continue
    lon, lat, z = float(f[0]), float(f[1]), float(f[2])
    r, a = (lon - 109.0) * 1000.0, (lat - 40.0) * 1000.0
    if -10.0 <= r <= nr + 10.0 and -20.0 <= a <= na + 20.0:
        out.write(struct.pack("<5d", r, a, z, lon, lat))
'''


@unittest.skipUnless(HAVE_SCIPY and HAVE_GMT and shutil.which("bash"), "需要 scipy、gmt、bash")
class DemHelper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="dem_helper_")
        cls.bin = os.path.join(cls.tmp, "bin")
        os.makedirs(cls.bin)
        p = os.path.join(cls.bin, "SAT_llt2rat")
        with open(p, "w") as f:
            f.write(FAKE_SAT_LLT2RAT)
        os.chmod(p, 0o755)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _workdir(self, name, nr, na, lon_shift=0.0):
        """工作目录: PRM (只含 SLC 尺寸) 和 dem.grd (z = 100 + 5*r + 3*a, 点距 0.5 像元, 覆盖整景外扩)。"""
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        with open(os.path.join(d, "m.PRM"), "w") as f:
            f.write("num_rng_bins = %d\nnum_lines = %d\n" % (nr, na))
        w = 109.0 - 0.004 + lon_shift
        e = 109.0 + (nr + 4) / 1000.0 + lon_shift
        s, n = 40.0 - 0.004, 40.0 + (na + 4) / 1000.0
        subprocess.run(["gmt", "grdmath", "-R%.6f/%.6f/%.6f/%.6f" % (w, e, s, n), "-I0.0005",
                        "X", "%.6f" % (109.0 + lon_shift), "SUB", "5000", "MUL",
                        "Y", "40", "SUB", "3000", "MUL", "ADD", "100", "ADD", "=", "dem.grd"],
                       cwd=d, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return d

    def _run(self, d, dec):
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ.get("PATH", ""))
        r = subprocess.run(["bash", HELPER, "m.PRM", "dem.grd", "log.txt", str(dec)], cwd=d, env=env,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=300)
        with open(os.path.join(d, "log.txt")) as f:
            return r.returncode, f.read()

    def _grid(self, d):
        """返回 (grdinfo -C 的字段, GMT 内存顺序的 [行, 列] 数据)。"""
        info = subprocess.run(["gmt", "grdinfo", "-C", "dem_ra.grd"], cwd=d, check=True, stdout=subprocess.PIPE,
                              universal_newlines=True).stdout.split()
        nx, ny = int(info[9]), int(info[10])
        raw = subprocess.run(["gmt", "grd2xyz", "dem_ra.grd", "-ZTLf"], cwd=d, check=True, stdout=subprocess.PIPE).stdout
        return info, np.frombuffer(raw, dtype="<f4").reshape(ny, nx)

    def _expected(self, nr, na, dec):
        rc, ac = np.meshgrid(block_centers(nr, dec), block_centers(na, dec))
        return 100.0 + 5.0 * rc + 3.0 * ac

    def test_grid_matches_phasediff_layout(self):
        d = self._workdir("layout", 24, 16)
        rc, log = self._run(d, 4)
        self.assertEqual(rc, 0, log)
        info, grid = self._grid(d)
        self.assertEqual([float(x) for x in info[1:5]], [0.0, 24.0, 0.0, 16.0])
        self.assertEqual([float(info[7]), float(info[8]), int(info[9]), int(info[10])], [4.0, 4.0, 6, 4])
        # 第 j 行是第 j 段方位线, 第 i 列是第 i 段距离向 (phasediff 读 topo 就是按 data[xt + ncols*yt] 这个顺序)
        np.testing.assert_allclose(grid, self._expected(24, 16, 4), atol=0.05)
        self.assertGreater(grid[-1, 0], grid[0, 0])
        for leftover in ("dem.xyz", "topo_ra.f32"):
            self.assertFalse(os.path.exists(os.path.join(d, leftover)), leftover)

    def test_auto_picks_largest_common_divisor_of_4_2_1(self):
        for name, nr, na, dec in (("auto4", 24, 16, 4), ("auto2", 26, 16, 2), ("auto1", 25, 16, 1)):
            d = self._workdir(name, nr, na)
            rc, log = self._run(d, "auto")
            self.assertEqual(rc, 0, log)
            info, grid = self._grid(d)
            self.assertEqual((grid.shape, float(info[7])), ((na // dec, nr // dec), float(dec)), name)
            np.testing.assert_allclose(grid, self._expected(nr, na, dec), atol=0.05)

    def test_dec_that_does_not_divide_falls_back_to_full_resolution(self):
        d = self._workdir("fallback", 24, 16)
        rc, log = self._run(d, 3)
        self.assertEqual(rc, 0, log)
        self.assertIn("不能被 3 整除", log)
        info, grid = self._grid(d)
        self.assertEqual(grid.shape, (16, 24))
        np.testing.assert_allclose(grid, self._expected(24, 16, 1), atol=0.05)

    def test_failure_returns_nonzero_and_leaves_no_grid(self):
        d = self._workdir("outside", 24, 16, lon_shift=1.0)    # DEM 整体在影像之外: 没有点落在影像范围内
        rc, log = self._run(d, 4)
        self.assertEqual(rc, 1, log)
        self.assertIn("DEM 投到雷达坐标失败", log)
        self.assertFalse(os.path.exists(os.path.join(d, "dem_ra.grd")))
        self.assertFalse(os.path.exists(os.path.join(d, "dem.xyz")))


@unittest.skipUnless(HAVE_SCIPY and HAVE_GDAL, "需要 scipy 与 osgeo")
class GeocodeDimensions(unittest.TestCase):
    """DEM 格点: r = (经度-109)*1000, a = (纬度-40)*1000, 点距 0.4 像元; 影像 64 x 48, DEM 只盖到 r<=59.6, a<=45.6。"""
    NR, NA, RES = 64, 48, 0.0001

    def _trans(self):
        rr, aa = np.meshgrid(np.arange(0.0, 59.7, 0.4), np.arange(0.0, 45.7, 0.4))
        r, a = rr.ravel(), aa.ravel()
        return np.column_stack((r, a, np.zeros_like(r), 109.0 + r / 1000.0, 40.0 + a / 1000.0))

    def _value_at(self, grid, lon_min, lat_max, r, a):
        lon, lat = 109.0 + r / 1000.0, 40.0 + a / 1000.0
        col = int(round((lon - lon_min) / self.RES - 0.5))
        row = int(round((lat_max - lat) / self.RES - 0.5))
        return grid[row, col]

    def test_full_resolution_input(self):
        val = np.add.outer(1000.0 * np.arange(self.NA), np.arange(self.NR))          # val[j, i] = 1000 j + i
        grid, lon_min, lat_max = gg.geocode(self._trans(), val, self.RES, (self.NR, self.NA))
        # 像元编号按"中心在整数"取整: r=57.6 属于第 58 个像元, 不是第 57 个
        self.assertEqual(self._value_at(grid, lon_min, lat_max, 57.6, 40.0), 1000.0 * 40 + 58)
        self.assertEqual(self._value_at(grid, lon_min, lat_max, 57.2, 39.6), 1000.0 * 40 + 57)

    def test_downsampled_input_uses_the_given_full_size_not_the_dem_extent(self):
        val = np.add.outer(100.0 * np.arange(self.NA // 4), np.arange(self.NR // 4))   # 16 x 12, val[j, i] = 100 j + i
        grid, lon_min, lat_max = gg.geocode(self._trans(), val, self.RES, (self.NR, self.NA))
        # 缩放比取 1/4: r=58 -> 第 floor(58.5/4)=14 列, a=40 -> 第 floor(40.5/4)=10 行
        self.assertEqual(self._value_at(grid, lon_min, lat_max, 58.0, 40.0), 100.0 * 10 + 14)
        # 不给整景尺寸时按 DEM 点的最大坐标推断 (60 x 46), 缩放比偏离 1/4, 远端会错一列
        inferred, lon_min2, lat_max2 = gg.geocode(self._trans(), val, self.RES)
        self.assertNotEqual(self._value_at(inferred, lon_min2, lat_max2, 58.0, 40.0), 100.0 * 10 + 14)

    def test_pixels_outside_the_dem_footprint_are_nodata(self):
        """影像足迹是斜的 (这里用一条对角带): 经纬度外包框的空角不能被最近邻外推填上假数据。"""
        val = np.ones((self.NA, self.NR))
        trans = self._trans()
        trans = trans[np.abs(trans[:, 0] - trans[:, 1]) < 5.0]
        grid, _, _ = gg.geocode(trans, val, self.RES, (self.NR, self.NA))
        self.assertTrue(np.isnan(grid[0, 0]))          # 左上角: 距离向最小、方位向最大, 不在带内
        self.assertTrue(np.isnan(grid[-1, -1]))        # 右下角
        self.assertEqual(grid[-1, 0], 1.0)             # 左下角: r=0, a=0, 在带内


if __name__ == "__main__":
    unittest.main()
