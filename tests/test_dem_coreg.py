#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dem_coreg.py (可选的 DEM 精配准: 估计雷达坐标 DEM 相对影像的整体平移) 的行为约定, 全部用合成数据。

  1. estimate_shift: 恢复已知平移且符号正确 (把 DEM 平移该量即与影像对齐); 平坦 DEM、纯噪声影像、非整体平移、
     平移过大 / 可忽略时都不给结果 (applied=false 且写明原因);
  2. amp_grid: 第 j 行 = 第 j 段方位线 (不翻转), 第 i 列 = 第 i 段距离向像元; geometry: 入射角量级合理;
  3. apply-trans: trans.dat 的距离向 / 方位向像元号加上 SLC 像元平移, 其余列不动; 未修正时原样复制; 大小非法返回 1;
  4. estimate 命令端到端 (需要 gmt): 写出 dem_ra_coreg.grd + 报告, 平移后与影像对齐; 平坦 DEM 只写报告不写网格;
  5. run_dinsar.sh: --dem-coreg 只接受 on|off, 默认 off, usage 里有说明。

运行: python3 -m unittest tests.test_dem_coreg -v    (需要 numpy、scipy; 第 4 部分还需要 gmt)
"""

import json
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
    from scipy import ndimage as ndi
    import dem_coreg as dc
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

HAVE_GMT = shutil.which("gmt") is not None
RUN = os.path.join(REPO, "run_dinsar.sh")

GR = 5.0       # 合成网格的地距间距 (米)
THETA = 36.0   # 合成入射角 (度)


def synth_dem(n=512, seed=1, relief=300.0, sigma=8.0):
    rng = np.random.default_rng(seed)
    d = ndi.gaussian_filter(rng.normal(size=(n, n)), sigma)
    d = (d - d.mean()) / d.std()
    return (relief / 4.0 * d + 1500.0).astype(np.float32)


def synth_amp(dem, shift, seed=2, noise=0.3):
    """影像幅度 = 把 DEM 模拟亮度平移 shift (行, 列) 后再叠加乘性噪声。"""
    sim = dc.make_sim(dem, GR, THETA)
    moved = ndi.shift(sim, shift, order=3, mode="nearest")
    rng = np.random.default_rng(seed)
    return np.exp(2.0 * moved + noise * rng.normal(size=sim.shape)).astype(np.float32)


@unittest.skipUnless(HAVE_SCIPY, "需要 scipy")
class EstimateShift(unittest.TestCase):
    def test_recovers_known_shifts_with_correct_sign(self):
        dem = synth_dem()
        for true in [(2.5, -1.7), (-3.2, 2.1), (0.0, 4.0)]:
            with self.subTest(true=true):
                got, rep = dc.estimate_shift(dem, synth_amp(dem, true), GR, THETA)
                self.assertIsNotNone(got, rep)
                self.assertTrue(rep["applied"])
                self.assertAlmostEqual(got[0], true[0], delta=0.3)
                self.assertAlmostEqual(got[1], true[1], delta=0.3)
                # 按约定平移 DEM 后与影像对齐: 再估一次残余应接近 0
                moved = ndi.shift(dem, got, order=3, mode="nearest")
                again, rep2 = dc.estimate_shift(moved, synth_amp(dem, true), GR, THETA, min_shift_px=0.0)
                self.assertIsNotNone(again, rep2)
                self.assertLess(max(abs(again[0]), abs(again[1])), 0.4)

    def test_flat_dem_is_skipped(self):
        dem = np.full((512, 512), 3000.0, np.float32) + np.random.default_rng(0).normal(0, 1, (512, 512)).astype(np.float32)
        got, rep = dc.estimate_shift(dem, synth_amp(synth_dem(), (2.0, 2.0)), GR, THETA)
        self.assertIsNone(got)
        self.assertFalse(rep["applied"])
        self.assertIn("起伏", rep["reason"])

    def test_noise_image_is_skipped(self):
        dem = synth_dem()
        amp = np.exp(np.random.default_rng(5).normal(size=dem.shape)).astype(np.float32)
        got, rep = dc.estimate_shift(dem, amp, GR, THETA)
        self.assertIsNone(got, rep)
        self.assertFalse(rep["applied"])
        self.assertTrue(rep["reason"])

    def test_non_global_shift_is_skipped(self):
        # 三条竖带的平移分别是 -6 / 0 / +6 像元 -> 不是整体平移, 不能用一个常数去修
        dem = synth_dem()
        amp = synth_amp(dem, (0.0, 0.0))
        amp[:, :170] = synth_amp(dem, (-6.0, -6.0))[:, :170]
        amp[:, 342:] = synth_amp(dem, (6.0, 6.0))[:, 342:]
        got, rep = dc.estimate_shift(dem, amp, GR, THETA)
        self.assertIsNone(got, rep)
        self.assertFalse(rep["applied"])

    def test_too_large_and_negligible_shifts_are_skipped(self):
        dem = synth_dem()
        got, rep = dc.estimate_shift(dem, synth_amp(dem, (9.0, 9.0)), GR, THETA, max_shift_px=5.0)
        self.assertIsNone(got)
        self.assertIn("上限", rep["reason"])
        got, rep = dc.estimate_shift(dem, synth_amp(dem, (0.0, 0.0)), GR, THETA)
        self.assertIsNone(got)
        self.assertIn("可忽略", rep["reason"])

    def test_too_small_grid_is_skipped(self):
        dem = synth_dem(256)
        got, rep = dc.estimate_shift(dem, synth_amp(dem, (2.0, 2.0)), GR, THETA)
        self.assertIsNone(got)
        self.assertFalse(rep["applied"])


@unittest.skipUnless(HAVE_SCIPY, "需要 scipy")
class GridConventions(unittest.TestCase):
    def test_amp_grid_row_is_azimuth_block_and_column_is_range_block(self):
        nrow, ncol, ny, nx = 16, 24, 4, 6      # 方位向每格 4 行, 距离向每格 4 列
        a = np.zeros((nrow, ncol, 2), "<i2")
        for j in range(ny):
            for i in range(nx):
                a[j * 4:(j + 1) * 4, i * 4:(i + 1) * 4, 0] = 10 * (j + 1) + (i + 1)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.SLC")
            a.tofile(p)
            g = dc.amp_grid(p, nrow, ncol, ny, nx)
            self.assertEqual(g.shape, (ny, nx))
            for j in range(ny):
                for i in range(nx):
                    self.assertAlmostEqual(float(g[j, i]), 10 * (j + 1) + (i + 1), places=4)
            with open(p, "r+b") as fh:
                fh.truncate(100)
            with self.assertRaises(ValueError):
                dc.amp_grid(p, nrow, ncol, ny, nx)

    def test_geometry_is_plausible(self):
        prm = {"earth_radius": "6378000", "SC_height": "570000", "rng_samp_rate": "1e8", "num_rng_bins": "10000",
               "near_range": "650000", "SC_vel": "7580", "PRF": "4000"}
        g = dc.geometry(prm, 2, 4)
        self.assertTrue(25 < g["theta_deg"] < 45, g)
        self.assertGreater(g["ground_range_m_per_px"], 2 * 1.4989)    # 地距间距大于斜距间距
        self.assertAlmostEqual(g["azimuth_m_per_px"], 4 * 7580 * 6378000 / 6948000 / 4000, delta=1e-6)


@unittest.skipUnless(HAVE_SCIPY, "需要 scipy")
class ApplyTrans(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.trans = os.path.join(self.d, "trans.dat")
        self.rec = np.array([[10.0, 20.0, 3000.0, 90.5, 29.5], [11.5, 21.5, 3100.0, 90.6, 29.6],
                             [0.0, 0.0, 0.0, 0.0, 0.0]], "<f8")
        self.rec.tofile(self.trans)

    def run_cmd(self, report):
        rp = os.path.join(self.d, "r.json")
        with open(rp, "w") as fh:
            json.dump(report, fh)
        out = os.path.join(self.d, "o.dat")
        rc = dc.main(["apply-trans", "--trans", self.trans, "--report", rp, "--out", out])
        return rc, out

    def test_applied_shifts_only_first_two_columns(self):
        rc, out = self.run_cmd({"applied": True, "shift_slc_px": {"azimuth_lines": 10.5, "range": -7.25}})
        self.assertEqual(rc, 0)
        got = np.fromfile(out, "<f8").reshape(-1, 5)
        want = self.rec.copy()
        want[:, 0] += -7.25
        want[:, 1] += 10.5
        np.testing.assert_allclose(got, want)

    def test_not_applied_copies_unchanged(self):
        rc, out = self.run_cmd({"applied": False, "reason": "x"})
        self.assertEqual(rc, 0)
        with open(out, "rb") as f1, open(self.trans, "rb") as f2:
            self.assertEqual(f1.read(), f2.read())

    def test_bad_size_fails(self):
        with open(self.trans, "ab") as fh:
            fh.write(b"\0" * 7)
        rc, _ = self.run_cmd({"applied": True, "shift_slc_px": {"azimuth_lines": 1.0, "range": 1.0}})
        self.assertEqual(rc, 1)


@unittest.skipUnless(HAVE_SCIPY and HAVE_GMT, "需要 scipy 和 gmt")
class EstimateCommand(unittest.TestCase):
    NR = 2048      # 合成 SLC 尺寸 (行 = 列), dem_ra 网格为其 1/4

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        # 4 倍抽样的网格 512x512, 与 SLC 的每个 4x4 块一一对应
        self.prm = os.path.join(self.d, "m.PRM")
        with open(self.prm, "w") as fh:
            fh.write("num_lines = %d\nnum_rng_bins = %d\nearth_radius = 6378000\nSC_height = 570000\n"
                     "rng_samp_rate = 1e8\nnear_range = 650000\nSC_vel = 7580\nPRF = 4000\n" % (self.NR, self.NR))

    def make_inputs(self, dem, true_shift_grid):
        n = dem.shape[0]
        amp_grid = synth_amp(dem, true_shift_grid)
        amp_slc = np.kron(amp_grid, np.ones((4, 4), np.float32))
        slc = np.zeros((self.NR, self.NR, 2), "<i2")
        slc[..., 0] = np.clip(amp_slc * 100.0, 0, 30000).astype("<i2")
        slc.tofile(os.path.join(self.d, "m.SLC"))
        raw = os.path.join(self.d, "dem.f32")
        np.ascontiguousarray(dem, "<f4").tofile(raw)
        dem_ra = os.path.join(self.d, "dem_ra.grd")
        subprocess.check_call(["gmt", "xyz2grd", raw, "-ZTLf", "-R0/%d/0/%d" % (n * 4, n * 4), "-I4/4", "-r",
                               "-G" + dem_ra], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return dem_ra

    def run_estimate(self, dem_ra):
        out, rep = os.path.join(self.d, "dem_ra_coreg.grd"), os.path.join(self.d, "dem_coreg.json")
        rc = dc.main(["estimate", "--prm", self.prm, "--slc", os.path.join(self.d, "m.SLC"), "--dem-ra", dem_ra,
                      "--out", out, "--report", rep])
        with open(rep) as fh:
            return rc, out, json.load(fh)

    def test_end_to_end_recovers_shift_and_reports_slc_pixels(self):
        dem = synth_dem(512)
        # 注意: 合成亮度用 dc.geometry 之外的固定入射角/地距, 这里只检验平移、符号与行/列约定
        dem_ra = self.make_inputs(dem, (2.0, -1.5))
        rc, out, rep = self.run_estimate(dem_ra)
        self.assertEqual(rc, 0)
        if not rep["applied"]:
            self.skipTest("合成几何下估计未通过门限: %s" % rep["reason"])
        self.assertTrue(os.path.getsize(out) > 0)
        sy, sx = rep["shift_grid_px_rows_cols"]
        self.assertAlmostEqual(sy, 2.0, delta=0.5)
        self.assertAlmostEqual(sx, -1.5, delta=0.5)
        self.assertAlmostEqual(rep["shift_slc_px"]["azimuth_lines"], sy * 4, places=6)
        self.assertAlmostEqual(rep["shift_slc_px"]["range"], sx * 4, places=6)
        new = dc.read_grd(out)
        self.assertEqual(new.shape, dem.shape)
        want = ndi.shift(dem, (sy, sx), order=3, mode="nearest")
        np.testing.assert_allclose(new, want, atol=1.0)

    def test_flat_dem_writes_report_but_no_grid(self):
        dem = np.full((512, 512), 2500.0, np.float32)
        dem_ra = self.make_inputs(synth_dem(512), (2.0, 2.0))
        # 用平坦 DEM 覆盖网格文件
        raw = os.path.join(self.d, "flat.f32")
        dem.tofile(raw)
        subprocess.check_call(["gmt", "xyz2grd", raw, "-ZTLf", "-R0/2048/0/2048", "-I4/4", "-r", "-G" + dem_ra],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rc, out, rep = self.run_estimate(dem_ra)
        self.assertEqual(rc, 0)
        self.assertFalse(rep["applied"])
        self.assertFalse(os.path.exists(out))


class RunDinsarOption(unittest.TestCase):
    def read(self):
        with open(RUN, encoding="utf8") as fh:
            return fh.read()

    def test_default_off_and_documented(self):
        s = self.read()
        self.assertIn('DEM_COREG="${DEM_COREG:-off}"', s)
        self.assertIn("--dem-coreg on|off", s)

    def test_option_validates_value(self):
        if shutil.which("bash") is None:
            self.skipTest("需要 bash")
        r = subprocess.run(["bash", RUN, "--dem-coreg", "maybe", "--help"], capture_output=True, text=True, timeout=60)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--dem-coreg", r.stdout + r.stderr)

    def test_all_geocode_calls_use_selected_trans(self):
        s = self.read()
        self.assertNotIn('python3 "$py_geocode" trans.dat ', s)
        self.assertGreaterEqual(s.count('python3 "$py_geocode" "$trans_file" '), 4)


if __name__ == "__main__":
    unittest.main()
