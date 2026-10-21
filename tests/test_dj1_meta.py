#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""电建一号新版产品: 元数据解析、干涉几何、像对筛选、h5 -> GMTSAR 转换器的单元测试。

全部使用合成的圆轨道和小尺寸 h5, 不依赖真实数据。几何量有解析真值:
  * 从影像轨道 = 主影像轨道整体平移一个垂直于速度的向量 b  =>  B⊥、B∥ 可精确算出;
  * 从影像只是把起始时间推后 d 秒                         =>  方位向偏移 = -d / 方位采样间隔。
h5 转换测试需要 h5py (没有则跳过)。

运行: python3 -m unittest tests.test_dj1_meta -v   或   python3 -m unittest discover -s tests -v
"""

import datetime as dt
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "gmtsar-master", "preproc", "DJ1_preproc", "src_swath")
sys.path.insert(0, SRC)
import dj1_meta as dm  # noqa: E402
import dj1_pair_select as sel  # noqa: E402
import make_slc_dj1h5 as conv  # noqa: E402

try:
    conv._import_h5py()
    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

HAVE_GMTSAR = shutil.which("calc_dop_orb") is not None and shutil.which("SAT_llt2rat") is not None

GM = 3.986004418e14
R_ORBIT = 6378137.0 + 531000.0
INCL = math.radians(97.5)
OMEGA = math.sqrt(GM / R_ORBIT ** 3)
T_REF = dm.to_sec(dt.datetime(2030, 1, 5, 4, 23, 45))   # 场景中心时刻


def _orbit_state(t, node_lon=0.0, tshift=0.0):
    """圆轨道 (惯性系近似为 ECEF, 测试不关心地球自转), 返回位置、速度。tshift 为轨道相位平移。"""
    u = OMEGA * (t - T_REF + tshift)
    p_orb = np.array([math.cos(u), math.sin(u), 0.0])
    v_orb = np.array([-math.sin(u), math.cos(u), 0.0])
    ci, si = math.cos(INCL), math.sin(INCL)
    rx = np.array([[1, 0, 0], [0, ci, -si], [0, si, ci]])
    cl, sl = math.cos(node_lon), math.sin(node_lon)
    rz = np.array([[cl, -sl, 0], [sl, cl, 0], [0, 0, 1]])
    m = rz @ rx
    return R_ORBIT * (m @ p_orb), R_ORBIT * OMEGA * (m @ v_orb)


def _fmt_time(sec):
    return dm.from_sec(sec).strftime("%Y-%m-%dT%H:%M:%S.%f")


def _zero_doppler_analytic(target, offset=(0.0, 0.0, 0.0)):
    """目标点在 (平移后的) 解析圆轨道上的零多普勒时刻: 对 (P(t)-T)·V(t) 二分求根。"""
    off = np.asarray(offset, dtype=float)

    def f(t):
        p, v = _orbit_state(t)
        return float(np.dot(p + off - target, v))

    lo, hi = T_REF - 20.0, T_REF + 20.0
    flo = f(lo)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def _center_target():
    """场景中心点: 取 T_REF 时刻卫星位置沿垂直于速度方向 (左侧下视) 前移 700 km 斜距的点。"""
    p, v = _orbit_state(T_REF)
    n = -p / np.linalg.norm(p)
    cross = np.cross(n, v / np.linalg.norm(v))
    d = math.cos(math.radians(37)) * n + math.sin(math.radians(37)) * cross
    return p + 700e3 * d


def _ecef_to_llh(x):
    lon = math.degrees(math.atan2(x[1], x[0]))
    p = math.hypot(x[0], x[1])
    lat = math.atan2(x[2], p * (1 - dm.WGS84_E2))
    for _ in range(8):
        n = dm.WGS84_A / math.sqrt(1 - dm.WGS84_E2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - n
        lat = math.atan2(x[2], p * (1 - dm.WGS84_E2 * n / (n + h)))
    return math.degrees(lat), lon, h


def make_xml(path, name="SAT_SM_SLC_1SVV_20300105t042344_000001_000001_000001", acq_id=1, mode="SM",
             pol="VV", orbit_dir="DESCENDING", look_side="left", look_angle=33.8, heading=192.8,
             n_az=40, n_rg=32, ati=2.358e-4, dr=0.7495, fs=200e6, bw=160e6, r0=640000.0,
             tshift=0.0, pos_offset=(0.0, 0.0, 0.0), zd_delay=0.0, half_width_deg=0.2, center_shift_deg=0.0,
             n_vec=161, vec_step=0.5, prf=None, orbit_fn=None, tag_quant=0.0, vel_noise=0.0,
             orbit_level="Rapid orbit", vel_bias=(0.0, 0.0, 0.0)):
    """写一份结构与真实产品一致的扁平 XML。zd_delay: 起始时间相对 T_REF 的延后秒数。
    orbit_fn(t) -> (位置, 速度) 可替换轨道模型; tag_quant 把 XML 里的矢量时间标签量化到该间隔 (模拟实测的标签抖动);
    vel_noise 给速度加高斯噪声 (模拟位置与速度不自洽); vel_bias 给全部速度矢量加同一个常量偏置 (模拟位置的导数
    与速度矢量系统性不一致的轨道)。"""
    ctr = _center_target()
    lat, lon, h = _ecef_to_llh(ctr)
    lon += center_shift_deg
    c = {"center": (n_az // 2, n_rg // 2, lat, lon),
         "first_near": (1, 1, lat + half_width_deg, lon - half_width_deg),
         "first_far": (n_rg, 1, lat + half_width_deg, lon + half_width_deg),
         "last_far": (n_rg, n_az, lat - half_width_deg, lon + half_width_deg),
         "last_near": (1, n_az, lat - half_width_deg, lon - half_width_deg)}
    t0 = T_REF + zd_delay - 0.5 * n_az * ati
    t1 = t0 + n_az * ati
    rows = ["<?xml version=\"1.0\" encoding=\"UTF-8\"?>", "<Metadata>", "  <Orbit_State_Vectors>",
            "    <count>%d</count>" % n_vec]
    for i in range(n_vec):
        t = T_REF + (i - n_vec // 2) * vec_step
        p, v = (orbit_fn or (lambda tt: _orbit_state(tt, tshift=tshift)))(t)
        p = p + np.array(pos_offset)
        if vel_noise:
            v = v + vel_noise * np.random.default_rng(i).standard_normal(3)
        v = v + np.array(vel_bias)
        tag = round((t - T_REF) / tag_quant) * tag_quant + T_REF if tag_quant else t
        rows += ["    <orbit_vector>", "      <number>%d</number>" % (i + 1), "      <time>%s</time>" % _fmt_time(tag)]
        # 位置/速度的小数位不能太少: 插值速度是位置样条的导数, 位置舍入 δ 会带来约 δ/步长 的速度误差,
        # 再乘上约 650 km 的斜距力臂就是零多普勒时刻误差
        rows += ["      <pos%s>%.7f</pos%s>" % (k, p["XYZ".index(k)], k) for k in "XYZ"]
        rows += ["      <vel%s>%.9f</vel%s>" % (k, v["XYZ".index(k)], k) for k in "XYZ"]
        rows += ["    </orbit_vector>"]
    rows += ["  </Orbit_State_Vectors>"]
    fields = dict(product_name=name, product_file=name + ".h5", satellite_name="Dianjian-01", acquisition_mode=mode,
                  polarization=pol, look_side=look_side, satellite_look_angle=look_angle,
                  orbit_direction=orbit_dir, orbit_absolute_number=1692, orbit_relative_number=22,
                  orbit_processing_level=orbit_level, acquisition_id=acq_id, carrier_frequency=9.6e9,
                  range_sampling_rate=fs, chirp_bandwidth=bw, chirp_duration=23.6e-6,
                  processing_prf=prf or 1.0 / ati, acquisition_prf=prf or 1.0 / ati, azimuth_time_interval=ati,
                  total_processed_bandwidth_azimuth=2727.5, number_of_azimuth_samples=n_az,
                  number_of_range_samples=n_rg, slant_range_to_first_pixel=r0, slant_range_spacing=dr,
                  azimuth_ground_spacing=1.74, acquisition_start_utc=_fmt_time(t0 - 0.5),
                  zerodoppler_start_utc=_fmt_time(t0), zerodoppler_end_utc=_fmt_time(t1), heading=heading,
                  incidence_near=35.8, incidence_center=37.0, incidence_far=38.2, avg_scene_height=h,
                  calibration_factor=5.8e-5)
    for k, v in fields.items():
        rows.append("  <%s>%s</%s>" % (k, v, k))
    for k, (s, l, la, lo) in c.items():
        rows.append("  <coord_%s>%d %d %.8f %.8f</coord_%s>" % (k, s, l, la, lo, k))
    rows += ["  <Doppler_Centroid_Coefficients>", "    <dc_coefficients_list>",
             "      <reference_pixel_time>0.00427</reference_pixel_time>"]
    for i, v in enumerate((-120.0, 3.5e5, 0.0)):
        rows += ["      <coefficient>", "        <number>%d</number>" % i, "        <value>%r</value>" % v,
                 "      </coefficient>"]
    rows += ["    </dc_coefficients_list>", "  </Doppler_Centroid_Coefficients>", "</Metadata>"]
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")
    return path


class TestParseAndOrbit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.xml = make_xml(os.path.join(cls.tmp.name, "a.XML"))
        cls.p = dm.Product(cls.xml)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_fields(self):
        p = self.p
        self.assertEqual((p.n_az, p.n_rg), (40, 32))
        self.assertEqual(p.mode, "SM")
        self.assertEqual(p.orbit_dir, "DESCENDING")
        self.assertEqual(p.look_side, "left")
        self.assertAlmostEqual(p.wavelength, 0.0312, places=3)
        self.assertEqual(p.dc_poly, [-120.0, 3.5e5, 0.0])
        self.assertEqual(len(p.sv_t), 161)
        self.assertAlmostEqual(p.t_zd_end - p.t_zd_start, 40 * 2.358e-4, places=6)
        self.assertEqual(p.date, "20300105")

    def test_rejects_old_format(self):
        bad = os.path.join(self.tmp.name, "old.xml")
        with open(bad, "w") as f:
            f.write("<product><adsHeader/></product>")
        with self.assertRaises(ValueError):
            dm.Product(bad)

    def test_hermite_matches_circle(self):
        """采样点之间的 Hermite 插值与解析圆轨道的差应小于 1 mm / 1 mm/s。"""
        for dtm in (0.0, 0.123, 0.4999, 7.77, -11.3):
            t = T_REF + dtm
            pos, vel = self.p.sat_state(t)
            pe, ve = _orbit_state(t)
            self.assertLess(np.linalg.norm(pos - pe), 1e-3)
            self.assertLess(np.linalg.norm(vel - ve), 1e-3)

    def test_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            self.p.sat_state(T_REF + 1e4)

    def test_zero_doppler_time(self):
        """目标放在 T_REF + 2.3 s 时刻卫星的正侧方 (垂直于速度), 零多普勒时刻必须还原这个时刻。"""
        t0 = T_REF + 2.3
        p, v = _orbit_state(t0)
        n = -p / np.linalg.norm(p)
        side = np.cross(n, v / np.linalg.norm(v))
        target = p + 650e3 * (0.8 * n + 0.6 * side)
        self.assertAlmostEqual(self.p.zero_doppler_time(target), t0, places=6)

    def test_target_geometry_line_col(self):
        t, _, r, line, col = self.p.target_geometry(_center_target())
        self.assertAlmostEqual(t, T_REF, places=5)
        self.assertAlmostEqual(r, 700e3, delta=1.0)
        self.assertAlmostEqual(line, 0.5 * 40, delta=0.05)
        self.assertAlmostEqual(col, (700e3 - 640000.0) / 0.7495, delta=2.0)


_DYN_P0, _DYN_V0 = _orbit_state(T_REF)


def _dyn_state(t):
    """与 dm.propagate_state 同一套动力学的"真实"轨道 (用比外推更细的步长积分)"""
    return dm.propagate_state(_DYN_P0, _DYN_V0, t - T_REF, step=0.01)


class TestOrbitTimeAxis(unittest.TestCase):
    """实时轨道的时间标签被量化 (抖动 0.1~0.4 ms), 但位置/速度本身等间隔、自洽: 应还原等间隔时间轴。"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _prod(self, name, **kw):
        kw.setdefault("orbit_fn", _dyn_state)
        return dm.Product(make_xml(os.path.join(self.tmp.name, name + ".XML"), **kw))

    def test_quantized_tags_are_regularized(self):
        step = 0.025
        p = self._prod("q", n_vec=161, vec_step=step, tag_quant=1.7e-4, orbit_level="Real Time orbit")
        self.assertTrue(p.orbit_info["regularized"], p.orbit_info)
        true_t = T_REF + (np.arange(161) - 80) * step
        self.assertGreater(np.abs(p.sv_t_tags - true_t).max(), 5e-5)          # 标签确实抖了
        self.assertLess(np.abs(p.sv_t - true_t).max(), 3e-5)                    # 规整后接近真值 (起点由平均得到)
        self.assertAlmostEqual(p.orbit_info["step"], step, delta=step * 1e-5)

    def test_exact_uniform_tags_untouched(self):
        p = self._prod("u", n_vec=161, vec_step=0.5)
        self.assertFalse(p.orbit_info["regularized"])
        np.testing.assert_array_equal(p.sv_t, p.sv_t_tags)

    def test_inconsistent_velocity_keeps_tags(self):
        p = self._prod("n", n_vec=161, vec_step=0.025, tag_quant=1.7e-4, vel_noise=8.0)
        self.assertFalse(p.orbit_info["regularized"])
        np.testing.assert_array_equal(p.sv_t, p.sv_t_tags)

    def test_too_few_vectors(self):
        t = np.arange(5) * 1.0
        _, info = dm.regularize_orbit_times(t, np.zeros((5, 3)), np.ones((5, 3)))
        self.assertFalse(info["regularized"])


class TestPropagation(unittest.TestCase):
    def test_forward_backward_returns_to_start(self):
        p1, v1 = dm.propagate_state(_DYN_P0, _DYN_V0, 12.0)
        p0, v0 = dm.propagate_state(p1, v1, -12.0)
        self.assertLess(np.linalg.norm(p0 - _DYN_P0), 1e-5)
        self.assertLess(np.linalg.norm(v0 - _DYN_V0), 1e-7)

    def test_step_convergence(self):
        a, _ = dm.propagate_state(_DYN_P0, _DYN_V0, 10.0, step=0.5)
        b, _ = dm.propagate_state(_DYN_P0, _DYN_V0, 10.0, step=0.01)
        self.assertLess(np.linalg.norm(a - b), 1e-3)

    def test_gravity_magnitude(self):
        """地面静止的状态受到的加速度应接近 g 的量级 (检验引力项符号与量纲)"""
        a = dm._accel(np.array([dm.WGS84_A + 531e3, 0.0, 0.0]), np.zeros(3))
        self.assertAlmostEqual(-a[0], 8.2, delta=0.3)
        self.assertAlmostEqual(a[1], 0.0, places=6)


class TestGmtsarOrbit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _prod(self, name, **kw):
        kw.setdefault("orbit_fn", _dyn_state)
        return dm.Product(make_xml(os.path.join(self.tmp.name, name + ".XML"), **kw))

    def test_window_only_orbit_is_extended(self):
        """轨道只覆盖成像时段 (实时轨道) -> 外推到前后各 10 s, 间隔 0.5 s, 与真值一致"""
        p = self._prod("w", n_vec=41, vec_step=0.025, tag_quant=1.7e-4, orbit_level="Real Time orbit")
        o = p.gmtsar_orbit()
        self.assertTrue(o["extended"])
        self.assertEqual(o["dt"], 0.5)
        np.testing.assert_allclose(np.diff(o["t"]), 0.5, atol=1e-9)
        self.assertLessEqual(o["t"][0], p.t_zd_start - 10.0 + 1e-9)
        self.assertGreaterEqual(o["t"][-1], p.t_zd_end + 10.0 - 1e-9)
        for k in (0, len(o["t"]) // 2, len(o["t"]) - 1):
            pt, vt = _dyn_state(o["t"][k])
            self.assertLess(np.linalg.norm(o["p"][k] - pt), 0.05)
            self.assertLess(np.linalg.norm(o["v"][k] - vt), 0.005)
        self.assertLess(o["max_dev_m"], 0.05)

    def test_xml_velocity_matches_vectors(self):
        """xml_velocity 在密集轨道与稀疏轨道 (窗口内不足 5 个矢量) 上都应还原速度矢量"""
        bias = np.array([0.0, 0.03, 0.05])
        for n, step in ((41, 0.025), (161, 0.5)):
            p = self._prod("xv%d" % n, n_vec=n, vec_step=step, vel_bias=tuple(bias), zd_delay=0.0125)
            for dt in (0.0, 0.0125, 0.3):
                t = p.t_mid + dt
                np.testing.assert_allclose(p.xml_velocity(t), _dyn_state(t)[1] + bias, atol=2e-4)

    def test_inconsistent_pv_center_velocity_follows_xml(self):
        """位置导数与速度矢量系统性不一致的轨道: LED 中心速度应取速度矢量, 不是位置样条的导数。
        场景中心落在两个矢量的中点 (样条导数在这里偏离速度矢量最多, 偏差 1.5 x 位置/速度差)。"""
        bias = np.array([0.0, 0.03, 0.05])
        p = self._prod("pv", n_vec=41, vec_step=0.025, zd_delay=0.0125, tag_quant=1.7e-4, vel_bias=tuple(bias),
                       orbit_level="Real Time orbit")
        self.assertTrue(p.orbit_info["regularized"], p.orbit_info)
        o = p.gmtsar_orbit()
        self.assertTrue(o["extended"])
        t_c = p.t_mid
        v_xml = _dyn_state(t_c)[1] + bias
        v_spline = p.sat_state(t_c)[1]
        # 样条导数确实偏了 (偏置的沿轨分量会被时间轴规整里的间隔估计吸收一部分, 所以小于 1.5 x |bias| = 0.087)
        self.assertGreater(np.linalg.norm(v_spline - v_xml), 0.04)
        # 窗口平均速度失配 = 偏置中垂直于速度的部分 (沿轨分量被时间轴规整的间隔估计吸收了)
        v_hat = _dyn_state(t_c)[1] / np.linalg.norm(_dyn_state(t_c)[1])
        b_perp = bias - (bias @ v_hat) * v_hat
        self.assertAlmostEqual(o["pv_mismatch_mps"], np.linalg.norm(b_perp), delta=1e-3)
        # LED 里离场景中心最近的点 = 以 (中心位置, XML 速度) 为起点外推的结果
        k = int(np.argmin(np.abs(o["t"] - t_c)))
        p_c = p.sat_state(t_c)[0]
        p_x, v_x = dm.propagate_state(p_c, v_xml, o["t"][k] - t_c)
        np.testing.assert_allclose(o["v"][k], v_x, atol=1e-5)
        np.testing.assert_allclose(o["p"][k], p_x, atol=1e-6)
        _, v_s = dm.propagate_state(p_c, v_spline, o["t"][k] - t_c)
        self.assertGreater(np.linalg.norm(o["v"][k] - v_s), 0.04)           # 与旧做法 (样条导数) 明显不同

    def test_consistent_pv_has_no_mismatch(self):
        """位置与速度自洽的轨道: 失配量接近 0, LED 与以样条速度为起点的外推一致"""
        p = self._prod("ok", n_vec=41, vec_step=0.025, zd_delay=0.0125, orbit_level="Real Time orbit")
        o = p.gmtsar_orbit()
        self.assertLess(o["pv_mismatch_mps"], 1e-3)
        t_c = p.t_mid
        pc, vs = p.sat_state(t_c)
        k = int(np.argmin(np.abs(o["t"] - t_c)))
        _, v_old = dm.propagate_state(pc, vs, o["t"][k] - t_c)
        np.testing.assert_allclose(o["v"][k], v_old, atol=2e-3)

    def test_max_dev_window_only_counts_imaging_vectors(self):
        p = self._prod("dw", n_vec=161, vec_step=0.025, zd_delay=0.0125, vel_bias=(0.0, 0.03, 0.05),
                       n_az=4000, ati=1.0 / 4000.0, orbit_level="Real Time orbit")
        o = p.gmtsar_orbit()
        self.assertTrue(o["extended"])
        self.assertLessEqual(o["max_dev_window_m"], o["max_dev_m"])
        self.assertGreater(o["max_dev_m"], 0.0)

    def test_long_uniform_orbit_kept_as_is(self):
        p = self._prod("l", n_vec=161, vec_step=0.5)
        o = p.gmtsar_orbit()
        self.assertFalse(o["extended"])
        self.assertEqual(len(o["t"]), 161)
        np.testing.assert_array_equal(o["p"], p.sv_p)
        self.assertLess(o["pv_mismatch_mps"], 1e-3)
        self.assertEqual(o["max_dev_window_m"], 0.0)

    def test_gmtsar_calc_dop_orb_coverage_requirement(self):
        """calc_dop_orb 在影像起止时刻前后各 ±2 s 取轨道: 生成的轨道必须全部覆盖"""
        p = self._prod("c", n_vec=41, vec_step=0.025, orbit_level="Real Time orbit")
        o = p.gmtsar_orbit()
        self.assertLess(o["t"][0], p.t_zd_start - 2.5)
        self.assertGreater(o["t"][-1], p.t_zd_end + 2.5)


class TestPairGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _prod(self, name, **kw):
        return dm.Product(make_xml(os.path.join(self.tmp.name, name + ".XML"), **kw))

    def test_identical_pair_is_zero_baseline(self):
        a, b = self._prod("m0"), self._prod("s0")
        g = dm.pair_geometry(a, b)
        self.assertAlmostEqual(g["b_perp"], 0.0, places=3)
        self.assertAlmostEqual(g["range_offset_px"], 0.0, places=3)
        self.assertAlmostEqual(g["azimuth_offset_px"], 0.0, places=3)

    def test_perpendicular_baseline_matches_analytic(self):
        """从影像轨道整体平移 b (垂直于速度): B⊥、B∥ 有解析真值。"""
        m = self._prod("m1")
        target = _center_target()
        p, v = _orbit_state(T_REF)
        u = (p - target) / np.linalg.norm(p - target)
        cross = np.cross(u, v / np.linalg.norm(v))
        along_los = u - np.dot(u, v / np.linalg.norm(v)) * v / np.linalg.norm(v)
        b = 80.0 * cross + 35.0 * along_los / np.linalg.norm(along_los)
        s = self._prod("s1", pos_offset=tuple(b))
        g = dm.pair_geometry(m, s)
        b_par = float(np.dot(b, u))
        b_perp = math.sqrt(float(np.dot(b, b)) - b_par ** 2)
        self.assertAlmostEqual(g["b_par"], b_par, delta=0.01)
        self.assertAlmostEqual(abs(g["b_perp"]), b_perp, delta=0.01)
        # 约定: 从影像坐标 = 主影像坐标 + 偏移。从影像卫星比主影像更远 (b_par > 0) 时斜距更长, 列号更大
        self.assertAlmostEqual(g["range_offset_px"], b_par / 0.7495, delta=0.2)

    def test_along_track_shift_does_not_change_baseline(self):
        """从影像整条轨道沿轨平移 30 m: 卫星只是提前/推迟经过目标 (零多普勒时刻随之移动, 由方位向配准吸收),
        B⊥/B∥ 不变; 报告的沿轨分量只剩速度方向随轨道曲率转动的几米。

        方位向偏移的真值用解析圆轨道上的二分求根算 (与被测代码的 Hermite 插值互相独立)。不能用 -30/|v|/ATI 估算:
        平移后零多普勒平面 (垂直于 v(t)) 还会随速度方向转动, 实际时间移动比 d/|v| 大约 9%, 算出来是 -18.2 px 而非 -16.7 px"""
        m = self._prod("m6")
        p, v = _orbit_state(T_REF)
        vhat = v / np.linalg.norm(v)
        u = (p - _center_target()) / np.linalg.norm(p - _center_target())
        cross = np.cross(u, vhat)
        g0 = dm.pair_geometry(m, self._prod("s6a", pos_offset=tuple(60.0 * cross)))
        shift = 60.0 * cross + 30.0 * vhat
        g1 = dm.pair_geometry(m, self._prod("s6b", pos_offset=tuple(shift)))
        self.assertAlmostEqual(g1["b_perp"], g0["b_perp"], delta=0.05)
        self.assertAlmostEqual(g1["b_par"], g0["b_par"], delta=0.05)
        self.assertAlmostEqual(g0["b_along"], 0.0, delta=0.05)
        self.assertLess(abs(g1["b_along"]), 5.0)
        target = dm.llh_to_ecef(*m.center_llh)
        expected = (_zero_doppler_analytic(target, shift) - _zero_doppler_analytic(target)) / 2.358e-4
        self.assertLess(expected, -16.0)                      # 先确认方向与量级: 卫星靠前, 更早经过目标
        self.assertAlmostEqual(g1["azimuth_offset_px"], expected, delta=0.02)
        self.assertAlmostEqual(g0["azimuth_offset_px"], (_zero_doppler_analytic(target, 60.0 * cross)
                                                         - _zero_doppler_analytic(target)) / 2.358e-4, delta=0.02)

    def test_azimuth_offset_from_start_time(self):
        m = self._prod("m2")
        s = self._prod("s2", zd_delay=0.05)
        g = dm.pair_geometry(m, s)
        self.assertAlmostEqual(g["azimuth_offset_px"], -0.05 / 2.358e-4, delta=0.05)
        self.assertAlmostEqual(g["b_perp"], 0.0, places=2)

    def test_sign_of_perp_baseline_flips(self):
        m = self._prod("m3")
        p, v = _orbit_state(T_REF)
        target = _center_target()
        u = (p - target) / np.linalg.norm(p - target)
        cross = np.cross(u, v / np.linalg.norm(v))
        gp = dm.pair_geometry(m, self._prod("sp", pos_offset=tuple(50 * cross)))
        gn = dm.pair_geometry(m, self._prod("sn", pos_offset=tuple(-50 * cross)))
        self.assertGreater(gp["b_perp"] * gn["b_perp"] * -1, 0)
        self.assertAlmostEqual(abs(gp["b_perp"]), 50.0, delta=0.01)

    def test_footprint_overlap(self):
        a = self._prod("fa")
        same = dm.footprint_overlap(a, self._prod("fb"))
        self.assertAlmostEqual(same[1], 1.0, places=3)
        half = dm.footprint_overlap(a, self._prod("fc", center_shift_deg=0.2))
        self.assertAlmostEqual(half[1], 0.5, delta=0.02)
        far = dm.footprint_overlap(a, self._prod("fd", center_shift_deg=2.0))
        self.assertEqual(far[0], 0.0)


class TestPairSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _prod(self, name, **kw):
        return dm.Product(make_xml(os.path.join(self.tmp.name, name + ".XML"), **kw))

    def test_good_pair_grade_a(self):
        m = self._prod("m", acq_id=1)
        s = self._prod("s", acq_id=2)
        s.t_zd_start += 8 * 86400.0       # 只改标注的成像时间 (时间基线 8 天), 轨道仍用同一条
        s.t_zd_end += 8 * 86400.0
        r = sel.evaluate_pair(m, s)
        self.assertEqual(r["grade"], "A", r["hard"] + r["soft"])

    def test_mode_mismatch_rejected(self):
        r = sel.evaluate_pair(self._prod("m2", acq_id=1), self._prod("s2", acq_id=2, mode="SP"))
        self.assertEqual(r["grade"], "X")
        self.assertTrue(any("模式" in h for h in r["hard"]))

    def test_orbit_direction_mismatch_rejected(self):
        r = sel.evaluate_pair(self._prod("m3", acq_id=1), self._prod("s3", acq_id=2, orbit_dir="ASCENDING"))
        self.assertEqual(r["grade"], "X")

    def test_large_baseline_rejected(self):
        m = self._prod("m4", acq_id=1)
        p, v = _orbit_state(T_REF)
        u = (p - _center_target()) / np.linalg.norm(p - _center_target())
        cross = np.cross(u, v / np.linalg.norm(v))
        s = self._prod("s4", acq_id=2, pos_offset=tuple(5000.0 * cross))
        r = sel.evaluate_pair(m, s)
        self.assertEqual(r["grade"], "X")
        self.assertTrue(any("垂直基线" in h for h in r["hard"]))

    def test_low_overlap_rejected(self):
        r = sel.evaluate_pair(self._prod("m5", acq_id=1), self._prod("s5", acq_id=2, center_shift_deg=0.33))
        self.assertEqual(r["grade"], "X")
        self.assertTrue(any("重叠" in h for h in r["hard"]))

    def test_select_scans_directory(self):
        d = os.path.join(self.tmp.name, "scan")
        os.makedirs(d)
        make_xml(os.path.join(d, "ICEYE_SM_SLC_a.XML"), name="ICEYE_SM_SLC_a", acq_id=11)
        make_xml(os.path.join(d, "ICEYE_SM_SLC_b.XML"), name="ICEYE_SM_SLC_b", acq_id=12)
        prods, pairs = sel.select(d, include_all=True)
        self.assertEqual(len(prods), 2)
        self.assertEqual(len(pairs), 1)


@unittest.skipUnless(HAVE_H5PY, "需要 h5py")
class TestConverter(unittest.TestCase):
    N_AZ, N_RG = 42, 30        # 都不是 4 的倍数: 转换后应裁成 40 x 28

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.xml = make_xml(os.path.join(cls.tmp.name, "p.XML"), n_az=cls.N_AZ, n_rg=cls.N_RG)
        rng = np.random.default_rng(3)
        cls.z = (rng.standard_normal((cls.N_AZ, cls.N_RG)) + 1j * rng.standard_normal((cls.N_AZ, cls.N_RG))) * 40
        cls.z = cls.z.astype(np.complex64)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _write_h5(self, name, layout):
        import h5py
        path = os.path.join(self.tmp.name, name + ".h5")
        with h5py.File(path, "w") as f:
            if layout == "split":
                f.create_dataset("s_i", data=self.z.real.copy())
                f.create_dataset("s_q", data=self.z.imag.copy())
            elif layout == "split_transposed":
                f.create_dataset("s_i", data=self.z.real.T.copy())
                f.create_dataset("s_q", data=self.z.imag.T.copy())
            elif layout == "complex":
                f.create_dataset("slc", data=self.z)
        return path

    def _convert(self, layout, scale=2.0):
        h5 = self._write_h5("d_" + layout, layout)
        out = os.path.join(self.tmp.name, "out_" + layout, "20300105")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        conv.convert(self.xml, h5, out, scale=scale, block_lines=16, quiet=True)
        return out

    def _check_slc(self, out, scale=2.0):
        raw = np.fromfile(out + ".SLC", dtype=np.int16)
        self.assertEqual(raw.size, 40 * 28 * 2)
        raw = raw.reshape(40, 56)
        got = raw[:, 0::2] + 1j * raw[:, 1::2]
        want = np.rint(self.z[:40, :28].real * scale) + 1j * np.rint(self.z[:40, :28].imag * scale)
        np.testing.assert_array_equal(got.real, want.real)
        np.testing.assert_array_equal(got.imag, want.imag)

    def test_split_layout(self):
        out = self._convert("split")
        self._check_slc(out)

    def test_transposed_layout_is_normalized(self):
        out = self._convert("split_transposed")
        self._check_slc(out)

    def test_complex_layout(self):
        out = self._convert("complex")
        self._check_slc(out)

    def test_prm_contents(self):
        out = self._convert("split")
        prm = {}
        with open(out + ".PRM") as f:
            for line in f:
                if "=" in line:
                    k, v = line.split("=", 1)
                    prm[k.strip()] = v.strip()
        self.assertEqual(int(prm["num_rng_bins"]), 28)
        self.assertEqual(int(prm["num_lines"]), 40)
        self.assertEqual(int(prm["bytes_per_line"]), 4 * 28)
        self.assertAlmostEqual(float(prm["PRF"]), 1.0 / 2.358e-4, places=2)
        self.assertAlmostEqual(float(prm["near_range"]), 640000.0, places=2)
        self.assertEqual(prm["lookdir"], "L")
        self.assertEqual(prm["orbdir"], "D")
        # 文件名里不能带目录 (GMTSAR 在 PRM 所在目录里按文件名读)
        self.assertEqual(prm["SLC_file"], "20300105.SLC")
        self.assertEqual(prm["led_file"], "20300105.LED")
        self.assertAlmostEqual(float(prm["radar_wavelength"]), dm.C_LIGHT / 9.6e9, places=7)
        # clock_start 为年内天数 + 当天秒/86400, 与 SC_clock_start 的小数部分一致
        sc = float(prm["SC_clock_start"])
        self.assertEqual(int(sc // 1000), 2030)
        self.assertAlmostEqual(sc - 2030000, float(prm["clock_start"]), places=8)
        # 零多普勒产品: PRM 的 fd1 必须是 0, XML 的多普勒常数项 (-120 Hz) 只记在转换信息里
        self.assertEqual(float(prm["fd1"]), 0.0)
        with open(out + ".convert.txt") as f:
            info = dict(ln.strip().split(" = ") for ln in f if " = " in ln)
        self.assertAlmostEqual(float(info["xml_doppler_fd1"]), -120.0, places=6)

    @unittest.skipUnless(HAVE_GMTSAR, "需要 GMTSAR 的 calc_dop_orb 和 SAT_llt2rat")
    def test_gmtsar_backgeocoding_matches_zero_doppler_geometry(self):
        """独立核对: 用 GMTSAR 自己的 SAT_llt2rat 把地面点投影回转换后 SLC 的 (距离, 方位) 像素, 应与 dj1_meta 的零多普勒
        几何一致。这检验了 PRM/LED 的时间、近距、PRF 约定能被 GMTSAR 正确解释。
        同时验证这个检验对 fd1 敏感: 若按 XML 的多普勒常数项 (-120 Hz) 写 fd1, 方位坐标会错开上百行。"""
        # 合成场景的目标点在卫星"向下 x 速度"一侧, 即右侧; SAT_llt2rat 的精确模式会按 PRM 的 lookdir 检查侧视方向, 必须一致
        xml = make_xml(os.path.join(self.tmp.name, "gm.XML"), n_az=self.N_AZ, n_rg=self.N_RG, look_side="right")
        out = os.path.join(self.tmp.name, "gm_out", "20300105")
        os.makedirs(os.path.dirname(out))
        conv.convert(xml, self._write_h5("d_gm", "split"), out, scale=2.0, block_lines=16, quiet=True)
        wd = os.path.dirname(out)
        prod = dm.Product(xml)
        ctr = _center_target()
        t0, p0, r0, _, _ = prod.target_geometry(ctr)
        vhat = prod.sat_state(t0)[1]
        speed = np.linalg.norm(vhat)
        vhat = vhat / speed
        pts = [ctr + vhat * dl * prod.ati * speed - (p0 - ctr) / r0 * dr * prod.dr
               for dl in (-15, 0, 15) for dr in (-20, 0, 20)]
        stdin = "".join("%.9f %.9f %.4f\n" % (lon, lat, h) for lat, lon, h in map(_ecef_to_llh, pts))

        def project(fd1_arg):
            """按 GMTSAR 惯例 (PRM0 + calc_dop_orb 的输出) 生成 PRM, 再用 SAT_llt2rat 精确模式反投影"""
            shutil.copy(out + ".PRM", os.path.join(wd, "g.PRM0"))
            subprocess.run(["calc_dop_orb", "g.PRM0", "g.log", "0", fd1_arg], cwd=wd, check=True, capture_output=True)
            with open(os.path.join(wd, "g.PRM"), "w") as f:
                for part in ("g.PRM0", "g.log"):
                    with open(os.path.join(wd, part)) as g:
                        f.write(g.read())
            res = subprocess.run(["SAT_llt2rat", "g.PRM", "1"], cwd=wd, input=stdin, check=True, capture_output=True,
                                 text=True).stdout
            return [list(map(float, ln.split())) for ln in res.splitlines() if ln.strip()]

        with open(out + ".PRM") as f:
            fd1_in_file = float(re.search(r"^fd1\s*=\s*(\S+)", f.read(), re.M).group(1))
        got = project("%.6f" % fd1_in_file)
        wrong = project("-120")
        self.assertEqual(len(got), len(pts))
        pixel = dm.C_LIGHT / (2.0 * prod.fs)    # GMTSAR 由 rng_samp_rate 推距离像元, 不读 XML 的 slant_range_spacing
        for p, g, w in zip(pts, got, wrong):
            _, _, r, line, _ = prod.target_geometry(p)
            # SAT_llt2rat 精确模式的输出: 距离、方位都从 0 起算, 与 dj1_meta.target_geometry 同一约定
            self.assertAlmostEqual(g[0], (r - prod.r0) / pixel, delta=0.1)
            self.assertAlmostEqual(g[1], line, delta=0.05)
            self.assertGreater(abs(w[1] - g[1]), 50.0)

    def test_led_contents(self):
        out = self._convert("split")
        with open(out + ".LED") as f:
            rows = [ln.split() for ln in f if ln.strip()]
        n = int(rows[0][0])
        self.assertEqual(n, 161)
        self.assertEqual(len(rows) - 1, n)
        self.assertAlmostEqual(float(rows[0][4]), 0.5, places=6)
        p = dm.Product(self.xml)
        self.assertAlmostEqual(float(rows[1][3]), p.sv_p[0][0], places=2)
        self.assertAlmostEqual(float(rows[-1][8]), p.sv_v[-1][2], places=5)
        # 起始时刻 (年, 年内日, 当天秒) 能还原为 XML 里的 UTC
        y, j, s = int(rows[1][0]), int(rows[1][1]), float(rows[1][2])
        t0 = dt.datetime(y, 1, 1) + dt.timedelta(days=j - 1, seconds=s)
        self.assertAlmostEqual((t0 - dm.from_sec(p.sv_t[0])).total_seconds(), 0.0, places=5)

    def _gmtsar_time_axis(self, out):
        """按 GMTSAR 的规则解释 LED 与 PRM: 返回 (轨道首点时间, dsec, 点数, 场景起止时间), 单位秒 (从 1 月 1 日起算, 日序不减 1)"""
        with open(out + ".LED") as f:
            head = f.readline().split()
            nrows = sum(1 for ln in f if ln.strip())
        nd, iy, idy, sec, dsec = int(head[0]), int(head[1]), int(head[2]), float(head[3]), float(head[4])
        self.assertEqual(nrows, nd)
        prm = {}
        with open(out + ".PRM") as f:
            for ln in f:
                if "=" in ln:
                    k, v = ln.split("=", 1)
                    prm[k.strip()] = v.strip()
        t_start = 86400.0 * float(prm["clock_start"])
        t_end = t_start + int(prm["num_lines"]) / float(prm["PRF"])
        return 86400.0 * idy + sec, dsec, nd, t_start, t_end

    def test_led_matches_prm_time_axis(self):
        out = self._convert("split")
        pt0, dsec, nd, t_start, t_end = self._gmtsar_time_axis(out)
        self.assertGreater(t_start, pt0)
        self.assertLess(t_end, pt0 + (nd - 1) * dsec)

    def test_window_only_orbit_covers_calc_dop_orb_margin(self):
        """实时轨道只覆盖成像时段, 转换后的 LED 必须延伸到场景前后至少 2.5 s (GMTSAR calc_dop_orb 的要求)"""
        xml = make_xml(os.path.join(self.tmp.name, "w.XML"), n_az=self.N_AZ, n_rg=self.N_RG, n_vec=41, vec_step=0.025,
                       orbit_fn=_dyn_state, tag_quant=1.7e-4, orbit_level="Real Time orbit")
        h5 = self._write_h5("d_w", "split")
        out = os.path.join(self.tmp.name, "w_out", "20300105")
        os.makedirs(os.path.dirname(out))
        conv.convert(xml, h5, out, scale=2.0, block_lines=16, quiet=True)
        pt0, dsec, nd, t_start, t_end = self._gmtsar_time_axis(out)
        self.assertAlmostEqual(dsec, 0.5, places=9)
        self.assertLess(pt0, t_start - 2.5)
        self.assertGreater(pt0 + (nd - 1) * dsec, t_end + 2.5)
        with open(out + ".convert.txt") as f:
            info = dict(ln.strip().split(" = ") for ln in f if " = " in ln)
        self.assertEqual(info["led_extended"], "True")
        self.assertEqual(info["orbit_time_regularized"], "True")

    def test_metadata_only_rewrites_prm_and_led(self):
        out = self._convert("split")
        for ext in (".PRM", ".LED"):
            os.remove(out + ext)
        conv.metadata_only(self.xml, out, quiet=True)
        self.assertTrue(os.path.isfile(out + ".PRM") and os.path.isfile(out + ".LED"))
        with open(out + ".convert.txt") as f:
            info = dict(ln.strip().split(" = ") for ln in f if " = " in ln)
        self.assertIn("scale", info)                     # 沿用之前转换时的比例因子
        os.truncate(out + ".SLC", 100)
        with self.assertRaises(ValueError):
            conv.metadata_only(self.xml, out, quiet=True)

    def test_auto_scale_and_clipping(self):
        h5 = self._write_h5("d_auto", "split")
        out = os.path.join(self.tmp.name, "auto")
        conv.convert(self.xml, h5, out, block_lines=16, quiet=True)
        with open(out + ".convert.txt") as f:
            txt = dict(ln.strip().split(" = ") for ln in f if " = " in ln)
        self.assertGreater(float(txt["scale"]), 0)
        self.assertEqual(int(txt["clipped_samples"]), 0)
        raw = np.fromfile(out + ".SLC", dtype=np.int16)
        self.assertLessEqual(np.abs(raw).max(), 32767)


if __name__ == "__main__":
    unittest.main()
