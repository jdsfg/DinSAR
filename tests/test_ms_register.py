#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ms_register 合成数据测试 (不需要 GMTSAR 二进制)

构造带已知仿射偏移与指定相干性的 SLC 对，检查:
  1. 相干性足够时, 恢复的偏移模型与真值在影像四角处相差 < 0.1 像素, 且质量门限通过;
  2. PRM 中给出严重错误的初始 rshift/ashift 不影响结果 (多尺度估计不依赖初始偏移);
  3. 相干性很低时, 要么门限不通过 (退出码 3), 要么结果仍然准确 (绝不允许"通过但不准");
  4. 两幅不相关的影像必须被拒绝 (退出码 3);
  5. 写回 PRM 的参数与 GMTSAR resamp 的坐标约定一致;
  6. verify 能区分已对齐与残余 0.6 像素偏差的影像。

运行: python3 -m unittest discover -s tests -v    (需要 numpy 与 scipy, 约 3-4 分钟)
"""

import json
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

try:
    from scipy import ndimage
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gmtsar-master", "gmtsar", "csh"))
import ms_register as ms  # noqa: E402

NX, NY = 1700, 2000
PRM_TEMPLATE = """num_valid_az = {ny}
nrows = {ny}
first_line = 1
st_rng_bin = 1
num_rng_bins = {nx}
bytes_per_line = {bpl}
good_bytes_per_line = {bpl}
num_lines = {ny}
num_patches = 1
PRF = 4105.09
near_range = 626552.19
rng_samp_rate = 120000000.0
radar_wavelength = 0.0555171
rshift = {rshift}
ashift = {ashift}
sub_int_r = 0
sub_int_a = 0
stretch_r = 0
stretch_a = 0
a_stretch_r = 0
a_stretch_a = 0
SLC_file = {slc}
"""

# 真值: 从影像坐标 = 主影像坐标 + 偏移
TRUTH = dict(cr=(12.37, 3.0e-4, -1.0e-4), ca=(-5.62, 2.0e-4, 1.5e-4))
ZERO = dict(cr=(0.0, 0.0, 0.0), ca=(0.0, 0.0, 0.0))


def _bandlimit_filter(shape, ovs, f_pass=0.34, f_stop=0.42):
    """可分离升余弦低通；ovs 为该网格相对 1 倍采样的过采样倍数，截止频率以 1 倍网格的周/像素计。
    真实 SLC 在距离/方位向都有约 1.2 倍以上的过采样，合成数据必须同样不混叠，
    否则幅度互相关的亚像素偏差会被人为放大。"""
    def axis(n):
        f = np.abs(np.fft.fftfreq(n)) * ovs
        t = np.clip((f - f_pass) / (f_stop - f_pass), 0.0, 1.0)
        return 0.5 * (1 + np.cos(np.pi * t))
    return (axis(shape[0])[:, None] * axis(shape[1])[None, :]).astype(np.float32)


def _complex_noise(rng, shape):
    return rng.standard_normal(shape, dtype=np.float32) + 1j * rng.standard_normal(shape, dtype=np.float32)


def _smooth_field(rng, shape, sigma):
    f = ndimage.gaussian_filter(rng.standard_normal(shape, dtype=np.float32), sigma, mode="wrap")
    return f / f.std()


_SCENES = {}


def make_scene(seed):
    """主影像坐标下 2 倍过采样的连续场: 复数散斑与反射率 (带限, 与真实 SLC 的采样条件相当)。"""
    if seed in _SCENES:
        return _SCENES[seed]
    rng = np.random.default_rng(seed)
    shape = (2 * NY, 2 * NX)
    speckle = np.fft.ifft2(np.fft.fft2(_complex_noise(rng, shape)) * _bandlimit_filter(shape, 2)).astype(np.complex64)
    speckle /= np.sqrt(np.mean(np.abs(speckle) ** 2))
    tex = (0.9 * _smooth_field(rng, shape, 3) + 0.9 * _smooth_field(rng, shape, 12)
           + 1.2 * _smooth_field(rng, shape, 60))
    refl = np.exp(0.7 * tex).astype(np.float32)
    _SCENES.clear()
    _SCENES[seed] = (speckle, refl)
    return _SCENES[seed]


def _sample(a, yy, xx):
    kw = dict(order=5, mode="constant", cval=0.0)
    if np.iscomplexobj(a):
        return (ndimage.map_coordinates(a.real, [yy, xx], **kw) + 1j * ndimage.map_coordinates(a.imag, [yy, xx], **kw))
    return ndimage.map_coordinates(a, [yy, xx], **kw)


def _write_slc(path, z):
    peak = np.percentile(np.abs(z), 99.5)
    z = z * (6000.0 / peak)
    out = np.empty((z.shape[0], 2 * z.shape[1]), np.int16)
    out[:, 0::2] = np.clip(np.rint(z.real), -32000, 32000)
    out[:, 1::2] = np.clip(np.rint(z.imag), -32000, 32000)
    out.tofile(path)


def write_prm(path, slc, rshift=0, ashift=0):
    with open(path, "w") as f:
        f.write(PRM_TEMPLATE.format(nx=NX, ny=NY, bpl=4 * NX, slc=slc, rshift=rshift, ashift=ashift))


def make_pair(tmpdir, rho, truth=TRUTH, seed=1, unrelated=False):
    """生成 m.SLC / s.SLC 及 PRM。rho 为主从散斑相干性 (复相关系数模)。"""
    speckle, refl = make_scene(seed)
    yy, xx = np.mgrid[0:NY, 0:NX].astype(np.float64)
    master = _sample(speckle, 2 * yy, 2 * xx) * _sample(refl, 2 * yy, 2 * xx)

    cr, ca = truth["cr"], truth["ca"]
    # 从影像像素 (x', y') 对应的主影像坐标: 解 [x'; y'] = [x; y] + [cr0 + cr1 x + cr2 y; ca0 + ca1 x + ca2 y]
    inv = np.linalg.inv(np.array([[1 + cr[1], cr[2]], [ca[1], 1 + ca[2]]]))
    xm = inv[0, 0] * (xx - cr[0]) + inv[0, 1] * (yy - ca[0])
    ym = inv[1, 0] * (xx - cr[0]) + inv[1, 1] * (yy - ca[0])
    del xx, yy
    rng = np.random.default_rng(seed + 100)
    indep = np.fft.ifft2(np.fft.fft2(_complex_noise(rng, (NY, NX))) * _bandlimit_filter((NY, NX), 1))
    indep /= np.sqrt(np.mean(np.abs(indep) ** 2))
    if unrelated:
        other, other_refl = make_scene(seed + 777)
        yy2, xx2 = np.mgrid[0:NY, 0:NX].astype(np.float64)
        sl = _sample(other, 2 * yy2, 2 * xx2) * _sample(other_refl, 2 * yy2, 2 * xx2)
    else:
        sl = (rho * _sample(speckle, 2 * ym, 2 * xm) + np.sqrt(max(1 - rho * rho, 0.0)) * indep) \
            * _sample(refl, 2 * ym, 2 * xm)
    _write_slc(os.path.join(tmpdir, "m.SLC"), master)
    _write_slc(os.path.join(tmpdir, "s.SLC"), sl)
    write_prm(os.path.join(tmpdir, "m.PRM"), "m.SLC")
    write_prm(os.path.join(tmpdir, "s.PRM"), "s.SLC")
    return os.path.join(tmpdir, "m.PRM"), os.path.join(tmpdir, "s.PRM")


def model_error(model, truth, nx=NX, ny=NY):
    """主影像四角及中心处模型与真值的最大偏差 (像素)。"""
    t = ms.Affine(truth["cr"], truth["ca"])
    worst = 0.0
    for (x, y) in [(0, 0), (nx, 0), (0, ny), (nx, ny), (nx / 2, ny / 2)]:
        a, b = model.predict(x, y)
        c, d = t.predict(x, y)
        worst = max(worst, abs(a - c), abs(b - d))
    return worst


def model_from_prm(path):
    prm = ms.read_prm(path)
    return ms.Affine(
        [float(prm["rshift"]) + float(prm["sub_int_r"]), float(prm["stretch_r"]), float(prm["a_stretch_r"])],
        [float(prm["ashift"]) + float(prm["sub_int_a"]), float(prm["stretch_a"]), float(prm["a_stretch_a"])])


def _estimate(d, out_name="s_out.PRM"):
    return ms.main(["estimate", os.path.join(d, "m.PRM"), os.path.join(d, "s.PRM"),
                    "-o", os.path.join(d, out_name), "--report", os.path.join(d, "rep.json"),
                    "--log", os.path.join(d, "log.txt"), "--points", os.path.join(d, "pts.dat"),
                    "--write-on-fail", "-q"])


def _load(path):
    with open(path) as f:
        return json.load(f)


@unittest.skipUnless(HAVE_SCIPY, "合成数据生成需要 scipy")
class SyntheticRegistration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ms_register_test_")
        ms._QUIET = True

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _pair(self, name, rho, **kw):
        d = os.path.join(self.tmp, name)
        os.makedirs(d, exist_ok=True)
        make_pair(d, rho, **kw)
        return d

    def test_recover_affine_and_ignore_wrong_prm_shift(self):
        d = self._pair("hi", 0.9)
        self.assertEqual(_estimate(d), ms.EXIT_OK)
        self.assertLess(model_error(model_from_prm(os.path.join(d, "s_out.PRM")), TRUTH), 0.1)
        # 初始 rshift/ashift 严重错误 (超出任何固定搜索窗) 不应影响结果
        write_prm(os.path.join(d, "s.PRM"), "s.SLC", rshift=-300, ashift=250)
        self.assertEqual(_estimate(d, "s_out2.PRM"), ms.EXIT_OK)
        self.assertLess(model_error(model_from_prm(os.path.join(d, "s_out2.PRM")), TRUTH), 0.1)

    def test_never_pass_with_wrong_answer(self):
        """相干性逐级降低: 通过则必须准确, 不准确则必须不通过。"""
        for rho in (0.5, 0.3, 0.15, 0.08):
            with self.subTest(rho=rho):
                d = self._pair("rho%02d" % int(rho * 100), rho)
                rc = _estimate(d)
                if rc == ms.EXIT_OK:
                    rep = _load(os.path.join(d, "rep.json"))
                    model = ms.Affine(rep["model"]["cr"], rep["model"]["ca"])
                    self.assertLess(model_error(model, TRUTH), 0.15, "rho=%.2f 门限通过但模型误差过大" % rho)
                else:
                    self.assertEqual(rc, ms.EXIT_GATE)

    def test_unrelated_images_rejected(self):
        d = self._pair("unrel", 0.9, unrelated=True)
        self.assertEqual(_estimate(d), ms.EXIT_GATE)

    def test_prm_updates_follow_resamp_convention(self):
        a = ms.Affine([-3.25, 1.0e-4, -2.0e-5], [7.75, 5.0e-6, -4.0e-5])
        u = a.prm_updates()
        self.assertEqual(u["rshift"], -4)
        self.assertAlmostEqual(float(u["sub_int_r"]), 0.75, places=6)
        self.assertEqual(u["ashift"], 7)
        self.assertAlmostEqual(float(u["sub_int_a"]), 0.75, places=6)
        x, y = 5000.0, 9000.0
        # resamp: ras = ram + (rshift + sub_int) + x*stretch + y*a_stretch
        dr = u["rshift"] + float(u["sub_int_r"]) + x * float(u["stretch_r"]) + y * float(u["a_stretch_r"])
        da = u["ashift"] + float(u["sub_int_a"]) + x * float(u["stretch_a"]) + y * float(u["a_stretch_a"])
        pr, pa = a.predict(x, y)
        self.assertAlmostEqual(dr, pr, places=5)
        self.assertAlmostEqual(da, pa, places=5)

    def test_verify_separates_aligned_from_misaligned(self):
        d = self._pair("aligned", 0.9, truth=ZERO)
        self.assertEqual(ms.main(["verify", os.path.join(d, "m.PRM"), os.path.join(d, "s.PRM"), "-q"]), ms.EXIT_OK)
        d2 = self._pair("shifted", 0.9, truth=dict(cr=(0.6, 0.0, 0.0), ca=(0.0, 0.0, 0.0)))
        self.assertEqual(ms.main(["verify", os.path.join(d2, "m.PRM"), os.path.join(d2, "s.PRM"), "-q"]),
                         ms.EXIT_GATE)


if __name__ == "__main__":
    unittest.main()
