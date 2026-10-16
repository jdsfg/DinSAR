#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""多尺度 SAR SLC 配准估计 (GMTSAR PRM 兼容)

不依赖轨道/PRM 中的初始偏移 (rshift/ashift)，按

    全局尺度 -> 512 -> 256 -> 128 (平原) -> 64 (山地)

逐级缩小匹配窗口、缩小搜索范围，估计从影像相对主影像的仿射偏移模型，
并把结果写回从影像 PRM (rshift/sub_int_r/stretch_r/a_stretch_r 及方位向同名参数)，
之后由 GMTSAR resamp 完成重采样。每一级都用 RANSAC + IRLS 做几何一致性检验，
质量不达标时以退出码 3 报告，而不是悄悄输出一个错误的配准结果。

偏移约定与 GMTSAR 一致:  从影像坐标 = 主影像坐标 + 偏移
    dr(x, y) = (rshift + sub_int_r) + stretch_r * x + a_stretch_r * y
    da(x, y) = (ashift + sub_int_a) + stretch_a * x + a_stretch_a * y
x 为距离向像素序号，y 为方位向行号 (均从 0 起，取主影像坐标)。

子命令:
    register MASTER.PRM SLAVE.PRM       一步完成: 估计 -> resamp 重采样 -> 独立复核 (平台/流水线入口)
    estimate MASTER.PRM SLAVE.PRM       只估计偏移模型, 写出估计后的 PRM, 输出 JSON 报告/点位/文本日志
    verify   MASTER.PRM SLAVE.PRM       对已配准 (重采样后) 的从影像做独立点位复核

退出码: 0 通过 | 2 输入错误 | 3 质量门限未通过 | 4 内部错误 (含 resamp 失败)

依赖: python3 + numpy (scipy 存在时自动用 scipy.fft 加速); register 另需 GMTSAR resamp
"""

import argparse
import copy
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

try:
    import scipy.fft as _fft
    _SCIPY_FFT = True
except Exception:  # scipy 只是加速项，不是硬依赖
    import numpy.fft as _fft
    _SCIPY_FFT = False

__version__ = "1.0"

EXIT_OK, EXIT_INPUT, EXIT_GATE, EXIT_INTERNAL = 0, 2, 3, 4

REG_KEYS = ("rshift", "sub_int_r", "stretch_r", "a_stretch_r",
            "ashift", "sub_int_a", "stretch_a", "a_stretch_a")

# 逐级匹配计划: 窗口边长 W、搜索半径 R (像素, 相对上一级模型预测)、
# 点阵 (距离向 x 方位向)、几何一致性内点阈值 thr (像素)、
# 通过一级所需的最少内点数与内点比例。
DEFAULT_SCHEDULE = (
    dict(W=512, R=64, grid=(7, 7), thr=2.5, min_inl=10, min_frac=0.20, ovs=(1, 1)),
    dict(W=256, R=16, grid=(12, 12), thr=1.2, min_inl=15, min_frac=0.25, ovs=(1, 1)),
    dict(W=128, R=8, grid=(20, 20), thr=0.7, min_inl=25, min_frac=0.25, ovs=(2, 2)),
    dict(W=64, R=6, grid=(28, 28), thr=0.5, min_inl=40, min_frac=0.25, ovs=(2, 2)),
)
CHECK_OVS = (2, 2)

GLOBAL_MIN_NCC = 0.10   # 全局多视相关峰低于此值视为不可信 (有效像对实测 0.77-0.84, 无关像对 0.03)
GLOBAL_MIN_PSR = 8.0

_WORKERS = 1
_QUIET = False


class InputError(Exception):
    pass


def log(msg):
    if not _QUIET:
        print("[ms_register] %s" % msg, file=sys.stderr, flush=True)


# ----------------------------------------------------------------------------
# FFT 封装
# ----------------------------------------------------------------------------
def _rfft2(a, s):
    if _SCIPY_FFT:
        return _fft.rfft2(a, s=s, workers=_WORKERS)
    return _fft.rfft2(a, s=s)


def _irfft2(a, s):
    if _SCIPY_FFT:
        return _fft.irfft2(a, s=s, workers=_WORKERS)
    return _fft.irfft2(a, s=s)


def _fft2(a):
    if _SCIPY_FFT:
        return _fft.fft2(a, workers=_WORKERS)
    return _fft.fft2(a)


def _ifft2(a):
    if _SCIPY_FFT:
        return _fft.ifft2(a, workers=_WORKERS)
    return _fft.ifft2(a)


def _fast_len(n):
    n = int(n)
    if _SCIPY_FFT:
        return _fft.next_fast_len(n, real=True)
    m = n
    while True:
        k = m
        for p in (2, 3, 5):
            while k % p == 0:
                k //= p
        if k == 1:
            return m
        m += 1


# ----------------------------------------------------------------------------
# PRM 与 SLC
# ----------------------------------------------------------------------------
def read_prm(path):
    prm = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                prm[k.strip()] = v.strip()
    return prm


def write_prm(src, dst, updates):
    """复制 src 到 dst，并让 updates 中的键只出现一次且位于文件末尾。"""
    with open(src, "r", errors="replace") as f:
        lines = f.read().splitlines()
    kept = [ln for ln in lines if ln.split("=", 1)[0].strip() not in updates]
    for k, v in updates.items():
        kept.append("%s = %s" % (k, v))
    tmp = dst + ".tmp"
    with open(tmp, "w") as f:
        f.write("\n".join(kept) + "\n")
    os.replace(tmp, dst)


class SLC:
    """GMTSAR SLC: 复数 int16 (I,Q 交错)，行主序，每行 num_rng_bins 个样点。"""

    def __init__(self, prm_path, slc_path=None):
        self.prm_path = os.path.abspath(prm_path)
        if not os.path.isfile(self.prm_path):
            raise InputError("PRM 不存在: %s" % prm_path)
        self.prm = read_prm(self.prm_path)
        try:
            self.nx = int(float(self.prm["num_rng_bins"]))
        except (KeyError, ValueError):
            raise InputError("PRM 缺少 num_rng_bins: %s" % prm_path)
        name = slc_path or self.prm.get("SLC_file")
        if not name:
            raise InputError("PRM 缺少 SLC_file: %s" % prm_path)
        cands = [name] if os.path.isabs(name) else [
            os.path.join(os.path.dirname(self.prm_path), name), name]
        self.path = next((c for c in cands if os.path.exists(c)), None)
        if self.path is None:
            raise InputError("找不到 SLC 文件: %s (来自 %s)" % (name, prm_path))
        rows_file = os.path.getsize(self.path) // (4 * self.nx)
        try:
            rows_prm = int(float(self.prm.get("num_patches", 1))) * int(float(self.prm["num_valid_az"]))
        except (KeyError, ValueError):
            rows_prm = rows_file
        self.ny = min(rows_file, rows_prm) if rows_prm > 0 else rows_file
        if self.ny < 64 or self.nx < 64:
            raise InputError("SLC 尺寸过小: %s (%d x %d)" % (self.path, self.nx, self.ny))
        self.mm = np.memmap(self.path, dtype=np.int16, mode="r", shape=(self.ny, 2 * self.nx))

    def read_cplx(self, y0, x0, h, w):
        """读取 (h, w) 复数块，越界部分补 0。"""
        out = np.zeros((h, w), np.complex64)
        ya, yb = max(y0, 0), min(y0 + h, self.ny)
        xa, xb = max(x0, 0), min(x0 + w, self.nx)
        if ya >= yb or xa >= xb:
            return out
        blk = np.asarray(self.mm[ya:yb, 2 * xa:2 * xb], dtype=np.float32)
        out[ya - y0:yb - y0, xa - x0:xb - x0] = blk[:, 0::2] + 1j * blk[:, 1::2]
        return out


def multilook_power(slc, looks, chunk_rows=512):
    nby, nbx = slc.ny // looks, slc.nx // looks
    out = np.empty((nby, nbx), np.float32)
    step = max(1, chunk_rows // looks)
    for b0 in range(0, nby, step):
        b1 = min(b0 + step, nby)
        blk = np.asarray(slc.mm[b0 * looks:b1 * looks, :2 * nbx * looks], dtype=np.float32)
        np.square(blk, out=blk)
        # 先沿内存连续的距离向求和，再沿方位向求和，比对四维数组做 mean(axis=(1, 3)) 快约 1.5 倍
        rng_sum = (blk[:, 0::2] + blk[:, 1::2]).reshape(-1, nbx, looks).sum(axis=2)
        out[b0:b1] = rng_sum.reshape(b1 - b0, looks, nbx).sum(axis=1) / float(looks * looks)
    return out


# ----------------------------------------------------------------------------
# 图像基本运算
# ----------------------------------------------------------------------------
def box_mean(img, k):
    k = max(1, int(k))
    if k == 1:
        return img.astype(np.float32, copy=True)
    lo = k // 2
    hi = k - 1 - lo
    p = np.pad(img.astype(np.float64), ((lo, hi), (lo, hi)), mode="edge")
    c = np.cumsum(p, axis=0)
    c = np.concatenate([np.zeros((1, c.shape[1])), c], axis=0)
    r = c[k:] - c[:-k]
    c = np.cumsum(r, axis=1)
    c = np.concatenate([np.zeros((c.shape[0], 1)), c], axis=1)
    r = c[:, k:] - c[:, :-k]
    return (r / (k * k)).astype(np.float32)


def window_sums(a, h, w):
    """所有 h x w 窗口的和，形状 (H-h+1, W-w+1)。"""
    c = np.cumsum(np.cumsum(a, axis=0, dtype=np.float64), axis=1)
    c = np.pad(c, ((1, 0), (1, 0)))
    return c[h:, w:] - c[:-h, w:] - c[h:, :-w] + c[:-h, :-w]


def prep_amp(a, mode):
    if mode == "sqrt":
        return np.sqrt(a)
    if mode == "log":
        return np.log1p(a / (float(a.mean()) + 1e-6))
    return a


def amp_upsampled(z, ovs_a, ovs_r, crop):
    """复数块经 FFT 补零过采样后取幅度，并去掉四周各 crop 个原始像素 (消除边缘振铃)。
    幅度谱带宽是复数据的两倍，直接对 1 倍采样的幅度做互相关会有"像素锁定"偏差。"""
    h, w = z.shape
    if ovs_a == 1 and ovs_r == 1:
        a = np.abs(z)
        return a[crop:h - crop, crop:w - crop] if crop else a
    zf = _fft2(z)
    zp = np.zeros((h * ovs_a, w * ovs_r), np.complex64)
    hh, ww = h // 2, w // 2
    zp[:hh, :ww] = zf[:hh, :ww]
    zp[:hh, zp.shape[1] - (w - ww):] = zf[:hh, ww:]
    zp[zp.shape[0] - (h - hh):, :ww] = zf[hh:, :ww]
    zp[zp.shape[0] - (h - hh):, zp.shape[1] - (w - ww):] = zf[hh:, ww:]
    a = np.abs(_ifft2(zp)) * float(ovs_a * ovs_r)
    return a[crop * ovs_a:(h - crop) * ovs_a, crop * ovs_r:(w - crop) * ovs_r]


# ----------------------------------------------------------------------------
# 全局尺度: 多视对数幅度 + 掩膜归一化互相关
# ----------------------------------------------------------------------------
def _log_highpass(power, hp_box):
    valid = power > 0
    if not valid.any():
        raise InputError("影像全为 0，无法配准")
    floor = 1e-3 * float(np.median(power[valid]))
    lg = np.log(np.maximum(power, floor)).astype(np.float32)
    lg[~valid] = 0.0
    hp = lg - box_mean(lg, hp_box)
    hp[~valid] = 0.0
    return hp, valid.astype(np.float32)


def masked_ncc(a, ma, b, mb, min_pixels):
    """相关函数 c[s] = sum_x a(x) b(x+s) 的掩膜归一化版本 (Padfield)。
    返回 (ncc, lag_y, lag_x, fft_shape)，ncc 在 FFT 环绕索引上，重叠不足处为 -2。"""
    fs = (_fast_len(a.shape[0] + b.shape[0] - 1), _fast_len(a.shape[1] + b.shape[1] - 1))

    def spec(x):
        return _rfft2(np.ascontiguousarray(x, dtype=np.float32), fs)

    def xc(p, q):
        return _irfft2(np.conj(p) * q, fs)

    fam = spec(ma)
    fbm = spec(mb)
    n_raw = xc(fam, fbm)
    fa = spec(a * ma)
    sa = xc(fa, fbm)
    fb = spec(b * mb)
    sb = xc(fam, fb)
    sab = xc(fa, fb)
    saa = xc(spec(a * a * ma), fbm)
    sbb = xc(fam, spec(b * b * mb))
    n = np.maximum(n_raw, 1.0)
    num = sab - sa * sb / n
    va = np.maximum(saa - sa * sa / n, 0.0)
    vb = np.maximum(sbb - sb * sb / n, 0.0)
    den = np.sqrt(va * vb)
    ncc = np.full(num.shape, -2.0, np.float32)
    ok = (n_raw >= min_pixels) & (den > 1e-6)
    ncc[ok] = (num[ok] / den[ok]).astype(np.float32)
    return ncc, fs


def _lag_index(i, n_b, fs):
    return i if i <= n_b - 1 else i - fs


def _parabola(m1, z, p1):
    den = m1 - 2.0 * z + p1
    return 0.0 if den == 0 else 0.5 * (m1 - p1) / den


def _second_peak(ncc, iy, ix, radius):
    """主峰 radius 像素之外的最高局部极大值 (环绕索引上取邻域)。"""
    valid = ncc > -1.5
    core = ncc
    is_max = valid.copy()
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                is_max &= core >= np.roll(np.roll(core, dy, axis=0), dx, axis=1)
    yy = np.abs(((np.arange(ncc.shape[0]) - iy + ncc.shape[0] // 2) % ncc.shape[0]) - ncc.shape[0] // 2)
    xx = np.abs(((np.arange(ncc.shape[1]) - ix + ncc.shape[1] // 2) % ncc.shape[1]) - ncc.shape[1] // 2)
    far = (yy[:, None] > radius) | (xx[None, :] > radius)
    cand = ncc[is_max & far]
    return float(cand.max()) if cand.size else float("-inf")


def global_stage(master, slave, target=1500, hp_box=15, looks=None, min_overlap=0.25):
    """全局尺度偏移 (SLC 像素)。返回 dict，失败时 ok=False。"""
    t0 = time.time()
    big = max(master.ny, master.nx, slave.ny, slave.nx)
    L = int(looks) if looks else max(1, int(math.ceil(big / float(target))))
    pa = multilook_power(master, L)
    pb = multilook_power(slave, L)
    a, ma = _log_highpass(pa, hp_box)
    b, mb = _log_highpass(pb, hp_box)
    t_read = time.time() - t0
    min_pix = min_overlap * min(a.size, b.size)
    ncc, fs = masked_ncc(a, ma, b, mb, min_pix)
    iy, ix = np.unravel_index(int(np.argmax(ncc)), ncc.shape)
    peak = float(ncc[iy, ix])
    res = dict(ok=False, looks=L, ncc=peak, time_s=0.0)
    if peak <= -1.5:
        res["reason"] = "重叠区域不足"
        return res
    ly = _lag_index(iy, b.shape[0], fs[0])
    lx = _lag_index(ix, b.shape[1], fs[1])
    # 相关面在 FFT 环绕索引上，滞后 0 位于索引 0，邻域须按模取
    ys = np.take(np.arange(fs[0]), np.arange(iy - 1, iy + 2), mode="wrap")
    xs = np.take(np.arange(fs[1]), np.arange(ix - 1, ix + 2), mode="wrap")
    nb = ncc[np.ix_(ys, xs)]
    sub_y = sub_x = 0.0
    if (nb > -1.5).all():
        sub_y = _parabola(nb[0, 1], nb[1, 1], nb[2, 1])
        sub_x = _parabola(nb[1, 0], nb[1, 1], nb[1, 2])
    valid = ncc > -1.5
    msk = valid.copy()
    for yy in range(-6, 7):
        msk[(iy + yy) % fs[0], np.take(np.arange(fs[1]), np.arange(ix - 6, ix + 7), mode="wrap")] = False
    side = ncc[msk]
    if side.size < 50:
        res["reason"] = "旁瓣样本不足"
        return res
    psr = (peak - float(side.mean())) / (float(side.std()) + 1e-12)
    second = _second_peak(ncc, iy, ix, radius=10)
    res.update(dr=(lx + sub_x) * L, da=(ly + sub_y) * L,
               psr=psr, second=second, peak_ratio=peak / second if second > 0 else float("inf"),
               time_read_s=t_read, time_s=time.time() - t0, tex_img=a)
    if peak < GLOBAL_MIN_NCC or psr < GLOBAL_MIN_PSR:
        res["reason"] = "全局相关峰不可信 (ncc=%.3f<%.2f 或 psr=%.1f<%.0f)，两幅影像可能不是同一区域/同轨重复观测" % (
            peak, GLOBAL_MIN_NCC, psr, GLOBAL_MIN_PSR)
        return res
    res["ok"] = True
    return res


# ----------------------------------------------------------------------------
# 模板匹配 (一个控制点)
# ----------------------------------------------------------------------------
def match_patch(tm, ss, up=16):
    """tm: 主影像 (h,w) 模板; ss: 从影像 (h+2R, w+2R) 搜索区。
    返回 dict(dy, dx, ncc, psr, edge) 或 None; dy/dx 为匹配位置相对窗口居中位置的亚像素位移。"""
    h, w = tm.shape
    H, W = ss.shape
    ly, lx = H - h + 1, W - w + 1
    if ly < 3 or lx < 3:
        return None
    n = float(h * w)
    tz = tm - tm.mean()
    tt = float(np.dot(tz.ravel().astype(np.float64), tz.ravel().astype(np.float64)))
    if tt < 1e-9:
        return None
    ss0 = ss - ss.mean()
    fs = (_fast_len(H), _fast_len(W))
    ft = _rfft2(tz.astype(np.float32), fs)
    fsr = _rfft2(ss0.astype(np.float32), fs)
    g = np.conj(ft) * fsr
    num = _irfft2(g, fs)[:ly, :lx].astype(np.float64)
    s1 = window_sums(ss0, h, w)
    s2 = window_sums(ss0.astype(np.float64) ** 2, h, w)
    var = s2 - s1 * s1 / n
    vfloor = max(1e-9, 1e-4 * float(np.median(var)))
    den = np.sqrt(tt * np.maximum(var, vfloor))
    ncc = num / den
    iy, ix = np.unravel_index(int(np.argmax(ncc)), ncc.shape)
    peak = float(ncc[iy, ix])
    msk = np.ones(ncc.shape, bool)
    msk[max(iy - 2, 0):iy + 3, max(ix - 2, 0):ix + 3] = False
    side = ncc[msk]
    psr = (peak - float(side.mean())) / (float(side.std()) + 1e-12) if side.size > 8 else 0.0
    cy, cx = (ly - 1) / 2.0, (lx - 1) / 2.0
    edge = iy in (0, ly - 1) or ix in (0, lx - 1)
    out = dict(dy=iy - cy, dx=ix - cx, ncc=peak, psr=psr, edge=edge)
    if edge:
        return out

    # 互相关分子用 DFT 局部上采样 (±1 像素, 步长 1/up)，分母由整数滞后双线性插值
    off = np.arange(-up, up + 1) / float(up)
    u = iy + off
    v = ix + off
    fy = np.fft.fftfreq(fs[0])
    fx = np.fft.rfftfreq(fs[1])
    wx = np.where((fx == 0) | (fx == 0.5), 1.0, 2.0)
    ey = np.exp(2j * np.pi * np.outer(u, fy)).astype(np.complex64)
    ex = np.exp(2j * np.pi * np.outer(fx, v)).astype(np.complex64)
    c_up = np.real(ey @ (g * wx[None, :]).astype(np.complex64) @ ex) / float(fs[0] * fs[1])
    u0 = np.clip(np.floor(u).astype(int), 0, ly - 2)
    v0 = np.clip(np.floor(v).astype(int), 0, lx - 2)
    tu = (u - u0)[:, None]
    tv = (v - v0)[None, :]
    d_up = ((1 - tu) * (1 - tv) * den[np.ix_(u0, v0)] + (1 - tu) * tv * den[np.ix_(u0, v0 + 1)]
            + tu * (1 - tv) * den[np.ix_(u0 + 1, v0)] + tu * tv * den[np.ix_(u0 + 1, v0 + 1)])
    n_up = c_up / d_up
    ky, kx = np.unravel_index(int(np.argmax(n_up)), n_up.shape)
    fy_, fx_ = float(off[ky]), float(off[kx])
    if 0 < ky < n_up.shape[0] - 1:
        fy_ += _parabola(n_up[ky - 1, kx], n_up[ky, kx], n_up[ky + 1, kx]) / up
    if 0 < kx < n_up.shape[1] - 1:
        fx_ += _parabola(n_up[ky, kx - 1], n_up[ky, kx], n_up[ky, kx + 1]) / up
    out.update(dy=iy + fy_ - cy, dx=ix + fx_ - cx, ncc=float(n_up[ky, kx]))
    return out


# ----------------------------------------------------------------------------
# 仿射偏移模型
# ----------------------------------------------------------------------------
class Affine:
    """dr = cr0 + cr1*x + cr2*y ;  da = ca0 + ca1*x + ca2*y (主影像绝对像素坐标)"""

    def __init__(self, cr, ca):
        self.cr = np.asarray(cr, dtype=np.float64)
        self.ca = np.asarray(ca, dtype=np.float64)

    @classmethod
    def translation(cls, dr, da):
        return cls([dr, 0.0, 0.0], [da, 0.0, 0.0])

    def predict(self, x, y):
        return (self.cr[0] + self.cr[1] * x + self.cr[2] * y,
                self.ca[0] + self.ca[1] * x + self.ca[2] * y)

    def prm_updates(self):
        r0, a0 = math.floor(self.cr[0]), math.floor(self.ca[0])
        return {"rshift": int(r0), "sub_int_r": "%.6f" % (self.cr[0] - r0),
                "stretch_r": "%.9g" % self.cr[1], "a_stretch_r": "%.9g" % self.cr[2],
                "ashift": int(a0), "sub_int_a": "%.6f" % (self.ca[0] - a0),
                "stretch_a": "%.9g" % self.ca[1], "a_stretch_a": "%.9g" % self.ca[2]}

    def as_dict(self):
        return dict(cr=[float(v) for v in self.cr], ca=[float(v) for v in self.ca])


def _residual_metrics(rr, ra):
    """只补充报告，不参与拟合/门禁；不足 4 点时残差统计为 null，而不是 0。"""
    metrics = dict(mad_definition="median(abs(residual - median(residual)))",
                   robust_sigma_definition="1.4826 * MAD (median-centered)",
                   rmse_definition="sqrt(mean(residual**2)); no mean subtraction",
                   metric_units="pixel", metric_min_points=4)
    for axis, residual in (("r", rr), ("a", ra)):
        mad = rmse = None
        if len(residual) >= 4:
            with np.errstate(over="ignore", invalid="ignore"):
                mad = float(np.median(np.abs(residual - np.median(residual))))
                rmse = float(np.sqrt(np.mean(residual ** 2)))
        metrics["mad_" + axis] = mad
        metrics["robust_sigma_" + axis] = 1.4826 * mad if mad is not None else None
        metrics["rmse_" + axis] = rmse
    return _clean(metrics)


def fit_affine_robust(x, y, dr, da, thr, rng, n_iter=1500, floor=0.2):
    """RANSAC(MSAC 打分) + Tukey IRLS。返回 (Affine, inlier_mask, info) 或 None。"""
    n = len(x)
    if n < 6:
        return None
    xc, yc = float(np.mean(x)), float(np.mean(y))
    sc = max(float(np.ptp(x)), float(np.ptp(y)), 1.0) / 2.0
    xn, yn = (x - xc) / sc, (y - yc) / sc
    X = np.column_stack([np.ones(n), xn, yn])
    D = np.column_stack([dr, da])

    idx = rng.integers(0, n, size=(n_iter, 3))
    A = X[idx]
    ok = np.abs(np.linalg.det(A)) > 0.05
    if ok.sum() < 20:
        return None
    coefs = np.linalg.solve(A[ok], D[idx[ok]])
    pred = np.einsum("nj,kjc->knc", X, coefs.astype(np.float64))
    r2 = ((D[None, :, :] - pred) ** 2).sum(-1)
    score = np.minimum(r2, thr * thr).sum(1)
    coef = coefs[int(np.argmin(score))]

    for _ in range(10):
        res = np.hypot(*(D - X @ coef).T)
        s = max(float(np.median(res)) / 1.1774, floor / 3.0)
        u = res / (4.685 * s)
        w = np.where(u < 1.0, (1.0 - u * u) ** 2, 0.0)
        if (w > 0).sum() < 4:
            break
        sw = np.sqrt(w)[:, None]
        new, *_ = np.linalg.lstsq(sw * X, sw * D, rcond=None)
        if np.max(np.abs(new - coef)) < 1e-9:
            coef = new
            break
        coef = new

    res = np.hypot(*(D - X @ coef).T)
    inl0 = res <= thr
    if inl0.sum() < 4:
        return None
    sigma = 1.4826 * float(np.median(np.abs((D - X @ coef)[inl0])))
    t_in = min(thr, max(3.0 * sigma, floor))
    inl = res <= t_in
    if inl.sum() < 4:
        return None
    coef, *_ = np.linalg.lstsq(X[inl], D[inl], rcond=None)
    res = np.hypot(*(D - X @ coef).T)
    inl = res <= t_in
    if inl.sum() < 4:
        return None
    coef, *_ = np.linalg.lstsq(X[inl], D[inl], rcond=None)
    resid = D[inl] - X[inl] @ coef

    cr = [coef[0, 0] - (coef[1, 0] * xc + coef[2, 0] * yc) / sc, coef[1, 0] / sc, coef[2, 0] / sc]
    ca = [coef[0, 1] - (coef[1, 1] * xc + coef[2, 1] * yc) / sc, coef[1, 1] / sc, coef[2, 1] / sc]
    model = Affine(cr, ca)
    pr, pa = model.predict(x, y)
    full = np.hypot(dr - pr, da - pa)
    info = dict(sigma=float(sigma), t_in=float(t_in),
                rms_r=float(np.sqrt(np.mean(resid[:, 0] ** 2))),
                rms_a=float(np.sqrt(np.mean(resid[:, 1] ** 2))))
    # 必须沿用最终 lstsq 的输入内点，不用返回前重新计算的 full <= t_in 掩膜。
    info.update(_residual_metrics(resid[:, 0], resid[:, 1]))
    info.update(rms_definition=info["rmse_definition"],
                metric_population="final least-squares input inliers; before returned mask recomputation",
                metric_n_points=int(inl.sum()),
                sigma_definition="1.4826 * median(abs(residual)); pooled r/a on inl0 before least-squares refinement; not median-centered MAD")
    return model, full <= t_in, info


# ----------------------------------------------------------------------------
# 点阵与逐点测量
# ----------------------------------------------------------------------------
def make_lattice(master, slave, model, W, R, grid, tex=None, phase=(0.5, 0.5), pick_frac=0.7):
    """主影像坐标下的点阵; 有纹理图时在每个格网内取纹理最强的位置。"""
    nxp, nyp = grid
    half = W // 2 + R + 4
    span_x, span_y = master.nx - 2 * half, master.ny - 2 * half
    if span_x <= 0 or span_y <= 0:
        return []
    cw, ch = span_x / float(nxp), span_y / float(nyp)
    pts = []
    if tex is not None:
        e, looks = tex
        e_med = float(np.median(e[e > 0])) if (e > 0).any() else 0.0
    for j in range(nyp):
        for i in range(nxp):
            x = half + (i + phase[0]) * cw
            y = half + (j + phase[1]) * ch
            if tex is not None:
                x0, x1 = x - 0.5 * pick_frac * cw, x + 0.5 * pick_frac * cw
                y0, y1 = y - 0.5 * pick_frac * ch, y + 0.5 * pick_frac * ch
                j0, j1 = max(int(y0 / looks), 0), min(int(y1 / looks) + 1, e.shape[0])
                i0, i1 = max(int(x0 / looks), 0), min(int(x1 / looks) + 1, e.shape[1])
                if j1 > j0 and i1 > i0:
                    sub = e[j0:j1, i0:i1]
                    k = int(np.argmax(sub))
                    jj, ii = divmod(k, sub.shape[1])
                    if sub[jj, ii] < 0.2 * e_med:
                        continue
                    x = (i0 + ii + 0.5) * looks
                    y = (j0 + jj + 0.5) * looks
            x = min(max(x, half), master.nx - half)
            y = min(max(y, half), master.ny - half)
            dr, da = model.predict(x, y)
            if not (half <= x + dr < slave.nx - half and half <= y + da < slave.ny - half):
                continue
            pts.append((x, y))
    return pts


def measure_points(master, slave, model, pts, W, R, prep="sqrt", ovs=(1, 1), zero_max=0.2):
    """逐点匹配，返回 ndarray 列: x, y, dr, da, ncc, psr, edge。
    ovs=(方位向, 距离向) 过采样倍数；>1 时先对复数块补零过采样再取幅度。"""
    rows = []
    osa, osr = ovs
    crop = 0 if (osa == 1 and osr == 1) else 3
    for (x, y) in pts:
        pdr, pda = model.predict(x, y)
        ox, oy = int(round(pdr)), int(round(pda))
        mx0, my0 = int(round(x)) - W // 2, int(round(y)) - W // 2
        zt = master.read_cplx(my0 - crop, mx0 - crop, W + 2 * crop, W + 2 * crop)
        if np.count_nonzero(zt == 0) > zero_max * zt.size:
            continue
        zs = slave.read_cplx(my0 + oy - R - crop, mx0 + ox - R - crop, W + 2 * R + 2 * crop, W + 2 * R + 2 * crop)
        if np.count_nonzero(zs == 0) > zero_max * zs.size:
            continue
        tm = prep_amp(amp_upsampled(zt, osa, osr, crop), prep)
        ss = prep_amp(amp_upsampled(zs, osa, osr, crop), prep)
        r = match_patch(tm, ss)
        if r is None:
            continue
        rows.append((mx0 + (W - 1) / 2.0, my0 + (W - 1) / 2.0,
                     ox + r["dx"] / osr, oy + r["dy"] / osa, r["ncc"], r["psr"], 1.0 if r["edge"] else 0.0))
    return np.array(rows, dtype=np.float64).reshape(-1, 7)


# ----------------------------------------------------------------------------
# 逐级估计
# ----------------------------------------------------------------------------
def _stage_stats(name, st, pts, meas, valid, inl, info, t):
    stats = dict(name=name, window=st["W"], search=st["R"], n_points=len(pts), n_measured=int(len(meas)),
                 n_valid=int(valid.sum()), n_inliers=int(inl.sum()) if inl is not None else 0,
                 inlier_frac=float(inl.sum()) / max(int(valid.sum()), 1) if inl is not None else 0.0,
                 rms_r=info.get("rms_r") if info else None, rms_a=info.get("rms_a") if info else None,
                 time_s=round(t, 2))
    stats.update(_residual_metrics([], []))
    stats.update(rms_definition=stats["rmse_definition"],
                 metric_population="final least-squares input inliers; before returned mask recomputation",
                 metric_n_points=0)
    if info:
        stats.update(info)
    return stats


def estimate(master_prm, slave_prm, schedule=DEFAULT_SCHEDULE, prep="sqrt", ncc_min=0.08,
             global_target=1500, seed=0, ovs_override=None):
    rng = np.random.default_rng(seed)
    t_all = time.time()
    master = SLC(master_prm)
    slave = SLC(slave_prm)
    report = dict(module="ms_register", version=__version__,
                  master=dict(prm=os.path.abspath(master_prm), slc=master.path, nx=master.nx, ny=master.ny),
                  slave=dict(prm=os.path.abspath(slave_prm), slc=slave.path, nx=slave.nx, ny=slave.ny),
                  stages=[], warnings=[])

    log("全局尺度: 多视幅度互相关 (不使用 PRM 中的初始偏移)")
    g = global_stage(master, slave, target=global_target)
    report["global"] = {k: v for k, v in g.items() if k not in ("tex_img", "tex_valid")}
    if not g["ok"]:
        report["status"] = "failed"
        report["fail_reason"] = "全局尺度失败: %s" % g.get("reason", "未知")
        return report, None, None
    log("  全局偏移 dr=%+.2f da=%+.2f (ncc=%.3f psr=%.1f looks=%d, %.1fs)"
        % (g["dr"], g["da"], g["ncc"], g["psr"], g["looks"], g["time_s"]))

    L = g["looks"]
    hp_img = g["tex_img"]
    model = Affine.translation(g["dr"], g["da"])
    last_ok = None
    final_points = None
    for st in schedule:
        t0 = time.time()
        W, R = st["W"], st["R"]
        name = "S%d" % W
        k = max(3, int(round(W / float(L))))
        tex = (box_mean(hp_img * hp_img, k), L)
        pts = make_lattice(master, slave, model, W, R, st["grid"], tex=tex)
        meas = measure_points(master, slave, model, pts, W, R, prep=prep, ovs=ovs_override or st["ovs"])
        valid = (meas[:, 4] >= ncc_min) & (meas[:, 6] < 0.5) if len(meas) else np.zeros(0, bool)
        fit = None
        if valid.sum() >= 6:
            m = meas[valid]
            fit = fit_affine_robust(m[:, 0], m[:, 1], m[:, 2], m[:, 3], st["thr"], rng)
        if fit is None:
            stats = _stage_stats(name, st, pts, meas, valid, None, None, time.time() - t0)
            stats.update(accepted=False, reason="有效点不足或几何拟合失败")
            report["stages"].append(stats)
            log("  %s: 有效点 %d/%d，几何拟合失败" % (name, int(valid.sum()), len(pts)))
            if last_ok is not None:
                report["warnings"].append("%s 尺度几何拟合失败，沿用 %s 的结果" % (name, last_ok["name"]))
            break
        new_model, inl_v, info = fit
        stats = _stage_stats(name, st, pts, meas, valid, inl_v, info, time.time() - t0)
        ok = (stats["n_inliers"] >= st["min_inl"]) and (stats["inlier_frac"] >= st["min_frac"])
        stats["accepted"] = bool(ok)
        cdr, cda = new_model.predict(master.nx / 2.0, master.ny / 2.0)
        stats["center_dr"], stats["center_da"] = float(cdr), float(cda)
        stats["model"] = new_model.as_dict()
        report["stages"].append(stats)
        log("  %s: 点 %d, 有效 %d, 内点 %d (%.0f%%), rmse=(%.3f,%.3f)px, robust_sigma=(%.3f,%.3f)px, 中心偏移=(%+.2f,%+.2f), %.1fs %s"
            % (name, len(pts), stats["n_valid"], stats["n_inliers"], 100 * stats["inlier_frac"],
               info["rms_r"], info["rms_a"], info["robust_sigma_r"], info["robust_sigma_a"],
               cdr, cda, stats["time_s"], "通过" if ok else "未通过"))
        if not ok:
            if last_ok is not None:
                report["warnings"].append("%s 尺度内点 %d (%.0f%%) 未达标，沿用 %s 的结果" % (
                    name, stats["n_inliers"], 100 * stats["inlier_frac"], last_ok["name"]))
            break
        model = new_model
        last_ok = stats
        m = meas[valid]
        final_points = np.column_stack([m, inl_v.astype(float)])
        final_points = (st, final_points)

    report["time_total_s"] = round(time.time() - t_all, 2)
    if last_ok is None:
        report["status"] = "failed"
        report["fail_reason"] = "没有任何一级通过几何一致性检验"
        return report, None, None
    report["status"] = "estimated"
    report["final_stage"] = last_ok["name"]
    report["model"] = model.as_dict()
    report["prm"] = {k: v for k, v in model.prm_updates().items()}
    cdr, cda = model.predict(master.nx / 2.0, master.ny / 2.0)
    report["center_offset"] = dict(range=float(cdr), azimuth=float(cda))
    return report, model, final_points


def check_model(master, slave, model, W=128, R=8, grid=(14, 14), prep="sqrt", ncc_min=0.08,
                tex=None, phase=(0.25, 0.75), ovs=CHECK_OVS):
    """在独立点阵上测量残余偏移 (相对 model 的预测)。返回统计 dict。"""
    pts = make_lattice(master, slave, model, W, R, grid, tex=tex, phase=phase)
    meas = measure_points(master, slave, model, pts, W, R, prep=prep, ovs=ovs)
    n_pts = len(pts)
    ok = (meas[:, 4] >= ncc_min) & (meas[:, 6] < 0.5) if len(meas) else np.zeros(0, bool)
    m = meas[ok]
    pr, pa = model.predict(m[:, 0], m[:, 1])
    rr, ra = m[:, 2] - pr, m[:, 3] - pa
    n_within = int(np.count_nonzero(np.hypot(rr, ra) <= 0.3))
    stats = dict(n_points=n_pts, n_valid=int(len(m)), valid_frac=len(m) / max(n_pts, 1),
                 n_within_0p3=n_within,
                 frac_within_0p3_all_candidates=n_within / n_pts if n_pts else None,
                 within_0p3_definition="hypot(residual_r, residual_a) <= 0.3 pixel; no centering",
                 frac_within_0p3_definition="n_within_0p3 / n_valid (legacy; requires >=4 valid points)",
                 frac_within_0p3_all_candidates_definition="n_within_0p3 / n_points; null if no candidates",
                 metric_population="valid independent check points (ncc >= ncc_min and edge < 0.5)",
                 metric_n_points=int(len(m)))
    stats.update(_residual_metrics(rr, ra))
    stats["rms_definition"] = stats["robust_sigma_definition"] + "; legacy rms_r/rms_a are not RMSE"
    if len(m) < 4:
        return stats
    mad = lambda v: 1.4826 * float(np.median(np.abs(v - np.median(v))))
    stats.update(median_r=float(np.median(rr)), median_a=float(np.median(ra)),
                 rms_r=mad(rr), rms_a=mad(ra),
                 max_abs_r=float(np.max(np.abs(rr - np.median(rr)))), max_abs_a=float(np.max(np.abs(ra - np.median(ra)))),
                 frac_within_0p3=float(np.mean(np.hypot(rr, ra) <= 0.3)),
                 residuals=np.column_stack([m[:, 0], m[:, 1], rr, ra, m[:, 4]]))
    return stats


# ----------------------------------------------------------------------------
# 门限
# ----------------------------------------------------------------------------
GATE_DEFAULTS = dict(
    min_stage="S256",      # 至少需要通过到该尺度
    min_inliers=20,
    min_inlier_frac=0.30,
    max_rms=0.35,          # 最终一级内点残差 (像素)
    min_coverage=5,        # 3x3 分区中至少有内点的分区数
    check_median=0.12,     # 独立点复核: |残余中位数|
    check_rms=0.35,        # 独立点复核: 稳健 rms
    check_valid_frac=0.30,
)


def coverage_cells(master, pts_xy):
    cx = np.clip((pts_xy[:, 0] / master.nx * 3).astype(int), 0, 2)
    cy = np.clip((pts_xy[:, 1] / master.ny * 3).astype(int), 0, 2)
    return int(len({(a, b) for a, b in zip(cx, cy)}))


def evaluate_gate(report, master, final_points, check, gate):
    reasons = []
    st = [s for s in report["stages"] if s.get("accepted")]
    if not st:
        return False, ["没有任何一级通过检验"]
    order = ["S512", "S256", "S128", "S64"]
    if order.index(st[-1]["name"]) < order.index(gate["min_stage"]):
        reasons.append("仅通过到 %s，低于要求的 %s" % (st[-1]["name"], gate["min_stage"]))
    last = st[-1]
    if last["n_inliers"] < gate["min_inliers"]:
        reasons.append("最终内点数 %d < %d" % (last["n_inliers"], gate["min_inliers"]))
    if last["inlier_frac"] < gate["min_inlier_frac"]:
        reasons.append("最终内点比例 %.2f < %.2f" % (last["inlier_frac"], gate["min_inlier_frac"]))
    if max(last["rms_r"], last["rms_a"]) > gate["max_rms"]:
        reasons.append("最终内点残差 rms=(%.3f,%.3f) > %.2f px" % (last["rms_r"], last["rms_a"], gate["max_rms"]))
    if final_points is not None:
        fp = final_points[1]
        cov = coverage_cells(master, fp[fp[:, 7] > 0.5][:, :2])
        report["coverage_cells"] = cov
        if cov < gate["min_coverage"]:
            reasons.append("内点只覆盖 %d/9 个分区 (< %d)，模型外推风险高" % (cov, gate["min_coverage"]))
    if check is not None:
        if check.get("n_valid", 0) < 8 or check.get("valid_frac", 0) < gate["check_valid_frac"]:
            reasons.append("独立点复核有效点不足 (%s)" % check.get("n_valid"))
        else:
            if max(abs(check["median_r"]), abs(check["median_a"])) > gate["check_median"]:
                reasons.append("独立点复核残余中位数 (%.3f,%.3f) > %.2f px"
                               % (check["median_r"], check["median_a"], gate["check_median"]))
            if max(check["rms_r"], check["rms_a"]) > gate["check_rms"]:
                reasons.append("独立点复核 rms (%.3f,%.3f) > %.2f px"
                               % (check["rms_r"], check["rms_a"], gate["check_rms"]))
    return (len(reasons) == 0), reasons


# ----------------------------------------------------------------------------
# 输出
# ----------------------------------------------------------------------------
def _clean(o):
    """转成可被严格 JSON 解析器读取的纯 Python 对象 (numpy 标量转原生类型, NaN/Inf 置 None)。"""
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items() if k != "residuals"}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
        return None
    return o


def write_report(path, report):
    with open(path, "w") as f:
        json.dump(_clean(report), f, ensure_ascii=False, indent=2)


def write_points(path, final_points):
    st, pts = final_points
    with open(path, "w") as f:
        f.write("# x_master dr y_master da ncc*100 inlier  (window=%d)  从影像坐标 = 主影像坐标 + 偏移\n" % st["W"])
        for p in pts:
            f.write(" %.1f %.3f %.1f %.3f %.2f %d\n" % (p[0], p[2], p[1], p[3], 100.0 * p[4], int(p[7])))


def write_text_log(path, report, master_prm, slave_prm):
    lines = ["# ms_register %s  多尺度配准结果" % __version__,
             "master_prm = %s" % master_prm, "slave_prm = %s" % slave_prm,
             "status = %s" % report.get("status")]
    if "center_offset" in report:
        lines += ["range_offset = %.3f pixel" % report["center_offset"]["range"],
                  "azimuth_offset = %.3f pixel" % report["center_offset"]["azimuth"],
                  "# 偏移取影像中心处, 从影像坐标 = 主影像坐标 + 偏移"]
        for k in REG_KEYS:
            lines.append("%s = %s" % (k, report["prm"][k]))
    g = report.get("global", {})
    if g.get("ok"):
        lines.append("global: dr=%.2f da=%.2f ncc=%.3f psr=%.1f" % (g["dr"], g["da"], g["ncc"], g["psr"]))
    for s in report.get("stages", []):
        lines.append("stage %s: points=%d valid=%d inliers=%d (%.0f%%) rmse=(%s,%s) robust_sigma=(%s,%s) accepted=%s"
                     % (s["name"], s["n_points"], s["n_valid"], s["n_inliers"], 100 * s["inlier_frac"],
                        "%.3f" % s["rms_r"] if s.get("rms_r") is not None else "-",
                        "%.3f" % s["rms_a"] if s.get("rms_a") is not None else "-",
                        "%.3f" % s["robust_sigma_r"] if s.get("robust_sigma_r") is not None else "-",
                        "%.3f" % s["robust_sigma_a"] if s.get("robust_sigma_a") is not None else "-", s["accepted"]))
    for w in report.get("warnings", []):
        lines.append("warning: %s" % w)
    for key in ("check", "gate", "resample", "verify", "outputs"):
        if key in report:
            lines.append("%s = %s" % (key, json.dumps(_clean(report[key]), ensure_ascii=False)))
    if report.get("fail_reason"):
        lines.append("fail_reason = %s" % report["fail_reason"])
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# 命令行
# ----------------------------------------------------------------------------
def _build_schedule(scales):
    if not scales:
        return DEFAULT_SCHEDULE
    by_w = {s["W"]: s for s in DEFAULT_SCHEDULE}
    out = []
    for w in scales:
        if w not in by_w:
            raise InputError("不支持的尺度 %d (可选 %s)" % (w, sorted(by_w, reverse=True)))
        out.append(by_w[w])
    return tuple(out)


def run_estimate(args):
    """估计 + 质量门限 + 落盘 (PRM/点位/报告/日志)。返回 (退出码, 报告 dict)。"""
    master = SLC(args.master)
    slave = SLC(args.slave)
    schedule = _build_schedule([int(v) for v in args.scales.split(",")] if args.scales else None)
    ovs = tuple(int(v) for v in args.ovs.split("x")) if args.ovs else None
    report, model, final_points = estimate(args.master, args.slave, schedule=schedule, prep=args.prep,
                                           global_target=args.global_target, seed=args.seed, ovs_override=ovs)
    check = None
    if model is not None:
        t0 = time.time()
        check = check_model(master, slave, model, prep=args.prep, ovs=ovs or CHECK_OVS)
        report["check"] = _clean(check)
        report["check_time_s"] = round(time.time() - t0, 2)
        log("  独立点复核(模型预测): 有效 %s/%s, 残余中位数=(%s,%s), robust_sigma=(%s,%s)" % (
            check.get("n_valid"), check.get("n_points"),
            "%+.3f" % check["median_r"] if "median_r" in check else "-",
            "%+.3f" % check["median_a"] if "median_a" in check else "-",
            "%.3f" % check["rms_r"] if "rms_r" in check else "-",
            "%.3f" % check["rms_a"] if "rms_a" in check else "-"))
    gate = dict(GATE_DEFAULTS)
    passed, reasons = (False, [report.get("fail_reason", "估计失败")]) if model is None else \
        evaluate_gate(report, master, final_points, check, gate)
    report["gate"] = dict(passed=passed, reasons=reasons, thresholds=gate)
    report["status"] = "passed" if passed else "failed"
    if model is not None and (passed or args.write_on_fail):
        out_prm = args.output or args.slave
        updates = model.prm_updates()
        if os.path.dirname(os.path.abspath(out_prm)) != os.path.dirname(slave.prm_path):
            updates["SLC_file"] = slave.path  # 输出 PRM 换了目录, 原 PRM 里的相对 SLC 路径会失效
        write_prm(args.slave, out_prm, updates)
        report["output_prm"] = os.path.abspath(out_prm)
        log("已写入配准参数: %s" % out_prm)
    if final_points is not None and args.points:
        write_points(args.points, final_points)
    write_report(args.report, report)
    write_text_log(args.log, report, args.master, args.slave)
    if passed:
        log("配准估计通过 (总耗时 %.1fs)" % report["time_total_s"])
        return EXIT_OK, report
    log("配准质量门限未通过: %s" % "; ".join(reasons))
    return EXIT_GATE, report


def cmd_estimate(args):
    return run_estimate(args)[0]


def verify_registered(master, slave, grid=(14, 14), prep="sqrt"):
    """主影像与已配准 (已重采样) 的从影像应当逐像素对齐: 在独立点阵上测残余偏移。
    master/slave 为尺寸相同的 SLC 对象。返回 (check dict, gate dict)。"""
    if (slave.nx, slave.ny) != (master.nx, master.ny):
        raise InputError("已配准从影像尺寸 %dx%d 与主影像 %dx%d 不一致" % (slave.nx, slave.ny, master.nx, master.ny))
    L = max(1, int(math.ceil(max(master.ny, master.nx) / 1500.0)))
    hp, _ = _log_highpass(multilook_power(master, L), 15)
    tex = (box_mean(hp * hp, max(3, int(round(128.0 / L)))), L)
    chk = check_model(master, slave, Affine.translation(0.0, 0.0), W=128, R=8, grid=grid, prep=prep, tex=tex)
    gate = dict(GATE_DEFAULTS)
    reasons = []
    if chk.get("n_valid", 0) < 8 or chk.get("valid_frac", 0) < gate["check_valid_frac"]:
        reasons.append("有效复核点不足 (%s/%s)" % (chk.get("n_valid"), chk.get("n_points")))
    else:
        if max(abs(chk["median_r"]), abs(chk["median_a"])) > gate["check_median"]:
            reasons.append("残余中位数 (%.3f,%.3f) > %.2f px" % (chk["median_r"], chk["median_a"], gate["check_median"]))
        if max(chk["rms_r"], chk["rms_a"]) > gate["check_rms"]:
            reasons.append("残余 rms (%.3f,%.3f) > %.2f px" % (chk["rms_r"], chk["rms_a"], gate["check_rms"]))
    if "median_r" in chk:
        log("复核: 有效 %d/%d, 残余中位数=(%+.3f,%+.3f)px, robust_sigma=(%.3f,%.3f)px, 0.3px内占有效点 %.0f%%" % (
            chk["n_valid"], chk["n_points"], chk["median_r"], chk["median_a"], chk["rms_r"], chk["rms_a"],
            100 * chk["frac_within_0p3"]))
    return chk, dict(passed=not reasons, reasons=reasons, thresholds=gate)


def cmd_verify(args):
    master = SLC(args.master)
    slave = SLC(args.slave, slc_path=args.slave_slc)
    t0 = time.time()
    chk, gate = verify_registered(master, slave, grid=args.grid, prep=args.prep)
    rep = dict(module="ms_register", version=__version__, mode="verify", master=master.path, slave=slave.path,
               check=chk, gate=gate, time_s=round(time.time() - t0, 2))
    if args.report:
        write_report(args.report, rep)
    if not gate["passed"]:
        log("复核未通过: %s" % "; ".join(gate["reasons"]))
        return EXIT_GATE
    log("复核通过")
    return EXIT_OK


class ResampError(Exception):
    pass


def run_resamp(exe, master, slave, est_prm, out_slc, out_prm, interp, prm_slc_file):
    """用 GMTSAR resamp 把从影像重采样到主影像网格。

    resamp 按 PRM 里的 SLC_file 相对当前目录取文件，所以在临时目录里放两份 SLC_file 改成绝对路径的 PRM，
    主、从影像不必在同一目录。输出先写到临时目录，成功后再改名，失败时不会留下半成品。
    resamp 输出的 PRM 沿用从影像原来的 SLC_file (GMTSAR 约定: 重采样后的 SLC 顶替原文件名)，
    这里改写成 prm_slc_file，避免出现 "配准后的 PRM 指向未配准 SLC" 的陷阱。"""
    tmp = tempfile.mkdtemp(prefix=".ms_register_", dir=os.path.dirname(out_slc))
    try:
        mprm, sprm = os.path.join(tmp, "m.PRM"), os.path.join(tmp, "s.PRM")
        write_prm(master.prm_path, mprm, {"SLC_file": master.path})
        write_prm(est_prm, sprm, {"SLC_file": slave.path})
        t_slc, t_prm = os.path.join(tmp, "reg.SLC"), os.path.join(tmp, "reg.PRM")
        t0 = time.time()
        proc = subprocess.run([exe, mprm, sprm, t_prm, t_slc, str(interp)], cwd=tmp,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        tail = "\n".join(proc.stdout.strip().splitlines()[-8:])
        if proc.returncode != 0 or not os.path.isfile(t_slc) or not os.path.isfile(t_prm):
            raise ResampError("resamp 退出码 %d %s" % (proc.returncode, tail))
        need = 4 * master.nx * master.ny
        if os.path.getsize(t_slc) < need:
            raise ResampError("resamp 输出 SLC 过小 (%d < %d 字节) %s" % (os.path.getsize(t_slc), need, tail))
        write_prm(t_prm, out_prm, {"SLC_file": prm_slc_file})
        os.replace(t_slc, out_slc)
        return dict(tool=exe, interp=int(interp), time_s=round(time.time() - t0, 2),
                    slc=out_slc, prm=out_prm, slc_file_in_prm=prm_slc_file)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _remove_quiet(*paths):
    for p in paths:
        try:
            os.remove(p)
        except OSError:
            pass


def cmd_register(args):
    t_all = time.time()
    master = SLC(args.master)
    slave = SLC(args.slave)
    exe = args.resamp or shutil.which("resamp")
    if not exe or not os.access(exe, os.X_OK):
        raise InputError("找不到 GMTSAR resamp (把 GMTSAR bin 目录加入 PATH, 或用 --resamp 指定)")
    outdir = os.path.abspath(args.outdir or os.path.dirname(slave.prm_path))
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(slave.prm_path))[0]

    def target(given, default_name):
        return os.path.abspath(given) if given else os.path.join(outdir, default_name)

    out_slc = target(args.out_slc, stem + "_reg.SLC")
    out_prm = target(args.out_prm, stem + "_reg.PRM")
    est_prm = target(args.est_prm, stem + "_est.PRM")
    for p in (out_slc, out_prm, est_prm):
        os.makedirs(os.path.dirname(p), exist_ok=True)
    inputs = {os.path.realpath(p) for p in (master.path, slave.path, master.prm_path, slave.prm_path)}
    outputs = [os.path.realpath(p) for p in (out_slc, out_prm, est_prm)]
    if len(set(outputs)) < 3 or inputs & set(outputs):
        raise InputError("输出文件 (配准后 SLC/PRM、估计 PRM) 不能互相重名，也不能覆盖输入的 SLC/PRM")
    if os.path.dirname(out_slc) == os.path.dirname(out_prm):
        prm_slc_file = args.prm_slc_file or os.path.basename(out_slc)
    else:
        prm_slc_file = args.prm_slc_file or out_slc
    need = 4 * master.nx * master.ny
    free = shutil.disk_usage(os.path.dirname(out_slc)).free
    if free < need * 1.05:
        raise InputError("输出目录剩余空间不足: 需要约 %.1f GB, 剩余 %.1f GB" % (need / 1e9, free / 1e9))

    ns = copy.copy(args)
    ns.output, ns.write_on_fail = est_prm, False
    ns.report = target(args.report, "coreg_report.json")
    ns.log = target(args.log, "coreg_offsets.log")
    ns.points = target(args.points, "coreg_points.dat")
    # 上一次运行遗留的结果在本次失败时不能被下游误当成这次的结果
    _remove_quiet(out_slc, out_prm, est_prm)

    rc, report = run_estimate(ns)
    report["mode"] = "register"
    report["outputs"] = dict(registered_slc=out_slc, registered_prm=out_prm, estimated_prm=est_prm,
                             report=ns.report, offsets_log=ns.log, points=ns.points)
    def finish(rc_, status, reason=None):
        report["status"] = status
        if reason:
            report["fail_reason"] = reason
        report["time_total_s"] = round(time.time() - t_all, 2)
        write_report(ns.report, report)
        write_text_log(ns.log, report, args.master, args.slave)
        return rc_

    if rc != EXIT_OK:
        return finish(rc, "failed")

    # run_estimate 已经把 "估计通过" 的报告写到了磁盘; 此后任何一步出错 (含意外异常) 都必须改写成失败,
    # 否则下游读到的是一份 status=passed 但并没有配准结果的报告
    try:
        log("重采样: resamp (interp=%d)" % args.interp)
        try:
            report["resample"] = run_resamp(exe, master, slave, est_prm, out_slc, out_prm, args.interp, prm_slc_file)
        except ResampError as e:
            log("重采样失败: %s" % e)
            _remove_quiet(out_slc, out_prm)
            return finish(EXIT_INTERNAL, "failed", "重采样失败: %s" % e)
        log("  已生成 %s (%.1fs)" % (out_slc, report["resample"]["time_s"]))

        if not args.skip_verify:
            t0 = time.time()
            chk, vgate = verify_registered(master, SLC(out_prm, slc_path=out_slc), grid=args.grid, prep=args.prep)
            report["verify"] = dict(check=_clean(chk), gate=vgate, time_s=round(time.time() - t0, 2))
            if not vgate["passed"]:
                reason = "重采样结果复核未通过: %s" % "; ".join(vgate["reasons"])
                log(reason)
                if not args.keep_failed:
                    _remove_quiet(out_slc, out_prm)
                return finish(EXIT_GATE, "failed", reason)
            log("重采样结果复核通过")
    except Exception as e:  # noqa: BLE001
        _remove_quiet(out_slc, out_prm)
        finish(EXIT_INTERNAL, "failed", "意外错误: %s: %s" % (type(e).__name__, e))
        raise
    log("配准完成: %s (总耗时 %.1fs)" % (out_slc, time.time() - t_all))
    return finish(EXIT_OK, "passed")


def main(argv=None):
    global _WORKERS, _QUIET
    ap = argparse.ArgumentParser(description="多尺度 SAR SLC 配准估计 (GMTSAR PRM 兼容)")
    ap.add_argument("--version", action="version", version="ms_register %s" % __version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    pe = sub.add_parser("estimate", help="估计配准参数并写回从影像 PRM")
    pe.add_argument("master")
    pe.add_argument("slave")
    pe.add_argument("-o", "--output", help="写入的 PRM (默认原地更新从影像 PRM)")
    pe.add_argument("--report", default="coreg_report.json")
    pe.add_argument("--log", default="coreg_offsets.log")
    pe.add_argument("--points", default="coreg_points.dat")
    pe.add_argument("--scales", help="逗号分隔的窗口序列, 默认 512,256,128,64")
    pe.add_argument("--prep", choices=["sqrt", "amp", "log"], default="sqrt")
    pe.add_argument("--global-target", type=int, default=1500, help="全局尺度多视后长边像素数")
    pe.add_argument("--ovs", help="强制所有尺度使用的 FFT 过采样倍数, 格式 方位x距离, 如 2x2 (默认 512/256 尺度 1x1, 128/64 尺度 2x2)")
    pe.add_argument("--seed", type=int, default=0)
    pe.add_argument("--write-on-fail", action="store_true", help="门限未通过时仍写回 PRM (仅调试)")
    pe.add_argument("--threads", type=int, default=1)
    pe.add_argument("-q", "--quiet", action="store_true")

    pv = sub.add_parser("verify", help="对已配准从影像做独立点位复核")
    pv.add_argument("master")
    pv.add_argument("slave", help="已配准从影像的 PRM (尺寸须与主影像一致)")
    pv.add_argument("--slave-slc", help="已配准从影像 SLC 路径 (若 PRM 里的 SLC_file 指向的不是配准后的 SLC，须显式给出)")
    pv.add_argument("--report")
    pv.add_argument("--grid", type=lambda s: tuple(int(v) for v in s.split("x")), default=(14, 14))
    pv.add_argument("--prep", choices=["sqrt", "amp", "log"], default="sqrt")
    pv.add_argument("--threads", type=int, default=1)
    pv.add_argument("-q", "--quiet", action="store_true")

    pr = sub.add_parser("register", help="估计 -> resamp 重采样 -> 独立复核, 输出配准后的从影像 SLC/PRM 与报告")
    pr.add_argument("master", help="主影像 PRM")
    pr.add_argument("slave", help="从影像 PRM (不读取其中的初始偏移, 也不会被改写)")
    pr.add_argument("--outdir", help="输出目录 (默认: 从影像 PRM 所在目录)")
    pr.add_argument("--out-slc", help="配准后从影像 SLC (默认 <outdir>/<从影像名>_reg.SLC)")
    pr.add_argument("--out-prm", help="配准后从影像 PRM (默认 <outdir>/<从影像名>_reg.PRM)")
    pr.add_argument("--prm-slc-file", help="写入配准后 PRM 的 SLC_file 值 (默认: 配准后 SLC 的文件名; "
                                           "流水线若按 GMTSAR 习惯用符号链接改名，可显式指定)")
    pr.add_argument("--est-prm", help="估计出配准参数、尚未重采样的从影像 PRM (默认 <outdir>/<从影像名>_est.PRM)")
    pr.add_argument("--report", help="JSON 报告 (默认 <outdir>/coreg_report.json)")
    pr.add_argument("--log", help="文本偏移日志 (默认 <outdir>/coreg_offsets.log)")
    pr.add_argument("--points", help="最终一级控制点 (默认 <outdir>/coreg_points.dat)")
    pr.add_argument("--interp", type=int, choices=[1, 2, 3, 4], default=4,
                    help="resamp 插值方式: 1 最近邻, 2 双线性, 3 双二次, 4 双 sinc (默认 4)")
    pr.add_argument("--resamp", help="resamp 可执行文件 (默认在 PATH 中查找)")
    pr.add_argument("--scales", help="逗号分隔的窗口序列, 默认 512,256,128,64")
    pr.add_argument("--prep", choices=["sqrt", "amp", "log"], default="sqrt")
    pr.add_argument("--global-target", type=int, default=1500, help="全局尺度多视后长边像素数")
    pr.add_argument("--ovs", help="强制所有尺度使用的 FFT 过采样倍数, 格式 方位x距离, 如 2x2")
    pr.add_argument("--grid", type=lambda s: tuple(int(v) for v in s.split("x")), default=(14, 14),
                    help="重采样结果复核点阵 (默认 14x14)")
    pr.add_argument("--seed", type=int, default=0)
    pr.add_argument("--skip-verify", action="store_true", help="跳过对重采样结果的独立复核 (不建议)")
    pr.add_argument("--keep-failed", action="store_true", help="复核未通过时保留重采样输出 (默认删除, 防止下游误用)")
    pr.add_argument("--threads", type=int, default=1)
    pr.add_argument("-q", "--quiet", action="store_true")

    args = ap.parse_args(argv)
    _WORKERS = max(1, args.threads)
    _QUIET = args.quiet
    handlers = {"estimate": cmd_estimate, "verify": cmd_verify, "register": cmd_register}
    try:
        return handlers[args.cmd](args)
    except InputError as e:
        print("[ms_register] 输入错误: %s" % e, file=sys.stderr)
        return EXIT_INPUT
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        print("[ms_register] 内部错误: %s" % e, file=sys.stderr)
        return EXIT_INTERNAL


if __name__ == "__main__":
    sys.exit(main())
