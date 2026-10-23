#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DEM 精配准 (可选步骤): 估计雷达坐标 DEM (dem_ra.grd) 相对 SAR 影像的整体平移并修正。

GMTSAR 的 phasediff 直接按轨道/DEM 几何把 DEM 投到雷达坐标, DEM 与影像之间常剩一定量的整体错位
(DEM 自身的水平定位误差 + 轨道时间/距离偏差)。错位会让山区的地形相位去不干净, 也会让地理编码结果整体偏移。
本工具像 Gamma 的 offset_pwrm 那样, 用 DEM 坡度模拟雷达亮度, 与主影像幅度做分块互相关, 取稳健中位数平移,
把 dem_ra.grd 平移后输出; 同一个平移也能作用到地理编码用的 trans.dat。

子命令:
  estimate     --prm master.PRM --slc master.SLC --dem-ra dem_ra.grd --out dem_ra_coreg.grd --report dem_coreg.json
               质量门限不过 (起伏太小 / 互相关峰太弱 / 各块不一致 / 平移可忽略 / 平移后复核不通过) 时不写 --out,
               报告里 applied=false 并给出原因, 退出码仍为 0 (主流程继续用原 dem_ra.grd)。
  apply-trans  --trans trans.dat --report dem_coreg.json --out trans_coreg.dat
               报告里 applied=true 时, 把 trans.dat (SAT_llt2rat -bod, 每条 5 个 float64) 的距离向、方位向像元编号加上
               同一平移; 否则原样复制。

约定: dem_ra.grd 第 j 行 = 第 j 段方位线 (与 fast_dem2topo_ra.sh / topo_ra_grid.py 一致, 读取用 gmt grd2xyz -ZTLf),
第 i 列 = 第 i 段距离向像元。dem_ra 的高程是相对 PRM earth_radius 球面的。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy import ndimage as ndi

C_LIGHT = 299792458.0


def read_prm(path):
    d = {}
    with open(path, errors="ignore") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
    return d


def geometry(prm, dec_r, dec_a):
    """中间距离处的入射角、网格像元对应的地距 / 方位向间距 (米)。"""
    re_, h = float(prm["earth_radius"]), float(prm["SC_height"])
    fs, nr = float(prm["rng_samp_rate"]), int(float(prm["num_rng_bins"]))
    ds = C_LIGHT / (2.0 * fs)
    rng = float(prm["near_range"]) + 0.5 * nr * ds
    rs = re_ + h
    cos_t = (rs ** 2 - re_ ** 2 - rng ** 2) / (2.0 * re_ * rng)
    theta = float(np.arccos(np.clip(cos_t, -1.0, 1.0)))
    vg = float(prm["SC_vel"]) * re_ / rs
    return {"theta_deg": float(np.degrees(theta)),
            "ground_range_m_per_px": dec_r * ds / np.sin(theta),
            "azimuth_m_per_px": dec_a * vg / float(prm["PRF"])}


def grd_info(path):
    o = subprocess.check_output(["gmt", "grdinfo", "-C", path]).decode().split()
    full = subprocess.check_output(["gmt", "grdinfo", path]).decode()
    return {"nx": int(o[9]), "ny": int(o[10]), "w": o[1], "e": o[2], "s": o[3], "n": o[4], "dx": o[7], "dy": o[8],
            "pixel": "Pixel node" in full}


def read_grd(path):
    i = grd_info(path)
    b = subprocess.check_output(["gmt", "grd2xyz", "-ZTLf", path])
    return np.frombuffer(b, dtype="<f4").reshape(i["ny"], i["nx"])


def write_grd(arr, like, out):
    i = grd_info(like)
    raw = out + ".f32"
    np.ascontiguousarray(arr, dtype="<f4").tofile(raw)
    try:
        cmd = ["gmt", "xyz2grd", raw, "-ZTLf", "-R%s/%s/%s/%s" % (i["w"], i["e"], i["s"], i["n"]),
               "-I%s/%s" % (i["dx"], i["dy"]), "-G" + out]
        if i["pixel"]:
            cmd.append("-r")
        subprocess.check_call(cmd)
    finally:
        if os.path.exists(raw):
            os.remove(raw)


def amp_grid(slc, nrow, ncol, ny, nx, step=64):
    """主影像 SLC (int16 实/虚交错) 的幅度, 按块平均成 ny x nx 网格 (第 j 行 = 第 j 段方位线)。"""
    if os.path.getsize(slc) < nrow * ncol * 4:
        raise ValueError("SLC 文件比 PRM 给的 %dx%d 小" % (nrow, ncol))
    la, lr = nrow // ny, ncol // nx
    mm = np.memmap(slc, dtype="<i2", mode="r", shape=(nrow, ncol, 2))
    out = np.zeros((ny, nx), np.float32)
    for r0 in range(0, ny, step):
        r1 = min(ny, r0 + step)
        a = mm[r0 * la:r1 * la, :nx * lr].astype(np.float32)
        p = a[..., 0] ** 2 + a[..., 1] ** 2
        out[r0:r1] = np.sqrt(p.reshape(r1 - r0, la, nx, lr).mean(axis=(1, 3)))
    return out


def reduce_grid(a, g):
    if g == 1:
        return a
    ny, nx = a.shape
    return a.reshape(ny // g, g, nx // g, g).mean(axis=(1, 3))


def make_sim(dem, gr_spacing, theta_deg):
    """DEM 坡度 -> 模拟亮度 cos(入射角 - 距离向坡角): 迎坡亮、背坡暗。"""
    g = np.gradient(dem, axis=1) / gr_spacing
    return np.cos(np.deg2rad(theta_deg) - np.arctan(g))


def prep(x, sigma=25):
    return x - ndi.gaussian_filter(x, sigma)


def tile_shifts(sim, amp, ts, maxs=25):
    """分块互相关: 返回 (块中心行, 列, 行向平移, 列向平移, 相关峰) 数组, 平移 = 影像相对 DEM 模拟图的位移 (像元)。"""
    ny, nx = sim.shape
    maxs = min(maxs, ts // 4)
    win = np.outer(np.hanning(ts), np.hanning(ts))
    out = []
    for r0 in range(0, ny - ts + 1, ts // 2):
        for c0 in range(0, nx - ts + 1, ts // 2):
            a = amp[r0:r0 + ts, c0:c0 + ts].astype(np.float64)
            s = sim[r0:r0 + ts, c0:c0 + ts].astype(np.float64)
            if not (np.isfinite(a).all() and np.isfinite(s).all()) or s.std() < 1e-6 or a.std() < 1e-6:
                continue
            a = (a - a.mean()) * win
            s = (s - s.mean()) * win
            cc = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(s))).real
            cc /= (np.sqrt((a ** 2).sum() * (s ** 2).sum()) + 1e-12)
            cc = np.fft.fftshift(cc)
            cy = cx = ts // 2
            sub = cc[cy - maxs:cy + maxs + 1, cx - maxs:cx + maxs + 1]
            iy, ix = np.unravel_index(np.argmax(sub), sub.shape)
            pk = sub[iy, ix]

            def par(m, c, p):
                d = m - 2 * c + p
                return 0.0 if d == 0 else 0.5 * (m - p) / d
            dy = iy - maxs + (par(sub[iy - 1, ix], pk, sub[iy + 1, ix]) if 0 < iy < 2 * maxs else 0)
            dx = ix - maxs + (par(sub[iy, ix - 1], pk, sub[iy, ix + 1]) if 0 < ix < 2 * maxs else 0)
            out.append((r0 + ts / 2, c0 + ts / 2, dy, dx, pk))
    return np.array(out).reshape(-1, 5)


def pick_tile_size(shape):
    for ts in (512, 256, 128):
        if min(shape) >= 4 * ts:
            return ts
    return None


def mad(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def estimate_shift(dem, amp, gr_spacing, theta_deg, min_peak=0.2, min_tiles=6, max_mad_px=1.0,
                   min_relief_m=40.0, min_shift_px=0.3, max_shift_px=15.0):
    """返回 (shift 或 None, 报告 dict)。shift = (行向, 列向) 像元, 把 DEM 平移该量即与影像对齐。"""
    rep = {"applied": False, "reason": "", "n_tiles_used": 0}
    relief = float(np.percentile(dem, 95) - np.percentile(dem, 5))
    rep["relief_p5_p95_m"] = relief
    if relief < min_relief_m:
        rep["reason"] = "DEM 起伏太小 (p5-p95 %.0f m < %.0f m), 互相关不可靠" % (relief, min_relief_m)
        return None, rep
    ts = pick_tile_size(dem.shape)
    if ts is None:
        rep["reason"] = "网格太小, 凑不出分块"
        return None, rep
    rep["tile_size_px"] = ts
    amp = np.where(np.isfinite(amp) & (amp > 0), amp, 0.0)
    if not (amp > 0).any():
        rep["reason"] = "主影像幅度全为 0"
        return None, rep
    a = prep(np.log(amp / np.median(amp[amp > 0]) + 0.02))
    best = None
    for sign in (1, -1):
        t = tile_shifts(prep(sign * make_sim(dem, gr_spacing, theta_deg)), a, ts)
        if len(t) and (best is None or np.median(t[:, 4]) > best[0]):
            best = (np.median(t[:, 4]), sign, t)
    if best is None:
        rep["reason"] = "没有可用分块"
        return None, rep
    _, sign, t = best
    rep["brightness_sign"] = sign
    good = t[t[:, 4] > max(0.1, np.percentile(t[:, 4], 40))]
    rep["n_tiles_total"] = int(len(t))
    rep["n_tiles_used"] = int(len(good))
    if len(good) < min_tiles:
        rep["reason"] = "可用分块只有 %d 个 (< %d)" % (len(good), min_tiles)
        return None, rep
    dy, dx = float(np.median(good[:, 2])), float(np.median(good[:, 3]))
    my, mx = mad(good[:, 2]), mad(good[:, 3])
    pk = float(np.median(good[:, 4]))
    rep.update(shift_px_rows_cols=[dy, dx], tile_mad_px=[my, mx], median_peak=pk)
    if pk < min_peak:
        rep["reason"] = "互相关峰太弱 (中位 %.2f < %.2f)" % (pk, min_peak)
        return None, rep
    if max(my, mx) > max_mad_px:
        rep["reason"] = "各块平移不一致 (MAD %.2f / %.2f 像元 > %.1f), 不是整体平移" % (my, mx, max_mad_px)
        return None, rep
    if max(abs(dy), abs(dx)) > max_shift_px:
        rep["reason"] = "平移 %.1f / %.1f 像元超出上限 %.0f, 不自动修正" % (dy, dx, max_shift_px)
        return None, rep
    if max(abs(dy), abs(dx)) < min_shift_px:
        rep["reason"] = "平移可忽略 (%.2f / %.2f 像元), 不重采样" % (dy, dx)
        return None, rep
    # 平移后复核: 残余平移应接近 0, 且互相关峰不变差
    new = ndi.shift(dem, (dy, dx), order=3, mode="nearest")
    t1 = tile_shifts(prep(sign * make_sim(new, gr_spacing, theta_deg)), a, ts)
    g1 = t1[t1[:, 4] > max(0.1, np.percentile(t1[:, 4], 40))] if len(t1) else t1
    if len(g1) < min_tiles:
        rep["reason"] = "平移后复核的可用分块不足"
        return None, rep
    ry, rx = float(np.median(g1[:, 2])), float(np.median(g1[:, 3]))
    pk1 = float(np.median(g1[:, 4]))
    rep["residual_after_px_rows_cols"] = [ry, rx]
    rep["median_peak_after"] = pk1
    if max(abs(ry), abs(rx)) > max(0.5, 0.5 * max(abs(dy), abs(dx))) or pk1 < 0.9 * pk:
        rep["reason"] = "平移后复核不通过 (残余 %.2f / %.2f 像元, 峰 %.2f -> %.2f)" % (ry, rx, pk, pk1)
        return None, rep
    rep["applied"] = True
    rep["reason"] = "ok"
    rep["max_abs_dem_change_m"] = float(np.max(np.abs(new - dem)))
    return (dy, dx), rep


def cmd_estimate(a):
    prm = read_prm(a.prm)
    nrow, ncol = int(float(prm["num_lines"])), int(float(prm["num_rng_bins"]))
    rep = {"applied": False, "reason": ""}

    def finish(msg=None):
        if msg:
            rep["reason"] = msg
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)
        print("[dem_coreg] %s: %s" % ("已修正" if rep["applied"] else "跳过", rep["reason"]))
        return 0

    try:
        info = grd_info(a.dem_ra)
        ny, nx = info["ny"], info["nx"]
        if nrow % ny or ncol % nx:
            return finish("dem_ra 网格 %dx%d 不能整除 SLC %dx%d" % (nx, ny, ncol, nrow))
        dec_a, dec_r = nrow // ny, ncol // nx
        g = 4 // max(dec_a, dec_r) if max(dec_a, dec_r) < 4 else 1
        if g < 1 or ny % g or nx % g:
            g = 1
        geo = geometry(prm, dec_r * g, dec_a * g)
        rep["geometry"] = geo
        rep["grid"] = {"dem_ra_dec_rng_az": [dec_r, dec_a], "estimation_block": g}
        full = read_grd(a.dem_ra).astype(np.float32)
        full = np.where(np.isfinite(full), full, 0.0).astype(np.float32)
        dem = reduce_grid(full, g)
        amp = amp_grid(a.slc, nrow, ncol, ny // g, nx // g)
        shift, est = estimate_shift(dem, amp, geo["ground_range_m_per_px"], geo["theta_deg"], min_peak=a.min_peak,
                                    min_tiles=a.min_tiles, max_mad_px=a.max_mad, min_relief_m=a.min_relief,
                                    min_shift_px=a.min_shift, max_shift_px=a.max_shift)
    except Exception as exc:  # 可选步骤: 任何失败都只报告, 不中断主流程
        return finish("估计失败: %s" % exc)
    rep.update(est)
    if shift is None:
        return finish()
    dy, dx = shift[0] * g, shift[1] * g           # 网格像元
    rep["shift_grid_px_rows_cols"] = [dy, dx]
    rep["shift_m_azimuth_range"] = [dy * geo["azimuth_m_per_px"] / g, dx * geo["ground_range_m_per_px"] / g]
    rep["shift_slc_px"] = {"azimuth_lines": dy * dec_a, "range": dx * dec_r}
    write_grd(ndi.shift(full, (dy, dx), order=3, mode="nearest"), a.dem_ra, a.out)
    return finish()


def cmd_apply_trans(a):
    with open(a.report) as fh:
        rep = json.load(fh)
    n = os.path.getsize(a.trans)
    if n % 40:
        print("trans.dat 大小不是 40 字节的整数倍", file=sys.stderr)
        return 1
    if not rep.get("applied"):
        shutil.copyfile(a.trans, a.out)
        print("[dem_coreg] 报告为未修正, trans.dat 原样复制")
        return 0
    d_rng, d_az = rep["shift_slc_px"]["range"], rep["shift_slc_px"]["azimuth_lines"]
    mm = np.memmap(a.trans, dtype="<f8", mode="r", shape=(n // 40, 5))
    out = np.memmap(a.out, dtype="<f8", mode="w+", shape=mm.shape)
    step = 2_000_000
    for i in range(0, len(mm), step):
        blk = np.array(mm[i:i + step])
        blk[:, 0] += d_rng
        blk[:, 1] += d_az
        out[i:i + step] = blk
    out.flush()
    print("[dem_coreg] trans.dat 距离向 %+.2f、方位向 %+.2f 像元" % (d_rng, d_az))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("estimate")
    e.add_argument("--prm", required=True)
    e.add_argument("--slc", required=True)
    e.add_argument("--dem-ra", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--report", required=True)
    e.add_argument("--min-peak", type=float, default=0.2)
    e.add_argument("--min-tiles", type=int, default=6)
    e.add_argument("--max-mad", type=float, default=1.0)
    e.add_argument("--min-relief", type=float, default=40.0)
    e.add_argument("--min-shift", type=float, default=0.3)
    e.add_argument("--max-shift", type=float, default=15.0)
    t = sub.add_parser("apply-trans")
    t.add_argument("--trans", required=True)
    t.add_argument("--report", required=True)
    t.add_argument("--out", required=True)
    a = p.parse_args(argv)
    return cmd_estimate(a) if a.cmd == "estimate" else cmd_apply_trans(a)


if __name__ == "__main__":
    sys.exit(main())
