#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""电建一号新版产品 (XML + h5 单精度复数 SLC) -> GMTSAR 的 SLC / PRM / LED。

用法:
    make_slc_dj1h5.py PRODUCT.XML PRODUCT.h5 NAME [--scale S] [--info]
    make_slc_dj1h5.py PRODUCT.XML - NAME --metadata-only      # SLC 已存在, 只重写 PRM/LED/转换信息

输出: NAME.SLC (int16 I/Q 交错, 行=方位向, 列=距离向) | NAME.PRM | NAME.LED | NAME.convert.txt (转换信息)
与旧版 make_slc_dj1 (类 Sentinel-1 的 XML + CInt16 TIFF) 的输出格式一致, 可直接进入后续流水线。
依赖 numpy 和 h5py; 系统 Python 里没有 h5py 时, 会在 $HOME/gmtsar/_tools/pylibs 和脚本旁的 pylibs 里找。

说明:
  * 新版 XML 的结构与旧版完全不同 (根节点 Metadata, 扁平字段), 解析见 dj1_meta.py。
  * GMTSAR 的 SLC 是 16 位整数, h5 是单精度复数, 需要一个整体比例因子。默认按幅度中位数映射到约 150,
    同时保证 99.99% 分位不超过 32000 (只有极少数最强散射点会被截断)。比例因子和截断样点数写入 NAME.convert.txt。
  * 与旧版一致, 行数和列数取 4 的整数倍 (丢弃多余的末尾行/列)。
  * LED 头里的时间间隔 dt 用轨道矢量实际间隔, 而不是旧版按 0.01 秒截断的值, 因为 GMTSAR 按 dt 等间隔插值轨道。
  * PRM 的 fd1 固定写 0 (理由见 write_prm)。PRM 里没有 SC_vel / earth_radius / SC_height*: 按 GMTSAR 惯例由 calc_dop_orb
    追加, 且从影像要沿用主影像的 earth_radius, 所以这一步属于"成对"处理, 不放在单景转换里:
        calc_dop_orb M.PRM M.log 0 0;  calc_dop_orb S.PRM S.log <M 的 earth_radius> 0;  cat M.log >> M.PRM; cat S.log >> S.PRM
"""

import argparse
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dj1_meta as dm  # noqa: E402

TARGET_MEDIAN = 150.0
CLIP_LEVEL = 32000.0
MIN_MEDIAN = 20.0


def _import_h5py():
    try:
        import h5py
        return h5py
    except ImportError:
        here = os.path.dirname(os.path.abspath(__file__))
        for cand in (os.path.join(os.environ.get("HOME", ""), "gmtsar", "_tools", "pylibs"),
                     os.path.join(here, "pylibs")):
            if os.path.isdir(cand):
                sys.path.insert(0, cand)
        import h5py
        return h5py


class SlcReader:
    """把 h5 里的复数数据统一成逻辑数组 S[line, sample] (行=方位向, 列=距离向)。"""

    def __init__(self, h5_path, prod):
        h5py = _import_h5py()
        self.h5 = h5py.File(h5_path, "r")
        self.prod = prod
        self.parts = None
        self.ds = None
        naz, nrg = prod.n_az, prod.n_rg
        cands = []

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and obj.ndim == 2 and sorted(obj.shape) == sorted((naz, nrg)):
                cands.append((name, obj))
        self.h5.visititems(visit)
        if not cands:
            raise ValueError("h5 里没有 %d x %d 的二维数据集: %s" % (naz, nrg, self.describe()))
        complex_ds = [c for c in cands if c[1].dtype.kind == "c" or (c[1].dtype.names and len(c[1].dtype.names) == 2)]
        if complex_ds:
            self.name, self.ds = complex_ds[0]
        else:
            real = {n.split("/")[-1].lower(): (n, d) for n, d in cands}
            pick = None
            for a, b in (("s_i", "s_q"), ("i", "q"), ("real", "imag"), ("re", "im")):
                if a in real and b in real:
                    pick = (real[a], real[b])
                    break
            if pick is None and len(cands) == 2:
                pick = (cands[0], cands[1])
            if pick is None:
                raise ValueError("无法识别复数数据集: %s" % self.describe())
            self.parts = (pick[0][1], pick[1][1])
            self.name = "%s + %s" % (pick[0][0], pick[1][0])
        first = self.ds if self.ds is not None else self.parts[0]
        self.stored_shape = first.shape
        self.transposed = first.shape == (nrg, naz) and naz != nrg
        if first.shape == (naz, nrg):
            self.transposed = False
        self.chunks = first.chunks

    def describe(self):
        lines = []

        def visit(name, obj):
            if hasattr(obj, "shape"):
                lines.append("%s shape=%s dtype=%s chunks=%s" % (name, obj.shape, obj.dtype, obj.chunks))
        self.h5.visititems(visit)
        return "; ".join(lines)

    def _read(self, l0, l1, c0, c1):
        """逻辑坐标 [l0:l1, c0:c1] -> complex64 (l1-l0, c1-c0)"""
        def get(d):
            a = d[c0:c1, l0:l1].T if self.transposed else d[l0:l1, c0:c1]
            return a
        if self.ds is not None:
            a = get(self.ds)
            if a.dtype.names:
                f = a.dtype.names
                return (a[f[0]].astype(np.float32) + 1j * a[f[1]].astype(np.float32)).astype(np.complex64)
            return np.ascontiguousarray(a, dtype=np.complex64)
        re_, im_ = get(self.parts[0]), get(self.parts[1])
        return (re_.astype(np.float32) + 1j * im_.astype(np.float32)).astype(np.complex64)

    def lines(self, l0, l1, ncols):
        return self._read(l0, l1, 0, ncols)

    def close(self):
        self.h5.close()


def choose_scale(reader, nlines, ncols, n_blocks=48, block_lines=64, seed=1):
    rng = np.random.default_rng(seed)
    starts = np.sort(rng.integers(0, max(nlines - block_lines, 1), size=n_blocks))
    amps = []
    for s in starts:
        blk = reader.lines(int(s), int(s) + block_lines, ncols)
        a = np.abs(blk).ravel()
        amps.append(a[::7][a[::7] > 0])
    a = np.concatenate(amps)
    med = float(np.median(a))
    p9999 = float(np.percentile(a, 99.99))
    scale = TARGET_MEDIAN / med
    if scale * p9999 > CLIP_LEVEL:
        scale = CLIP_LEVEL / p9999
    if scale * med < MIN_MEDIAN:
        scale = MIN_MEDIAN / med
    return scale, med, p9999


def ydsod(sec):
    """自 2000-01-01 起的秒 -> (年, 年内日, 当天秒), 当天秒保留双精度, 不经过 datetime 的微秒截断"""
    day0 = math.floor(sec / 86400.0) * 86400.0
    d = dm.from_sec(day0)
    return d.year, d.timetuple().tm_yday, sec - day0


def write_led(prod, path):
    """写 GMTSAR 的 LED 轨道文件, 返回 prod.gmtsar_orbit() 的结果 (供日志记录)。

    GMTSAR 只用表头的 起始时间 + 序号 x dsec 取时间, 所以 dsec 要给足有效数字
    (160 个间隔累积, 6 位小数就会带来米级的沿轨误差); 轨道的覆盖范围与外推见 Product.gmtsar_orbit。"""
    orb = prod.gmtsar_orbit()
    t, pos, vel = orb["t"], orb["p"], orb["v"]
    y0, j0, s0 = ydsod(t[0])
    with open(path, "w") as f:
        f.write("%d %d %d %.9f %.9f \n" % (len(t), y0, j0, s0, orb["dt"]))
        for i in range(len(t)):
            y, j, s = ydsod(t[i])
            f.write("%d %d %.9f %.6f %.6f %.6f %.9f %.9f %.9f \n" % (
                y, j, s, pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2]))
    return orb


def write_prm(prod, name, nrg, naz, path):
    # GMTSAR 按"PRM 与 SLC/LED 在同一目录"的习惯读文件, PRM 里只能写文件名, 不能带目录
    name = os.path.basename(name)
    d0 = dm.from_sec(prod.t_zd_start)
    doy = d0.timetuple().tm_yday
    sod = d0.hour * 3600 + d0.minute * 60 + d0.second + d0.microsecond * 1e-6
    clock_start = doy + sod / 86400.0
    prf = 1.0 / prod.ati
    clock_stop = clock_start + naz / prf / 86400.0
    sc_start = d0.year * 1000 + clock_start
    sc_stop = d0.year * 1000 + clock_stop
    dur = prod.chirp_dur or 2.68e-5
    slope = (prod.bw / dur) if prod.bw else 0.0
    # 产品的方位向时间是零多普勒时间 (XML zerodoppler_start_utc), 而 GMTSAR 的 SAT_llt2rat 把非零 fd1 当作"方位向以多普勒
    # 中心为基准"并据此平移方位坐标 (-lambda*fd1/(2*rdd) 秒): 写入 XML 的多普勒常数项 (SP 景 -151 Hz) 会让 DEM 地形相位
    # 和地理编码在方位向错开约 1200 行。GMTSAR 对 Sentinel-1 这类零多普勒产品也是强制 fd1=0 (calc_dop_orb ... 0)。
    # XML 里的原值记在 NAME.convert.txt 的 xml_doppler_fd1。
    fd1 = 0.0
    lookdir = "R" if prod.look_side.startswith("r") else "L"
    orbdir = "A" if prod.orbit_dir.startswith("A") else "D"
    rows = [
        ("num_valid_az", "%d" % naz), ("nrows", "%d" % naz), ("first_line", "1"), ("deskew", "n"),
        ("caltone", "0.000000"), ("st_rng_bin", "1"), ("Flip_iq", "n"), ("offset_video", "n"),
        ("az_res", "0.000000"), ("nlooks", "1"), ("chirp_ext", "0"), ("scnd_rng_mig", "0"),
        ("rng_spec_wgt", "1.000000"), ("rm_rng_band", "0.200000"), ("rm_az_band", "0.000000"),
        ("rshift", "0"), ("ashift", "0"), ("stretch_r", "0"), ("stretch_a", "0"), ("a_stretch_r", "0"),
        ("a_stretch_a", "0"), ("first_sample", "1"), ("SC_identity", "14"),
        ("rng_samp_rate", "%.6f" % prod.fs), ("input_file", "%s.raw" % name),
        ("num_rng_bins", "%d" % nrg), ("bytes_per_line", "%d" % (4 * nrg)), ("good_bytes_per_line", "%d" % (4 * nrg)),
        ("PRF", "%.6f" % prf), ("pulse_dur", "%.6e" % dur), ("near_range", "%.6f" % prod.r0),
        ("num_lines", "%d" % naz), ("num_patches", "1"),
        ("SC_clock_start", "%.10f" % sc_start), ("SC_clock_stop", "%.10f" % sc_stop),
        ("clock_start", "%.12f" % clock_start), ("clock_stop", "%.12f" % clock_stop),
        ("led_file", "%s.LED" % name), ("orbdir", orbdir), ("lookdir", lookdir),
        ("radar_wavelength", "%.7f" % prod.wavelength), ("chirp_slope", "%.6g" % slope),
        ("rng_samp_rate", "%.6f" % prod.fs), ("I_mean", "1"), ("Q_mean", "1"),
        ("equatorial_radius", "6378137.000000"), ("polar_radius", "6356752.310000"),
        ("fd1", "%.6f" % fd1), ("fdd1", "0.000000"), ("fddd1", "0.000000"),
        ("sub_int_r", "0.000000"), ("sub_int_a", "0.000000"), ("SLC_file", "%s.SLC" % name),
        ("dtype", "a"), ("SLC_scale", "1.000000"),
    ]
    with open(path, "w") as f:
        for k, v in rows:
            f.write("%s = %s\n" % (k, v))


def convert(xml, h5, name, scale=None, block_lines=256, quiet=False):
    t0 = time.time()
    prod = dm.Product(xml)
    reader = SlcReader(h5, prod)
    naz = prod.n_az - prod.n_az % 4
    nrg = prod.n_rg - prod.n_rg % 4
    if not quiet:
        print("[make_slc_dj1h5] %s  %s  %d 行 x %d 列 (取 4 的倍数: %d x %d)  数据集 %s %s 转置=%s 分块=%s" % (
            prod.name, prod.mode, prod.n_az, prod.n_rg, naz, nrg, reader.name, reader.stored_shape,
            reader.transposed, reader.chunks), flush=True)
    if scale is None:
        scale, med, p9999 = choose_scale(reader, naz, nrg)
        if not quiet:
            print("[make_slc_dj1h5] 幅度中位数 %.4g, 99.99%% 分位 %.4g -> 比例因子 %.5g (中位数映射到 %.0f)" % (
                med, p9999, scale, med * scale), flush=True)
    slc_path = name + ".SLC"
    tmp = slc_path + ".tmp"
    clipped = 0
    with open(tmp, "wb") as out:
        for l0 in range(0, naz, block_lines):
            l1 = min(l0 + block_lines, naz)
            z = reader.lines(l0, l1, nrg)
            iq = np.empty((l1 - l0, 2 * nrg), np.int16)
            re_ = np.rint(z.real * scale)
            im_ = np.rint(z.imag * scale)
            clipped += int(np.count_nonzero(np.abs(re_) > 32767) + np.count_nonzero(np.abs(im_) > 32767))
            iq[:, 0::2] = np.clip(re_, -32767, 32767)
            iq[:, 1::2] = np.clip(im_, -32767, 32767)
            out.write(iq.tobytes())
            if not quiet and (l0 // block_lines) % 20 == 0:
                print("\r[make_slc_dj1h5] %5.1f%%" % (100.0 * l1 / naz), end="", flush=True)
    os.replace(tmp, slc_path)
    reader.close()
    write_metadata(prod, name, naz, nrg, dict(scale="%.8g" % scale, clipped_samples=clipped,
                                              source=os.path.basename(h5)), quiet)
    if not quiet:
        print("[make_slc_dj1h5] 完成: %s.{SLC,PRM,LED}  截断样点 %d  耗时 %.0f s" % (name, clipped, time.time() - t0))
    return prod


def write_metadata(prod, name, naz, nrg, extra, quiet=False):
    """写 PRM、LED 和转换信息文件 NAME.convert.txt (extra 里的键值先写入, 再写轨道处理信息)。"""
    write_prm(prod, name, nrg, naz, name + ".PRM")
    orb = write_led(prod, name + ".LED")
    oi = prod.orbit_info
    with open(name + ".convert.txt", "w") as f:
        for k, v in extra.items():
            f.write("%s = %s\n" % (k, v))
        f.write("orbit_level = %s\norbit_time_regularized = %s\norbit_max_tag_residual_s = %.3e\n" % (
            prod.orbit_level, oi["regularized"], oi["max_tag_residual_s"]))
        f.write("led_vectors = %d\nled_step_s = %.9f\nled_extended = %s\nled_max_dev_m = %.3f\n" % (
            len(orb["t"]), orb["dt"], orb["extended"], orb["max_dev_m"]))
        f.write("xml_doppler_fd1 = %.6f\n" % (prod.dc_poly[0] if prod.dc_poly else 0.0))
    if not quiet:
        print("\n[make_slc_dj1h5] 轨道 (%s): %s; LED %d 个矢量, 间隔 %.3f s, %s, 与原始矢量最大位置差 %.2f m" % (
            prod.orbit_level, oi["reason"], len(orb["t"]), orb["dt"],
            "已外推到场景前后" if orb["extended"] else "沿用原始矢量", orb["max_dev_m"]))
    return orb


def _read_kv(path):
    out = {}
    if os.path.isfile(path):
        with open(path) as f:
            for ln in f:
                if " = " in ln:
                    k, v = ln.rstrip("\n").split(" = ", 1)
                    out[k.strip()] = v.strip()
    return out


def metadata_only(xml, name, quiet=False):
    """只重写 PRM/LED (SLC 已存在): 轨道处理逻辑更新后不必重读几 GB 的 h5。
    先核对已有 SLC 的大小与 XML 给出的尺寸一致, 避免 PRM 与 SLC 对不上。"""
    prod = dm.Product(xml)
    naz = prod.n_az - prod.n_az % 4
    nrg = prod.n_rg - prod.n_rg % 4
    slc = name + ".SLC"
    if not os.path.isfile(slc) or os.path.getsize(slc) != 4 * naz * nrg:
        raise ValueError("%s 不存在或大小不是 %d x %d x 4 字节" % (slc, naz, nrg))
    old = _read_kv(name + ".convert.txt") or _read_kv(name + ".scale")
    extra = {k: old[k] for k in ("scale", "clipped_samples", "source") if k in old}
    write_metadata(prod, name, naz, nrg, extra, quiet)
    return prod


def main(argv=None):
    ap = argparse.ArgumentParser(description="电建一号新版产品 -> GMTSAR SLC/PRM/LED")
    ap.add_argument("xml")
    ap.add_argument("h5", help="产品 h5 (用 --metadata-only 时不读取, 可写 -)")
    ap.add_argument("name", help="输出前缀, 例如 20260904 -> 20260904.SLC/.PRM/.LED")
    ap.add_argument("--scale", type=float, help="int16 比例因子 (默认按幅度统计自动选择)")
    ap.add_argument("--block-lines", type=int, default=256)
    ap.add_argument("--info", action="store_true", help="只打印 h5 结构和元数据, 不转换")
    ap.add_argument("--metadata-only", action="store_true", help="SLC 已存在, 只重写 PRM/LED/转换信息")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)
    if args.metadata_only:
        metadata_only(args.xml, args.name, args.quiet)
        return 0
    if args.info:
        prod = dm.Product(args.xml)
        r = SlcReader(args.h5, prod)
        print(r.describe())
        print("数据集: %s 存储形状 %s 需转置: %s" % (r.name, r.stored_shape, r.transposed))
        return 0
    convert(args.xml, args.h5, args.name, args.scale, args.block_lines, args.quiet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
