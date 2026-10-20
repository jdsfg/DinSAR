#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""电建一号 (Dianjian-01) 新版 SLC 产品的元数据解析与干涉几何计算。

产品: <name>.XML (根节点 Metadata, 扁平字段) + <name>.h5 (单精度复数 SLC) + <name>.tiff (8 位幅度快视图)。
本模块只依赖 python3 标准库和 numpy, 供 make_slc_dj1h5.py (转 GMTSAR SLC/PRM/LED) 与
dj1_pair_select.py (像对筛选) 共用。

约定: 时间一律用 "自 2000-01-01T00:00:00 UTC 起的秒" (float64, 微秒精度足够);
      轨道位置/速度为 ECEF (米, 米/秒), 与 GMTSAR 的 LED 文件一致。
"""

import bisect
import datetime as dt
import math
import xml.etree.ElementTree as ET

import numpy as np

C_LIGHT = 299792458.0
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
_EPOCH = dt.datetime(2000, 1, 1)


def parse_time(s):
    s = s.strip()
    fmt = "%Y-%m-%dT%H:%M:%S.%f" if "." in s else "%Y-%m-%dT%H:%M:%S"
    return dt.datetime.strptime(s, fmt)


def to_sec(t):
    return (t - _EPOCH).total_seconds()


def from_sec(sec):
    return _EPOCH + dt.timedelta(seconds=float(sec))


def llh_to_ecef(lat_deg, lon_deg, h=0.0):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    return np.array([(n + h) * math.cos(lat) * math.cos(lon),
                     (n + h) * math.cos(lat) * math.sin(lon),
                     (n * (1.0 - WGS84_E2) + h) * math.sin(lat)])


def _text(root, tag, cast=str, default=None):
    e = root.find(tag)
    if e is None or e.text is None or not e.text.strip():
        return default
    try:
        return cast(e.text.strip())
    except ValueError:
        return default


def _corner(root, tag):
    """'<sample> <line> <lat> <lon>' -> (sample, line, lat, lon)"""
    v = root.find(tag).text.split()
    return int(float(v[0])), int(float(v[1])), float(v[2]), float(v[3])


# ---------------------------------------------------------------------- 轨道: 时间轴规整与动力学外推
GM_EARTH = 3.986004418e14
J2_EARTH = 1.08262668e-3
OMEGA_EARTH = 7.2921150e-5


def regularize_orbit_times(t, p, v, exact_tol=2e-6, ratio_spread_tol=5e-4, step_agree_tol=5e-3, fit_rms_tol=0.05):
    """还原轨道矢量的等间隔时间轴。

    实测的电建一号实时轨道: 161 个矢量的位置对序号是严格的光滑多项式 (拟合残差 < 0.1 mm), 速度与位置的导数严格自洽
    (位置对序号的导数 / 速度 = 采样间隔, 全窗口内离散 1e-9 量级), 说明采样本来是等间隔的; 而 XML 里的时间标签
    被量化过, 抖动 0.1~0.4 ms (沿轨 1~3 m), 直接用会引入沿轨位置误差。
    做法: 采样间隔取速度反推值, 起点取 (标签 - 序号 x 间隔) 的平均。

    只在证据充分时才改标签: 标签本来就等间隔 (精密/快速轨道) 时原样返回; 位置不是光滑多项式、速度反推的间隔
    离散过大、或与标签给出的间隔相差超过 0.5% 时也原样返回。返回 (t_new, info)。"""
    t = np.asarray(t, dtype=float)
    n = len(t)
    info = dict(n=n, regularized=False, step=None, max_tag_residual_s=0.0, reason="")
    if n < 8:
        info["reason"] = "矢量太少"
        return t, info
    i = np.arange(n, dtype=float)
    span_step = (t[-1] - t[0]) / (n - 1)
    info["step"] = float(span_step)
    if np.max(np.abs(t - (t[0] + i * span_step))) <= exact_tol:
        info["reason"] = "时间标签本身等间隔"
        return t, info
    ic = i - i.mean()
    dpdi = np.empty_like(p)
    fit_rms = 0.0
    for k in range(3):
        c = np.polyfit(ic, p[:, k], 5)
        dpdi[:, k] = np.polyval(np.polyder(c), ic)
        fit_rms = max(fit_rms, float(np.sqrt(np.mean((np.polyval(c, ic) - p[:, k]) ** 2))))
    if fit_rms > fit_rms_tol:
        info["reason"] = "位置不是光滑多项式 (拟合 rms %.3f m), 保留原时间标签" % fit_rms
        return t, info
    ratio = np.linalg.norm(dpdi, axis=1) / np.linalg.norm(v, axis=1)
    step = float(np.sum(dpdi * v) / np.sum(v * v))
    if np.ptp(ratio) / step > ratio_spread_tol or abs(step - span_step) / span_step > step_agree_tol:
        info["reason"] = "位置与速度不自洽, 保留原时间标签"
        return t, info
    t_new = float(np.mean(t - i * step)) + i * step
    info.update(regularized=True, step=step, max_tag_residual_s=float(np.max(np.abs(t - t_new))),
                reason="按速度反推的间隔重建等间隔时间轴")
    return t_new, info


def _accel(p, v):
    """地固系 (旋转系) 加速度: 地球引力 (含 J2) + 科里奥利力 + 离心力。"""
    r = math.sqrt(float(p @ p))
    k = 1.5 * J2_EARTH * (WGS84_A / r) ** 2
    z2 = (p[2] / r) ** 2
    f = -GM_EARTH / r ** 3
    return np.array([f * p[0] * (1 + k * (1 - 5 * z2)) + 2 * OMEGA_EARTH * v[1] + OMEGA_EARTH ** 2 * p[0],
                     f * p[1] * (1 + k * (1 - 5 * z2)) - 2 * OMEGA_EARTH * v[0] + OMEGA_EARTH ** 2 * p[1],
                     f * p[2] * (1 + k * (3 - 5 * z2))])


def propagate_state(p, v, dt, step=0.05):
    """把地固系状态 (位置, 速度) 向前 (dt>0) 或向后 (dt<0) 传播 |dt| 秒, 四阶 Runge-Kutta。

    对照有 204 s 长弧的快速精密轨道: 以中心状态外推, ±12 s 内位置误差约 0.1 m, ±50 s 约 1 m,
    远小于实时轨道本身的米级不确定度。"""
    p = np.asarray(p, dtype=float)
    v = np.asarray(v, dtype=float)
    if dt == 0:
        return p, v
    n = max(1, int(math.ceil(abs(dt) / step)))
    h = dt / n
    for _ in range(n):
        k1p, k1v = v, _accel(p, v)
        k2p, k2v = v + 0.5 * h * k1v, _accel(p + 0.5 * h * k1p, v + 0.5 * h * k1v)
        k3p, k3v = v + 0.5 * h * k2v, _accel(p + 0.5 * h * k2p, v + 0.5 * h * k2v)
        k4p, k4v = v + h * k3v, _accel(p + h * k3p, v + h * k3v)
        p = p + h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        v = v + h / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
    return p, v


class Product:
    """一景 SLC 的元数据 (只含干涉/转换需要的字段)。"""

    def __init__(self, xml_path):
        self.xml_path = xml_path
        root = ET.parse(xml_path).getroot()
        if root.tag != "Metadata" or root.find("Orbit_State_Vectors") is None:
            raise ValueError("不是电建一号新版产品的 XML (根节点应为 Metadata): %s" % xml_path)
        self.name = _text(root, "product_name")
        self.h5_name = _text(root, "product_file")
        self.satellite = _text(root, "satellite_name", default="?")
        self.mode = _text(root, "acquisition_mode", default="?")
        self.pol = _text(root, "polarization", default="?")
        self.look_side = (_text(root, "look_side", default="right") or "right").lower()
        self.look_angle = _text(root, "satellite_look_angle", float)
        self.orbit_dir = (_text(root, "orbit_direction", default="?") or "?").upper()
        self.orbit_abs = _text(root, "orbit_absolute_number", int)
        self.orbit_rel = _text(root, "orbit_relative_number", int)
        self.orbit_level = _text(root, "orbit_processing_level", default="?")
        self.acq_id = _text(root, "acquisition_id", int)
        self.fc = _text(root, "carrier_frequency", float)
        self.fs = _text(root, "range_sampling_rate", float)
        self.bw = _text(root, "chirp_bandwidth", float)
        self.chirp_dur = _text(root, "chirp_duration", float)
        self.prf = _text(root, "processing_prf", float) or _text(root, "acquisition_prf", float)
        self.ati = _text(root, "azimuth_time_interval", float)
        self.bw_az = _text(root, "total_processed_bandwidth_azimuth", float)
        self.n_az = _text(root, "number_of_azimuth_samples", int)
        self.n_rg = _text(root, "number_of_range_samples", int)
        self.r0 = _text(root, "slant_range_to_first_pixel", float)
        self.dr = _text(root, "slant_range_spacing", float)
        self.daz_ground = _text(root, "azimuth_ground_spacing", float)
        self.t_acq_start = to_sec(parse_time(_text(root, "acquisition_start_utc")))
        self.t_zd_start = to_sec(parse_time(_text(root, "zerodoppler_start_utc")))
        self.t_zd_end = to_sec(parse_time(_text(root, "zerodoppler_end_utc")))
        self.heading = _text(root, "heading", float)
        self.inc_near = _text(root, "incidence_near", float)
        self.inc_center = _text(root, "incidence_center", float)
        self.inc_far = _text(root, "incidence_far", float)
        self.avg_height = _text(root, "avg_scene_height", float, 0.0)
        self.calibration = _text(root, "calibration_factor", float)
        self.corners = {k: _corner(root, "coord_" + k)
                        for k in ("first_near", "first_far", "last_far", "last_near", "center")}
        self.dc_poly = self._read_poly(root.find("Doppler_Centroid_Coefficients/dc_coefficients_list"))
        self.dc_ref_time = _text(root, "Doppler_Centroid_Coefficients/dc_coefficients_list/reference_pixel_time", float)
        sv = root.find("Orbit_State_Vectors").findall("orbit_vector")
        self.sv_t = np.array([to_sec(parse_time(v.find("time").text)) for v in sv])
        self.sv_p = np.array([[float(v.find(k).text) for k in ("posX", "posY", "posZ")] for v in sv])
        self.sv_v = np.array([[float(v.find(k).text) for k in ("velX", "velY", "velZ")] for v in sv])
        order = np.argsort(self.sv_t)
        self.sv_t, self.sv_p, self.sv_v = self.sv_t[order], self.sv_p[order], self.sv_v[order]
        self.sv_t_tags = self.sv_t.copy()
        self.sv_t, self.orbit_info = regularize_orbit_times(self.sv_t, self.sv_p, self.sv_v)

    @staticmethod
    def _read_poly(node):
        if node is None:
            return [0.0]
        out = []
        for k in sorted(node.findall("coefficient"), key=lambda e: int(e.find("number").text)):
            out.append(float(k.find("value").text))
        return out or [0.0]

    # ------------------------------------------------------------------ 基本量
    @property
    def date(self):
        return from_sec(self.t_zd_start).strftime("%Y%m%d")

    @property
    def t_mid(self):
        return 0.5 * (self.t_zd_start + self.t_zd_end)

    @property
    def wavelength(self):
        return C_LIGHT / self.fc

    @property
    def center_llh(self):
        _, _, lat, lon = self.corners["center"]
        return lat, lon, self.avg_height

    def footprint(self):
        """足迹四边形 [(lat, lon)] (近-远-远-近顺序)"""
        return [self.corners[k][2:4] for k in ("first_near", "first_far", "last_far", "last_near")]

    # ------------------------------------------------------------------ 轨道
    def sat_state(self, t):
        """三次 Hermite 插值的卫星位置和速度 (ECEF)。t 超出轨道时间范围会抛出 ValueError。"""
        ts = self.sv_t
        # 实时轨道的首尾矢量与成像起止时刻只差 0.1 ms 量级, 允许 2 ms 内的三次外推
        if t < ts[0] - 2e-3 or t > ts[-1] + 2e-3:
            raise ValueError("时刻超出轨道矢量范围")
        i = min(max(bisect.bisect_right(ts, t) - 1, 0), len(ts) - 2)
        h = ts[i + 1] - ts[i]
        s = (t - ts[i]) / h
        p0, p1, v0, v1 = self.sv_p[i], self.sv_p[i + 1], self.sv_v[i], self.sv_v[i + 1]
        h00, h10 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s
        h01, h11 = -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
        pos = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1
        d00, d10 = 6 * s ** 2 - 6 * s, 3 * s ** 2 - 4 * s + 1
        d01, d11 = -6 * s ** 2 + 6 * s, 3 * s ** 2 - 2 * s
        vel = (d00 * p0 + d10 * h * v0 + d01 * p1 + d11 * h * v1) / h
        return pos, vel

    def xml_velocity(self, t, half_window=0.5, min_pts=5):
        """XML 速度矢量在 t 时刻的值 (ECEF): 取 t 前后 half_window 秒内的速度矢量 (不足 min_pts 个就取最近的 min_pts 个),
        各分量做二次拟合。

        与 sat_state() 返回的速度不同: 后者是位置 Hermite 样条的导数。位置序列与速度序列自洽时两者一致 (在矢量节点上
        总是一致); 不自洽时, 样条导数在两个矢量之间偏离速度矢量 (中点处偏差为 1.5 x (位置有限差分 - 速度矢量))。
        实时轨道常见这种情形: 速度矢量很光滑, 位置只取整到亚毫米, 但 (位置 - 速度的积分) 随时间漂移; 以样条导数作外推起点,
        中心速度就会带上这份偏差, 外推轨道在景内累积成随方位向线性增长的视线向误差, 干涉图上表现为方位向的平面斜坡。
        所以外推起点的速度取速度矢量 (对应的不自洽程度用 pv_mismatch 监控)。"""
        ts = self.sv_t
        idx = np.nonzero(np.abs(ts - t) <= half_window)[0]
        if len(idx) < min_pts:
            idx = np.sort(np.argsort(np.abs(ts - t))[:min_pts])
        x = ts[idx] - t
        deg = min(2, len(idx) - 1)
        return np.array([np.polyfit(x, self.sv_v[idx, k], deg)[-1] for k in range(3)])

    def gmtsar_orbit(self, margin=10.0, step=0.5, min_margin=6.0):
        """生成 GMTSAR 能直接使用的等间隔轨道, 覆盖成像时段前后至少 margin 秒。

        GMTSAR 不读 LED 每行的时间列, 而是用"起始时间 + 序号 x dsec"取时间, 所以必须严格等间隔; calc_dop_orb
        还要在影像起止时刻前后各 2 s 以上取轨道 (差分求速度), 地理编码迭代也可能越过影像边界。
        电建一号实时轨道只覆盖成像时段本身 (SP 景首尾余量为 0), 直接用会在 calc_dop_orb 报
        "hermite interpolation point outside of data constraints"。

        原始轨道 (时间轴规整后) 等间隔且两端余量都不小于 min_margin 时原样使用; 否则以场景中心时刻的状态为起点,
        用 propagate_state 向两侧外推, 重建 step 间隔的等间隔轨道。
        外推起点: 位置取 sat_state(t_c) (三次 Hermite, 矢量间隔很小时位置误差可忽略), 速度取 xml_velocity(t_c)
        (XML 速度矢量的局部拟合), 而不是位置样条的导数, 见 xml_velocity 的说明。
        返回 dict(t0, dt, t, p, v, extended, max_dev_m, max_dev_window_m, pv_mismatch_mps):
        max_dev_m 是外推轨道与全部原始位置矢量 (含成像时段之外的) 的最大位置差; max_dev_window_m 只统计成像时段内的矢量;
        pv_mismatch_mps 是成像时段附近 位置给出的平均速度 与 速度矢量给出的平均速度 之差 (m/s, 见 pv_mismatch):
        位置与速度自洽时接近 0; 乘以成像时长约等于该景内位置与速度两种说法的相对漂移, 超过几厘米就说明轨道的
        位置与速度不自洽 (已按速度矢量定中心速度, 但这一景的轨道精度需要留意)。"""
        t_lo, t_hi = self.t_zd_start - margin, self.t_zd_end + margin
        ts = self.sv_t
        n0 = len(ts)
        dts = np.diff(ts)
        uniform = n0 >= 10 and np.ptp(dts) <= 2e-6
        if uniform and ts[0] <= self.t_zd_start - min_margin and ts[-1] >= self.t_zd_end + min_margin:
            return dict(t0=float(ts[0]), dt=float((ts[-1] - ts[0]) / (n0 - 1)), t=ts.copy(), p=self.sv_p.copy(),
                        v=self.sv_v.copy(), extended=False, max_dev_m=0.0, max_dev_window_m=0.0,
                        pv_mismatch_mps=self.pv_mismatch())
        t_c = self.t_mid
        p_c = self.sat_state(t_c)[0]
        v_c = self.xml_velocity(t_c)
        pv_mismatch = self.pv_mismatch()
        n = int(math.ceil((t_hi - t_lo) / step)) + 1
        t = t_lo + step * np.arange(n)
        P, V = np.empty((n, 3)), np.empty((n, 3))
        k_c = int(np.searchsorted(t, t_c))
        for ks in (range(k_c, n), range(k_c - 1, -1, -1)):
            pp, vv, tt = p_c, v_c, t_c
            for k in ks:
                pp, vv = propagate_state(pp, vv, t[k] - tt)
                tt = t[k]
                P[k], V[k] = pp, vv
        dev = dev_win = 0.0
        for tk, pk in zip(ts, self.sv_p):
            d = float(np.linalg.norm(propagate_state(p_c, v_c, tk - t_c)[0] - pk))
            dev = max(dev, d)
            if self.t_zd_start <= tk <= self.t_zd_end:
                dev_win = max(dev_win, d)
        return dict(t0=float(t[0]), dt=float(step), t=t, p=P, v=V, extended=True, max_dev_m=dev,
                    max_dev_window_m=dev_win, pv_mismatch_mps=pv_mismatch)

    def pv_mismatch(self, pad=1.0, min_pts=5):
        """成像时段前后各 pad 秒内, 位置给出的平均速度与速度矢量给出的平均速度之差 (m/s):
        |(位置净位移 - 速度矢量的梯形积分) / 时间跨度|。位置与速度自洽时接近 0; 不依赖场景中心落在矢量的哪个位置,
        也不受位置取整噪声影响 (净位移只用首尾两个位置)。窗口内矢量不足 min_pts 个时取场景中心附近最近的 min_pts 个。"""
        ts = self.sv_t
        idx = np.nonzero((ts >= self.t_zd_start - pad) & (ts <= self.t_zd_end + pad))[0]
        if len(idx) < min_pts:
            idx = np.sort(np.argsort(np.abs(ts - self.t_mid))[:min_pts])
        if len(idx) < 2:
            return float("nan")
        dt = np.diff(ts[idx])
        integ = np.sum(0.5 * (self.sv_v[idx[1:]] + self.sv_v[idx[:-1]]) * dt[:, None], axis=0)
        return float(np.linalg.norm(self.sv_p[idx[-1]] - self.sv_p[idx[0]] - integ) / (ts[idx[-1]] - ts[idx[0]]))

    def zero_doppler_time(self, target, t_guess=None):
        """目标点 (ECEF) 的零多普勒时刻: 使 (P(t)-T)·V(t)=0 的 t。二分法, 不外推。"""
        lo, hi = self.sv_t[0], self.sv_t[-1]

        def f(t):
            p, v = self.sat_state(t)
            return float(np.dot(p - target, v))

        flo, fhi = f(lo), f(hi)
        if flo * fhi > 0:
            raise ValueError("目标点的零多普勒时刻不在轨道矢量覆盖范围内")
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if flo * fm <= 0:
                hi, fhi = mid, fm
            else:
                lo, flo = mid, fm
        return 0.5 * (lo + hi)

    def target_geometry(self, target):
        """目标点在本景中的 (零多普勒时刻, 卫星位置, 斜距, 方位向行号(从0起), 距离向列号(从0起))"""
        t = self.zero_doppler_time(target)
        p, _ = self.sat_state(t)
        r = float(np.linalg.norm(p - target))
        line = (t - self.t_zd_start) / self.ati
        col = (r - self.r0) / self.dr
        return t, p, r, line, col


# ---------------------------------------------------------------------- 像对几何
def _local_xy(latlon, lat0, lon0):
    k = math.cos(math.radians(lat0))
    return [((lon - lon0) * k * 111.32, (lat - lat0) * 110.57) for lat, lon in latlon]


def _clip(subject, clipper):
    """Sutherland-Hodgman: 凸多边形裁剪, 点序须同为逆时针。"""
    def inside(p, a, b):
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= 0

    def cross(p1, p2, a, b):
        x1, y1, x2, y2 = p1[0], p1[1], p2[0], p2[1]
        x3, y3, x4, y4 = a[0], a[1], b[0], b[1]
        den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if abs(den) < 1e-12:
            return p2
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))

    out = subject
    for i in range(len(clipper)):
        a, b = clipper[i], clipper[(i + 1) % len(clipper)]
        inp, out = out, []
        if not inp:
            break
        s = inp[-1]
        for e in inp:
            if inside(e, a, b):
                if not inside(s, a, b):
                    out.append(cross(s, e, a, b))
                out.append(e)
            elif inside(s, a, b):
                out.append(cross(s, e, a, b))
            s = e
    return out


def _area(poly):
    return 0.5 * sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1]
                     for i in range(len(poly)))


def _ccw(poly):
    return poly if _area(poly) > 0 else poly[::-1]


def footprint_overlap(a, b):
    """返回 (重叠面积 km², 占 a 的比例, 占 b 的比例)。足迹都是凸四边形。"""
    lat0, lon0 = a.center_llh[0], a.center_llh[1]
    pa = _ccw(_local_xy(a.footprint(), lat0, lon0))
    pb = _ccw(_local_xy(b.footprint(), lat0, lon0))
    inter = _clip(pa, pb)
    ov = abs(_area(inter)) if len(inter) >= 3 else 0.0
    return ov, ov / abs(_area(pa)), ov / abs(_area(pb))


def pair_geometry(m, s):
    """以 m 为主影像、s 为从影像, 在主影像场景中心点处计算干涉几何。

    返回 dict: 时间基线(天)、垂直/平行基线(米)、临界基线(米)、预测的从影像相对主影像的偏移
    (距离向/方位向像素, 按轨道和 XML 起始时间推算, 受实时轨道误差影响, 只作参考)、
    视角/入射角/航向差、多普勒中心差。

    与 GMTSAR SAT_baseline 对照时注意两点 (基线矢量本身与它一致, alpha 角和基线长度相差 1 mm 以内):
      * b_par > 0 表示从影像更远; GMTSAR 的 B_parallel 符号相反 (正值 = 从影像更近)。
      * GMTSAR 把"入射角"当作视角代入 B∥=B·sin(θ-α) (比真实视角大地球曲率角, 约 3~4°), 且取场景起始行的基线,
        所以它的 B⊥ 与这里的相差约 B∥×0.06 (B∥=33 m 时差 2 m), B∥ 相差约 B⊥×0.06。这里用卫星处的真实视角。

    预测偏移与两对真实像对的多尺度配准实测值对比 (主影像同一像素处, 预测减实测): 距离向 SM -1.0~-1.6 像素、SP -4 像素,
    都约 -1 m; 方位向 SM -0.1~+3 行、SP -9~-14 行 (3~5 m)。原因是实时轨道:
      * 零多普勒时刻对速度误差极敏感 (力臂 650~700 km): 视线向 1 cm/s 约折合 SM 0.3~0.5 行、SP 3~13 行 (行距 0.19 m);
      * XML 实时轨道的速度与位置导数有 1~20 cm/s 的平滑系统差 (相邻矢量的差高度自相关, 不是白噪声), 对它做 Hermite
        插值会产生 70~100 m/s² 的假加速度 (真实约 8 m/s²), 零多普勒函数的斜率在相邻采样间隔里相差近 10 倍。改用转换器
        写进 LED 的动力学自洽轨道 (gmtsar_orbit) 预测时, 方位向差变成近似常数 (SM +2.7 行、SP -15.3 行, 沿场景变化
        < 0.15 行), 对应约 5 m / 3 m 的沿轨绝对位置误差。
    所以只用来排序和估计搜索范围, 不能当配准初值。方位向拉伸 (配准模型的行间斜率) 应取 ATI 之比, 与实测吻合到 1% 以内。"""
    lat, lon, h = m.center_llh
    target = llh_to_ecef(lat, lon, h)
    tm, pm, rm, line_m, col_m = m.target_geometry(target)
    _, ps, _, line_s, col_s = s.target_geometry(target)
    b = ps - pm
    u = (pm - target) / rm
    vm = m.sat_state(tm)[1]
    vhat = vm / np.linalg.norm(vm)
    # 沿轨分量不属于干涉基线 (方位向配准会吸收它), 要先去掉再分解。两条轨道的速度方向只差 0.001° 时,
    # 在 700 km 斜距上就会出现约 8 m 的沿轨分量, 若不去掉会把 B⊥ 高估数米 (与 GMTSAR SAT_baseline 对照得到)
    b_along = float(np.dot(b, vhat))
    b_t = b - b_along * vhat
    b_par = float(np.dot(b_t, u))
    b_perp = math.sqrt(max(float(np.dot(b_t, b_t)) - b_par ** 2, 0.0))
    # 符号: 取从影像相对主影像视线的"向上/向下"方向 (沿 u x v 方向为正), 便于区分两侧基线
    side = np.cross(u, vhat)
    sign = 1.0 if float(np.dot(b_t, side)) >= 0 else -1.0
    inc = math.radians(m.inc_center)
    res_r = C_LIGHT / (2.0 * m.bw) if m.bw else m.dr
    b_crit = m.wavelength * rm * math.tan(inc) * math.sin(inc) / (2.0 * res_r)
    f_dc_m = m.dc_poly[0]
    f_dc_s = s.dc_poly[0]
    return dict(
        dt_days=(s.t_mid - m.t_mid) / 86400.0,
        b_perp=sign * b_perp, b_par=b_par, b_along=b_along, b_crit=b_crit,
        range_offset_px=col_s - col_m, azimuth_offset_px=line_s - line_m,
        d_look=(s.look_angle or 0) - (m.look_angle or 0),
        d_inc=s.inc_center - m.inc_center, d_heading=((s.heading - m.heading + 180) % 360) - 180,
        d_fdc=f_dc_s - f_dc_m, dc_decorr=abs(f_dc_s - f_dc_m) / (m.bw_az or m.prf),
    )
