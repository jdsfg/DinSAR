#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""电建一号新版 SLC 产品: 干涉像对自动筛选 (只读 XML 元数据, 不需要解压 h5)。

用法:
    dj1_pair_select.py ROOT [--max-dt 40] [--json pairs.json] [--md pairs.md] [--all]

ROOT 下递归查找 *_SLC_*.XML。对每一对影像 (较早的做主影像):
  1. 硬性条件 (任一不满足即判 X, 不适合常规干涉): 同成像模式、同极化、同升降轨、同侧视、
     雷达参数一致 (载频/带宽/采样率/PRF 偏差不超过 1%)、视角差 <= 0.3°、航向差 <= 1°、
     足迹重叠 >= 30% (取较小一景为基准)。
  2. 几何/时间质量: 垂直基线 B⊥ (在主影像场景中心按轨道计算)、临界基线、多普勒中心差、时间基线。
  3. 先验相干性估计 = (1 - B⊥/B临界) x (1 - |Δf_dc|/B_az) x exp(-Δt/24天) x 重叠比例,
     只用于排序, 不代替配准后实际测得的相干性。
评级: A 推荐 | B 可用 | X 不适合。
"""

import argparse
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dj1_meta as dm  # noqa: E402

TAU_DAYS = 24.0


def load_products(root):
    files = sorted(glob.glob(os.path.join(root, "**", "*_SLC_*.XML"), recursive=True)
                   + glob.glob(os.path.join(root, "**", "*_SLC_*.xml"), recursive=True))
    out = []
    for f in files:
        try:
            out.append(dm.Product(f))
        except Exception as e:  # noqa: BLE001
            print("跳过 %s: %s" % (f, e), file=sys.stderr)
    out.sort(key=lambda p: p.t_mid)
    return out


def _rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-12)


def evaluate_pair(m, s):
    """m 为较早的主影像。返回 dict (含评级、原因、几何量)。"""
    hard, soft = [], []
    if m.mode != s.mode:
        hard.append("成像模式不同 (%s/%s)" % (m.mode, s.mode))
    if m.pol != s.pol:
        hard.append("极化不同 (%s/%s)" % (m.pol, s.pol))
    if m.orbit_dir != s.orbit_dir:
        hard.append("升降轨不同 (%s/%s)" % (m.orbit_dir, s.orbit_dir))
    if m.look_side != s.look_side:
        hard.append("侧视方向不同 (%s/%s)" % (m.look_side, s.look_side))
    for name, a, b in (("载频", m.fc, s.fc), ("距离向带宽", m.bw, s.bw), ("距离向采样率", m.fs, s.fs),
                       ("PRF", m.prf, s.prf)):
        if a and b and _rel(a, b) > 0.01:
            hard.append("%s相差 %.1f%%" % (name, 100 * _rel(a, b)))
    ov = dm.footprint_overlap(m, s)
    ov_min = min(ov[1], ov[2])
    if ov_min < 0.30:
        hard.append("足迹重叠仅占较小一景的 %.0f%%" % (100 * ov_min))
    g = {}
    try:
        g = dm.pair_geometry(m, s)
    except ValueError as e:
        hard.append("几何计算失败: %s" % e)
    if g:
        if abs(g["d_look"]) > 0.3:
            hard.append("视角差 %.2f° > 0.3°" % g["d_look"])
        if abs(g["d_heading"]) > 1.0:
            hard.append("航向差 %.2f° > 1°" % g["d_heading"])
    grade, score = "X", 0.0
    if g:
        bp, bc, dt = abs(g["b_perp"]), g["b_crit"], abs(g["dt_days"])
        geo = max(0.0, 1.0 - bp / bc)
        dopp = max(0.0, 1.0 - g["dc_decorr"])
        temporal = math.exp(-dt / TAU_DAYS)
        score = geo * dopp * temporal * ov_min
        if bp > 300 or bp > 0.10 * bc:
            soft.append("垂直基线 %.0f m 偏大 (临界基线 %.0f m)" % (bp, bc))
        if bp > 600 or bp > 0.25 * bc:
            hard.append("垂直基线 %.0f m 过大" % bp)
        if dt > 16:
            soft.append("时间基线 %.0f 天较长" % dt)
        if ov_min < 0.70:
            soft.append("重叠 %.0f%% 偏低" % (100 * ov_min))
        if g["dc_decorr"] > 0.10:
            soft.append("多普勒中心差占方位带宽 %.0f%%" % (100 * g["dc_decorr"]))
        if not hard:
            grade = "B" if soft else "A"
    return dict(master=m.acq_id, slave=s.acq_id, master_name=m.name, slave_name=s.name,
                master_date=m.date, slave_date=s.date, mode=m.mode, grade=grade, score=score,
                hard=hard, soft=soft, overlap_km2=ov[0], overlap_frac_master=ov[1], overlap_frac_slave=ov[2], **g)


def select(root, max_dt=40.0, include_all=False):
    prods = load_products(root)
    pairs = []
    for i in range(len(prods)):
        for j in range(i + 1, len(prods)):
            m, s = prods[i], prods[j]
            if (s.t_mid - m.t_mid) / 86400.0 > max_dt:
                continue
            r = evaluate_pair(m, s)
            if include_all or r["grade"] != "X":
                pairs.append(r)
    pairs.sort(key=lambda r: (r["grade"], -r["score"]))
    return prods, pairs


def describe(p):
    return "%d | %s %s | %s | %s %s | 视角 %.2f° | 航向 %.1f° | %d×%d | 重访轨道 %s" % (
        p.acq_id, p.date, p.name.split("_")[1], p.orbit_dir[:3], p.look_side, p.satellite,
        p.look_angle, p.heading, p.n_az, p.n_rg, p.orbit_level)


def to_markdown(prods, pairs):
    lines = ["# 电建一号像对筛选结果", "", "## 场景 (%d 景)" % len(prods), "",
             "| 编号 | 日期 | 模式 | 轨道 | 侧视 | 视角 | 航向 | 行×列 | 距离带宽 | 中心 (纬,经) | 轨道等级 |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for p in prods:
        lat, lon, _ = p.center_llh
        lines.append("| %d | %s | %s | %s | %s | %.2f° | %.1f° | %d×%d | %.0f MHz | %.3f, %.3f | %s |" % (
            p.acq_id, p.date, p.mode, p.orbit_dir[:3], p.look_side, p.look_angle, p.heading, p.n_az, p.n_rg,
            p.bw / 1e6, lat, lon, p.orbit_level))
    lines += ["", "## 可用像对 (按评级和先验相干性排序)", "",
              "| 评级 | 主影像 | 从影像 | 模式 | 间隔(天) | B⊥(m) | 临界基线(m) | 重叠(km²) | 重叠比例 | 先验相干性 | 预测偏移 距离/方位(px) | 备注 |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in pairs:
        note = "; ".join(r["hard"] + r["soft"]) or "-"
        lines.append("| %s | %d (%s) | %d (%s) | %s | %.1f | %.1f | %.0f | %.1f | %.0f%% | %.2f | %.0f / %.0f | %s |" % (
            r["grade"], r["master"], r["master_date"], r["slave"], r["slave_date"], r["mode"], r["dt_days"],
            r["b_perp"], r["b_crit"], r["overlap_km2"], 100 * min(r["overlap_frac_master"], r["overlap_frac_slave"]),
            r["score"], r["range_offset_px"], r["azimuth_offset_px"], note))
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description="电建一号干涉像对筛选")
    ap.add_argument("root")
    ap.add_argument("--max-dt", type=float, default=40.0, help="最大时间基线 (天)")
    ap.add_argument("--all", action="store_true", help="同时列出评级为 X 的组合")
    ap.add_argument("--json")
    ap.add_argument("--md")
    args = ap.parse_args(argv)
    prods, pairs = select(args.root, args.max_dt, args.all)
    if not prods:
        print("在 %s 下没有找到 *_SLC_*.XML" % args.root, file=sys.stderr)
        return 2
    md = to_markdown(prods, pairs)
    print(md)
    if args.md:
        with open(args.md, "w") as f:
            f.write(md)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(pairs, f, ensure_ascii=False, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
