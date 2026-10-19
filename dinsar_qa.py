#!/usr/bin/env python3
"""Generate lightweight QA metrics and panels for a DInSAR product."""

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

try:
    from osgeo import gdal, osr
except ImportError:
    gdal = None
    osr = None


THRESHOLDS = {
    "coherence": {
        "default_corr_thresh": 0.3,
        "coverage_warn_below": 0.5,
        "coverage_fail_below": 0.2,
    },
    "congruence_max_rad": 0.01,
    "unwrap_discontinuity": {
        "high_corr": 0.5,
        "low_corr": 0.3,
        "warn_above": 0.01,
        "fail_above": 0.05,
        "jump_rad": math.pi,
    },
    "tile_seam": {
        "warn_ratio_above": 1.5,
        "fail_ratio_above": 3.0,
        "corr_min": 0.3,
        "minimum_seam_fraction": 0.005,
        "reference_fraction_floor": 1e-4,
        "reference_offsets": [5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    },
    "edge_strip": {
        "invalid_line_fraction_min": 0.9,
        "residue_density_above": 0.2,
        "noisy_gap_tolerance_lines": 2,
        "warn_width_fraction_above": 0.01,
    },
    "dem": {
        "min_elevation_m": -500.0,
        "max_elevation_m": 9000.0,
        "warn_invalid_fraction_above": 0.01,
        "fail_invalid_fraction_above": 0.10,
    },
    "snaphu_config": {
        "lambda_relative_warn_above": 0.05,
        "near_range_relative_warn_above": 0.10,
    },
}

PI32 = np.float32(math.pi)
TWO_PI32 = np.float32(2.0 * math.pi)


def _check(check_id, name_zh, value, threshold, status, note_zh):
    return {
        "id": check_id,
        "name_zh": name_zh,
        "value": value,
        "threshold": threshold,
        "status": status,
        "note_zh": note_zh,
    }


def _read_grid(path):
    path = Path(path)
    info = subprocess.run(
        ["gmt", "grdinfo", "-C", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip().split()
    if len(info) < 11:
        raise ValueError(f"无法解析 GMT 网格尺寸: {path}")
    nx, ny = int(float(info[9])), int(float(info[10]))
    raw = subprocess.run(
        ["gmt", "grd2xyz", str(path), "-ZTLf"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    values = np.frombuffer(raw, dtype=np.float32)
    if values.size != nx * ny:
        raise ValueError(f"{path} 数据点数异常: {values.size} != {nx}x{ny}")
    return values.reshape((ny, nx))


def _wrap_pi(values):
    return np.remainder(values + PI32, TWO_PI32) - PI32


def _finite_percentile(values, percentile):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    return float(np.percentile(values, percentile)) if values.size else None


def _coreg_check(project_dir):
    report_path = Path(project_dir) / "SLC" / "coreg_report.json"
    if not report_path.is_file():
        return _check(
            "coreg",
            "配准质量",
            None,
            {"estimate_gate": "passed", "verify_gate": "passed"},
            "SKIP",
            "未找到 coreg_report.json，跳过配准门限检查。",
        )
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        estimate_gate = report.get("gate", {})
        verify = report.get("verify", {})
        verify_gate = verify.get("gate", {})
        verify_check = verify.get("check", {})
        estimate_passed = estimate_gate.get("passed") is True
        verify_passed = verify_gate.get("passed") is True
        robust_r = verify_check.get("robust_sigma_r", verify_check.get("sigma_r"))
        robust_a = verify_check.get("robust_sigma_a", verify_check.get("sigma_a"))
        value = {
            "estimate_gate_passed": estimate_passed,
            "verify_gate_passed": verify_passed,
            "verify_valid": verify_check.get("n_valid"),
            "verify_total": verify_check.get("n_points"),
            "robust_sigma_range_px": robust_r,
            "robust_sigma_azimuth_px": robust_a,
            "frac_within_0p3": verify_check.get("frac_within_0p3"),
        }
        passed = estimate_passed and verify_passed
        note = "估计与独立复核门限均通过。" if passed else "估计或独立复核门限未通过。"
        return _check(
            "coreg",
            "配准质量",
            value,
            {"estimate_gate": "passed", "verify_gate": "passed"},
            "PASS" if passed else "FAIL",
            note,
        )
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        return _check(
            "coreg",
            "配准质量",
            {"path": str(report_path), "error": str(exc)},
            {"estimate_gate": "passed", "verify_gate": "passed"},
            "FAIL",
            "配准报告无法读取或结构无效。",
        )


def _coherence_check(corr, corr_thresh):
    valid = np.isfinite(corr)
    samples = corr[valid]
    if not samples.size:
        return _check(
            "coherence",
            "相干性",
            {"valid_pixels": 0},
            {"coverage_warn_below": 0.5, "coverage_fail_below": 0.2},
            "FAIL",
            "相干性网格没有有效像元。",
        )
    coverage = float(np.count_nonzero(samples >= corr_thresh) / samples.size)
    frac_03 = float(np.count_nonzero(samples >= 0.3) / samples.size)
    frac_05 = float(np.count_nonzero(samples >= 0.5) / samples.size)
    counts, edges = np.histogram(samples, bins=10, range=(0.0, 1.0))
    if coverage < THRESHOLDS["coherence"]["coverage_fail_below"]:
        status = "FAIL"
    elif coverage < THRESHOLDS["coherence"]["coverage_warn_below"]:
        status = "WARN"
    else:
        status = "PASS"
    value = {
        "valid_pixels": int(samples.size),
        "mean": float(np.mean(samples)),
        "frac_ge_corr_thresh": coverage,
        "corr_thresh": float(corr_thresh),
        "frac_ge_0p3": frac_03,
        "frac_ge_0p5": frac_05,
        "histogram_10_bins_0_to_1": counts.astype(int).tolist(),
        "histogram_edges": edges.tolist(),
    }
    note = (
        "相干性覆盖率按 corr≥阈值统计；本产品估计量为 |<I>|/<|I|>，"
        "相较标准估计会偏高。"
    )
    return _check(
        "coherence",
        "相干性",
        value,
        {
            "corr_thresh": float(corr_thresh),
            "coverage_warn_below": THRESHOLDS["coherence"]["coverage_warn_below"],
            "coverage_fail_below": THRESHOLDS["coherence"]["coverage_fail_below"],
        },
        status,
        note,
    )


def _residue_metrics(phase, corr):
    p00, p01 = phase[:-1, :-1], phase[:-1, 1:]
    p10, p11 = phase[1:, :-1], phase[1:, 1:]
    c00, c01 = corr[:-1, :-1], corr[:-1, 1:]
    c10, c11 = corr[1:, :-1], corr[1:, 1:]
    valid = (
        np.isfinite(p00)
        & np.isfinite(p01)
        & np.isfinite(p10)
        & np.isfinite(p11)
        & np.isfinite(c00)
        & np.isfinite(c01)
        & np.isfinite(c10)
        & np.isfinite(c11)
    )
    loop = _wrap_pi(p01 - p00)
    loop += _wrap_pi(p11 - p01)
    loop += _wrap_pi(p10 - p11)
    loop += _wrap_pi(p00 - p10)
    residues = np.rint(loop / TWO_PI32).astype(np.int8)
    residue_mask = (residues != 0) & valid
    corr_mean = (c00 + c01 + c10 + c11) * np.float32(0.25)
    classes = (
        ("0_to_0p3", (corr_mean >= 0.0) & (corr_mean < 0.3)),
        ("0p3_to_0p5", (corr_mean >= 0.3) & (corr_mean < 0.5)),
        ("0p5_to_1", (corr_mean >= 0.5) & (corr_mean <= 1.0)),
    )

    def density(mask):
        mask = mask & valid
        count = int(np.count_nonzero(mask))
        return {
            "residues": int(np.count_nonzero(residue_mask & mask)),
            "valid_loops": count,
            "density": float(np.count_nonzero(residue_mask & mask) / count) if count else None,
        }

    value = {"overall": density(valid)}
    for name, mask in classes:
        value[name] = density(mask)
    return value, residue_mask, valid


def _residue_check(residue_value):
    return _check(
        "residues",
        "相位残差点",
        residue_value,
        None,
        "INFO",
        "残差点由 2×2 相位闭合环路计算，报告总体及相干性分档密度。",
    )


def _congruence_check(phase, unwrap):
    valid = np.isfinite(phase) & np.isfinite(unwrap)
    if not np.any(valid):
        maximum = None
    else:
        difference = _wrap_pi(unwrap[valid] - phase[valid])
        maximum = float(np.max(np.abs(difference)))
    limit = THRESHOLDS["congruence_max_rad"]
    return _check(
        "congruence",
        "解缠相位一致性",
        {"max_abs_wrapped_difference_rad": maximum, "valid_pixels": int(np.count_nonzero(valid))},
        {"max_abs_rad": limit},
        "FAIL" if maximum is None or maximum > limit else "PASS",
        "比较 angle(exp(i·(unwrap−phasefilt))) 的最大绝对值。",
    )


def _unwrap_discontinuity_check(unwrap, corr):
    high_limit = THRESHOLDS["unwrap_discontinuity"]["high_corr"]
    low_limit = THRESHOLDS["unwrap_discontinuity"]["low_corr"]
    jump_limit = THRESHOLDS["unwrap_discontinuity"]["jump_rad"]
    counts = {"high": [0, 0], "low": [0, 0]}
    bad_pixels = np.zeros(unwrap.shape, dtype=bool)
    for first, second, cf, cs, locations in (
        (unwrap[:, :-1], unwrap[:, 1:], corr[:, :-1], corr[:, 1:], (slice(None), slice(1, None))),
        (unwrap[:-1, :], unwrap[1:, :], corr[:-1, :], corr[1:, :], (slice(1, None), slice(None))),
    ):
        valid = np.isfinite(first) & np.isfinite(second) & np.isfinite(cf) & np.isfinite(cs)
        jump = valid & (np.abs(second - first) > jump_limit)
        high = valid & (cf >= high_limit) & (cs >= high_limit)
        low = valid & (cf < low_limit) & (cs < low_limit)
        counts["high"][0] += int(np.count_nonzero(jump & high))
        counts["high"][1] += int(np.count_nonzero(high))
        counts["low"][0] += int(np.count_nonzero(jump & low))
        counts["low"][1] += int(np.count_nonzero(low))
        bad_pixels[locations] |= jump

    high_fraction = counts["high"][0] / counts["high"][1] if counts["high"][1] else None
    low_fraction = counts["low"][0] / counts["low"][1] if counts["low"][1] else None
    if high_fraction is None:
        status = "WARN"
    elif high_fraction > THRESHOLDS["unwrap_discontinuity"]["fail_above"]:
        status = "FAIL"
    elif high_fraction > THRESHOLDS["unwrap_discontinuity"]["warn_above"]:
        status = "WARN"
    else:
        status = "PASS"
    value = {
        "high_corr_fraction": high_fraction,
        "high_corr_jumps": counts["high"][0],
        "high_corr_neighbor_pairs": counts["high"][1],
        "low_corr_fraction": low_fraction,
        "low_corr_jumps": counts["low"][0],
        "low_corr_neighbor_pairs": counts["low"][1],
    }
    return (
        _check(
            "unwrap_discontinuity",
            "解缠不连续",
            value,
            {
                "high_corr": high_limit,
                "jump_rad": jump_limit,
                "warn_fraction_above": THRESHOLDS["unwrap_discontinuity"]["warn_above"],
                "fail_fraction_above": THRESHOLDS["unwrap_discontinuity"]["fail_above"],
            },
            status,
            "统计高相干（corr≥0.5）和低相干（corr<0.3）邻接像元间 |Δunwrap|>π 的比例。",
        ),
        bad_pixels,
    )


def _line_jump_fraction(first, second, eligible=None):
    valid = np.isfinite(first) & np.isfinite(second)
    if eligible is not None:
        valid &= eligible
    count = int(np.count_nonzero(valid))
    fraction = float(np.count_nonzero(valid & (np.abs(second - first) > math.pi)) / count) if count else 0.0
    return fraction, count


def _seam_axis(unwrap, corr, tile_count, axis):
    size = unwrap.shape[axis]
    seams = [int(k * size // tile_count) for k in range(1, tile_count)]
    seams = [seam for seam in seams if 0 < seam < size]
    thresholds = THRESHOLDS["tile_seam"]

    def line_pair_stat(index):
        if axis == 0:
            first, second = unwrap[index - 1, :], unwrap[index, :]
            first_corr, second_corr = corr[index - 1, :], corr[index, :]
        else:
            first, second = unwrap[:, index - 1], unwrap[:, index]
            first_corr, second_corr = corr[:, index - 1], corr[:, index]
        eligible = (first_corr >= thresholds["corr_min"]) & (second_corr >= thresholds["corr_min"])
        return _line_jump_fraction(first, second, eligible)

    seam_values = []
    for seam in seams:
        fraction, neighbor_pairs = line_pair_stat(seam)
        reference_samples = []
        for offset in thresholds["reference_offsets"]:
            for signed_offset in (-offset, offset):
                index = seam + signed_offset
                if 0 < index < size:
                    reference_fraction, reference_pairs = line_pair_stat(index)
                    reference_samples.append(
                        {
                            "index": index,
                            "offset": signed_offset,
                            "fraction": reference_fraction,
                            "neighbor_pairs": reference_pairs,
                        }
                    )
        local_reference = (
            float(np.median([sample["fraction"] for sample in reference_samples]))
            if reference_samples
            else 0.0
        )
        ratio = float(fraction / max(local_reference, thresholds["reference_fraction_floor"]))
        if fraction < thresholds["minimum_seam_fraction"]:
            status = "PASS"
        elif ratio > thresholds["fail_ratio_above"]:
            status = "FAIL"
        elif ratio > thresholds["warn_ratio_above"]:
            status = "WARN"
        else:
            status = "PASS"
        seam_values.append(
            {
                "index": seam,
                "fraction": fraction,
                "neighbor_pairs": neighbor_pairs,
                "local_reference_fraction": local_reference,
                "local_reference_samples": reference_samples,
                "ratio": ratio,
                "status": status,
            }
        )
    reference_fraction = (
        float(np.median([item["local_reference_fraction"] for item in seam_values]))
        if seam_values
        else None
    )
    return {
        "seams": seam_values,
        "reference_fraction": reference_fraction,
    }


def _tile_seam_check(unwrap, corr, tiles, overlap):
    row_tiles, col_tiles = tiles
    rows = _seam_axis(unwrap, corr, row_tiles, 0)
    cols = _seam_axis(unwrap, corr, col_tiles, 1)
    ratios = [
        seam["ratio"]
        for axis in (rows, cols)
        for seam in axis["seams"]
        if seam["fraction"] >= THRESHOLDS["tile_seam"]["minimum_seam_fraction"]
    ]
    maximum = max(ratios) if ratios else None
    seam_statuses = [
        seam["status"]
        for axis in (rows, cols)
        for seam in axis["seams"]
    ]
    if not seam_statuses:
        status = "INFO"
    elif "FAIL" in seam_statuses:
        status = "FAIL"
    elif "WARN" in seam_statuses:
        status = "WARN"
    else:
        status = "PASS"
    value = {
        "row_seams": rows["seams"],
        "row_reference_fraction": rows["reference_fraction"],
        "column_seams": cols["seams"],
        "column_reference_fraction": cols["reference_fraction"],
        "max_seam_to_reference_ratio": maximum,
        "overlap_pixels": int(overlap),
        "tile_rows_columns": [int(row_tiles), int(col_tiles)],
    }
    return _check(
        "tile_seams",
        "分块接缝",
        value,
        {
            "warn_ratio_above": THRESHOLDS["tile_seam"]["warn_ratio_above"],
            "fail_ratio_above": THRESHOLDS["tile_seam"]["fail_ratio_above"],
            "corr_min": THRESHOLDS["tile_seam"]["corr_min"],
            "minimum_seam_fraction": THRESHOLDS["tile_seam"]["minimum_seam_fraction"],
            "reference_fraction_floor": THRESHOLDS["tile_seam"]["reference_fraction_floor"],
            "reference_offsets": THRESHOLDS["tile_seam"]["reference_offsets"],
        },
        status,
        "仅统计接缝两侧 corr≥0.3 的邻接像元，并与接缝附近 ±5 至 ±50 行/列的局部中位数比较；绝对跃变率低于 0.5% 时不告警。",
    )


def _edge_width(profile, threshold):
    width = 0
    for value in profile:
        if value is None or value < threshold:
            break
        width += 1
    return width


def _noisy_band(profile, start, threshold, max_gap):
    band_start = next(
        (
            index
            for index in range(start, min(start + max_gap + 1, len(profile)))
            if profile[index] is not None and profile[index] > threshold
        ),
        None,
    )
    if band_start is None:
        return None, 0

    width = 0
    gap = 0
    for value in profile[band_start:]:
        if value is not None and value > threshold:
            width += gap + 1
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break
    return band_start, width


def _noisy_width(profile, start, threshold, max_gap):
    return _noisy_band(profile, start, threshold, max_gap)[1]


def _edge_strip_check(residue_mask, valid_loops, phasefilt):
    thresholds = THRESHOLDS["edge_strip"]
    density_threshold = thresholds["residue_density_above"]

    def density_by_axis(axis):
        total = np.count_nonzero(valid_loops, axis=axis)
        residues = np.count_nonzero(residue_mask, axis=axis)
        values = np.divide(
            residues,
            total,
            out=np.full(total.shape, np.nan, dtype=np.float32),
            where=total > 0,
        )
        return [float(value) if np.isfinite(value) else None for value in values]

    invalid = ~np.isfinite(phasefilt) | (phasefilt == 0.0)
    row_invalid_fraction = np.mean(invalid, axis=1).tolist()
    column_invalid_fraction = np.mean(invalid, axis=0).tolist()
    row_density = density_by_axis(1)
    column_density = density_by_axis(0)
    profiles = {
        "top": (row_density, row_invalid_fraction),
        "bottom": (row_density[::-1], row_invalid_fraction[::-1]),
        "left": (column_density, column_invalid_fraction),
        "right": (column_density[::-1], column_invalid_fraction[::-1]),
    }
    dimensions = {
        "top": residue_mask.shape[0],
        "bottom": residue_mask.shape[0],
        "left": residue_mask.shape[1],
        "right": residue_mask.shape[1],
    }
    widths = {}
    fractions = {}
    band_starts = {}
    for edge, (density_profile, invalid_profile) in profiles.items():
        invalid_width = _edge_width(invalid_profile, thresholds["invalid_line_fraction_min"])
        noisy_start, noisy_width = _noisy_band(
            density_profile,
            invalid_width,
            density_threshold,
            thresholds["noisy_gap_tolerance_lines"],
        )
        dimension = dimensions[edge]
        widths[edge] = {
            "invalid_width": invalid_width,
            "noisy_width": noisy_width,
        }
        band_starts[edge] = noisy_start
        fractions[edge] = {
            "invalid_fraction": invalid_width / dimension if dimension else 0.0,
            "noisy_fraction": noisy_width / dimension if dimension else 0.0,
        }
    maximum = max((item["noisy_fraction"] for item in fractions.values()), default=0.0)
    status = (
        "WARN"
        if maximum > thresholds["warn_width_fraction_above"]
        else "PASS"
    )
    return _check(
        "edge_strip",
        "边缘残差带",
        {
            "width_pixels": widths,
            "width_fraction": fractions,
            "noisy_band_start_offset": band_starts,
        },
        {
            "invalid_line_fraction_min": thresholds["invalid_line_fraction_min"],
            "residue_density_above": density_threshold,
            "noisy_gap_tolerance_lines": thresholds["noisy_gap_tolerance_lines"],
            "warn_width_fraction_above": thresholds["warn_width_fraction_above"],
        },
        status,
        "记录边缘零值或无效填充带宽（本身不触发告警）；远/近距向噪声带通常表示从影像在该范围无重叠（可能由距离向偏移造成），应予掩膜。",
    )


def _product_mask_check(intf_dir, shape, edge_check):
    intf_dir = Path(intf_dir)
    mask_path = intf_dir / "product_mask.grd"
    if not mask_path.is_file():
        return None
    mask = _read_grid(mask_path)
    if mask.shape != shape:
        raise ValueError(f"product_mask 网格尺寸 {mask.shape} 与 corr {shape} 不一致")
    valid = np.isfinite(mask) & (mask > 0.5)
    summary_path = intf_dir / "product_mask_summary.json"
    summary = {}
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            summary = {}
    fractions = summary.get("fractions", {})
    edge_widths = edge_check["value"]["width_pixels"]
    edge_starts = edge_check["value"].get("noisy_band_start_offset", {})
    coverage = {}
    ny, nx = shape
    for edge in ("top", "bottom", "left", "right"):
        width = int(edge_widths[edge]["noisy_width"])
        start = edge_starts.get(edge)
        if width < 1 or start is None:
            coverage[edge] = {"pixels": 0, "masked_fraction": None}
            continue
        start = int(start)
        if edge == "top":
            region = valid[start : start + width, :]
        elif edge == "bottom":
            region = valid[ny - start - width : ny - start, :]
        elif edge == "left":
            region = valid[:, start : start + width]
        else:
            region = valid[:, nx - start - width : nx - start]
        count = int(region.size)
        coverage[edge] = {
            "pixels": count,
            "masked_fraction": float(np.count_nonzero(~region) / count) if count else None,
        }
    return _check(
        "product_mask",
        "产品有效性掩膜",
        {
            "masked_fraction_by_overlap": fractions.get("masked_fraction_by_overlap"),
            "masked_fraction_by_corr": fractions.get("masked_fraction_by_corr"),
            "masked_fraction_total": float(np.count_nonzero(~valid) / valid.size),
            "noisy_band_coverage_by_edge": coverage,
        },
        None,
        "INFO",
        "报告影像覆盖/相干性掩膜比例，以及各边缘残差噪声带中被掩膜的像元比例。",
    )


def _detrend_check(ramp, unwrap, detrended, corr, wavelength):
    finite_ramp = ramp[np.isfinite(ramp)]
    ramp_ptp = float(np.max(finite_ramp) - np.min(finite_ramp)) if finite_ramp.size else None
    ramp_mm = ramp_ptp * wavelength / (4.0 * math.pi) * 1000.0 if ramp_ptp is not None else None
    mask = (corr >= 0.5) & np.isfinite(corr) & np.isfinite(unwrap) & np.isfinite(detrended)
    value = {
        "orbital_ramp_peak_to_peak_rad": ramp_ptp,
        "orbital_ramp_peak_to_peak_mm": ramp_mm,
        "unwrap_std_corr_ge_0p5_rad": float(np.std(unwrap[mask])) if np.any(mask) else None,
        "unwrap_detrended_std_corr_ge_0p5_rad": float(np.std(detrended[mask])) if np.any(mask) else None,
        "valid_high_corr_pixels": int(np.count_nonzero(mask)),
    }
    return _check(
        "detrend",
        "去趋势效果",
        value,
        None,
        "INFO",
        "报告轨道坡度峰峰值及高相干区解缠前后标准差；仅作诊断。",
    )


def _los_check(los, corr):
    finite = np.isfinite(los)
    high = finite & np.isfinite(corr) & (corr >= 0.5)
    samples = los[high]
    all_values = los[finite]
    if samples.size:
        median = float(np.median(samples))
        mad = float(np.median(np.abs(samples - median)))
        robust_std = 1.4826 * mad
        percentiles = {
            str(p): float(value)
            for p, value in zip((1, 5, 25, 50, 75, 95, 99), np.percentile(samples, (1, 5, 25, 50, 75, 95, 99)))
        }
        standard_deviation = float(np.std(samples))
        outlier_fraction = (
            float(np.count_nonzero(finite & (np.abs(los) > 5.0 * robust_std)) / los.size)
            if robust_std > 0
            else 0.0
        )
    else:
        robust_std = None
        percentiles = {}
        standard_deviation = None
        outlier_fraction = None
    block_medians = []
    if samples.size:
        y_edges = np.linspace(0, los.shape[0], 9, dtype=int)
        x_edges = np.linspace(0, los.shape[1], 9, dtype=int)
        for y0, y1 in zip(y_edges[:-1], y_edges[1:]):
            for x0, x1 in zip(x_edges[:-1], x_edges[1:]):
                block_mask = high[y0:y1, x0:x1]
                if np.count_nonzero(block_mask) > 500:
                    block_medians.append(float(np.median(los[y0:y1, x0:x1][block_mask])))
    value = {
        "corr_ge_0p5_percentiles_mm": percentiles,
        "corr_ge_0p5_std_mm": standard_deviation,
        "corr_ge_0p5_robust_std_mm": robust_std,
        "overall_min_mm": float(np.min(all_values)) if all_values.size else None,
        "overall_max_mm": float(np.max(all_values)) if all_values.size else None,
        "fraction_all_pixels_abs_gt_5_robust_std": outlier_fraction,
        "block_medians_8x8_mm": block_medians,
        "block_median_std_mm": float(np.std(block_medians)) if block_medians else None,
        "valid_corr_ge_0p5_pixels": int(samples.size),
    }
    return _check(
        "los",
        "视线向位移",
        value,
        {"corr_min": 0.5, "block_grid": "8x8", "min_block_valid_pixels": 501},
        "INFO",
        "位移以 mm 表示；单景对的精度受此级别长波残差（大气/轨道）限制。",
    )


def _dem_check(dem):
    invalid = (
        ~np.isfinite(dem)
        | (dem < THRESHOLDS["dem"]["min_elevation_m"])
        | (dem > THRESHOLDS["dem"]["max_elevation_m"])
    )
    fraction = float(np.count_nonzero(invalid) / dem.size) if dem.size else 1.0
    if fraction > THRESHOLDS["dem"]["fail_invalid_fraction_above"]:
        status = "FAIL"
    elif fraction > THRESHOLDS["dem"]["warn_invalid_fraction_above"]:
        status = "WARN"
    else:
        status = "PASS"
    return _check(
        "dem",
        "DEM 有效性",
        {"invalid_fraction": fraction, "invalid_pixels": int(np.count_nonzero(invalid)), "total_pixels": int(dem.size)},
        {
            "valid_elevation_m": [
                THRESHOLDS["dem"]["min_elevation_m"],
                THRESHOLDS["dem"]["max_elevation_m"],
            ],
            "warn_invalid_fraction_above": THRESHOLDS["dem"]["warn_invalid_fraction_above"],
            "fail_invalid_fraction_above": THRESHOLDS["dem"]["fail_invalid_fraction_above"],
        },
        status,
        "无效像元包括 NaN、低于 -500 m 或高于 9000 m。",
    )


def _epsg_4326(wkt):
    if not wkt:
        return False
    if osr is not None:
        reference = osr.SpatialReference()
        if reference.ImportFromWkt(wkt) == 0:
            reference.AutoIdentifyEPSG()
            code = reference.GetAuthorityCode(None)
            if code is None:
                code = reference.GetAttrValue("GEOGCS|AUTHORITY", 1)
            if code is None:
                code = reference.GetAttrValue("GEOGCRS|ID", 1)
            return code == "4326"
    return bool(re.search(r'(?:AUTHORITY|ID)\s*\[\s*"EPSG"\s*,\s*"?4326"?', wkt, re.IGNORECASE))


def _gdal_tiff(path):
    if gdal is not None:
        dataset = gdal.Open(str(path), gdal.GA_ReadOnly)
        if dataset is None:
            raise ValueError("GDAL 无法读取")
        band = dataset.GetRasterBand(1)
        if band is None:
            raise ValueError("GeoTIFF 不包含栅格波段")
        nodata = band.GetNoDataValue()
        total = dataset.RasterXSize * dataset.RasterYSize
        invalid_count = 0
        for y0 in range(0, dataset.RasterYSize, 512):
            array = band.ReadAsArray(0, y0, dataset.RasterXSize, min(512, dataset.RasterYSize - y0))
            invalid = ~np.isfinite(array)
            if nodata is not None:
                invalid |= array == nodata
            invalid_count += int(np.count_nonzero(invalid))
        result = {
            "readable": True,
            "width": int(dataset.RasterXSize),
            "height": int(dataset.RasterYSize),
            "epsg_4326": _epsg_4326(dataset.GetProjection()),
            "nodata_fraction": float(invalid_count / total) if total else None,
            "nodata_value": nodata,
        }
        dataset = None
        return result

    result = subprocess.run(
        ["gdalinfo", "-json", "-stats", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    metadata = json.loads(result.stdout)
    wkt = metadata.get("coordinateSystem", {}).get("wkt", "")
    band_metadata = metadata.get("bands", [{}])[0].get("metadata", {}).get("", {})
    valid_percent = band_metadata.get("STATISTICS_VALID_PERCENT")
    return {
        "readable": True,
        "width": metadata.get("size", [None, None])[0],
        "height": metadata.get("size", [None, None])[1],
        "epsg_4326": _epsg_4326(wkt),
        "nodata_fraction": 1.0 - float(valid_percent) / 100.0 if valid_percent is not None else None,
        "nodata_value": metadata.get("bands", [{}])[0].get("noDataValue"),
    }


def _geotiff_check(project_dir):
    names = ("unwrap_detrended", "corr", "final_data_values")
    results = {}
    failures = []
    for name in names:
        path = Path(project_dir) / "results" / f"{name}.tif"
        if not path.is_file():
            results[name] = {"readable": False, "error": "missing"}
            failures.append(f"{name}: missing")
            continue
        try:
            results[name] = _gdal_tiff(path)
            if not results[name]["epsg_4326"]:
                failures.append(f"{name}: CRS is not EPSG:4326")
        except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
            results[name] = {"readable": False, "error": str(exc)}
            failures.append(f"{name}: unreadable")
    status = "FAIL" if failures else "PASS"
    note = "；".join(failures) if failures else "三个 GeoTIFF 均可读且 CRS 为 EPSG:4326。"
    return _check(
        "geotiff",
        "地理编码 GeoTIFF",
        results,
        {"required": [f"{name}.tif" for name in names], "epsg": 4326},
        status,
        note,
    )


def _read_prm_value(path, key):
    key = key.lower()
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        name, value = line.split("=", 1)
        if name.strip().lower() == key:
            try:
                return float(value.strip().split()[0])
            except (ValueError, IndexError):
                return None
    return None


def _snaphu_config_check(config_path, master_prm, wavelength):
    if not config_path:
        return _check(
            "snaphu_config",
            "SNAPHU 参数",
            None,
            {
                "lambda_relative_warn_above": THRESHOLDS["snaphu_config"]["lambda_relative_warn_above"],
                "near_range_relative_warn_above": THRESHOLDS["snaphu_config"]["near_range_relative_warn_above"],
            },
            "SKIP",
            "未提供 SNAPHU 配置文件。",
        )
    config_path = Path(config_path)
    if not config_path.is_file():
        return _check(
            "snaphu_config",
            "SNAPHU 参数",
            {"path": str(config_path)},
            None,
            "WARN",
            "SNAPHU 配置文件不存在。",
        )
    config = {}
    for line in config_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^\s*(LAMBDA|NEARRANGE)\s+([-+0-9.eE]+)", line, re.IGNORECASE)
        if match:
            config[match.group(1).upper()] = float(match.group(2))
    prm_near_range = _read_prm_value(master_prm, "near_range")
    config_lambda = config.get("LAMBDA")
    config_near = config.get("NEARRANGE")
    lambda_relative = abs(config_lambda - wavelength) / abs(wavelength) if config_lambda and wavelength else None
    near_relative = (
        abs(config_near - prm_near_range) / abs(prm_near_range)
        if config_near is not None and prm_near_range
        else None
    )
    warning = (
        lambda_relative is None
        or near_relative is None
        or lambda_relative > THRESHOLDS["snaphu_config"]["lambda_relative_warn_above"]
        or near_relative > THRESHOLDS["snaphu_config"]["near_range_relative_warn_above"]
    )
    value = {
        "config_path": str(config_path),
        "config_lambda_m": config_lambda,
        "prm_wavelength_m": float(wavelength),
        "lambda_relative_difference": lambda_relative,
        "config_near_range_m": config_near,
        "prm_near_range_m": prm_near_range,
        "near_range_relative_difference": near_relative,
    }
    thresholds = {
        "lambda_relative_warn_above": THRESHOLDS["snaphu_config"]["lambda_relative_warn_above"],
        "near_range_relative_warn_above": THRESHOLDS["snaphu_config"]["near_range_relative_warn_above"],
    }
    note = "参数与主影像 PRM 不匹配或缺失。" if warning else "LAMBDA 和 NEARRANGE 与主影像参数一致。"
    return _check("snaphu_config", "SNAPHU 参数", value, thresholds, "WARN" if warning else "PASS", note)


def _palette(values, stops, valid):
    clipped = np.clip(np.where(valid, values, 0.0), 0.0, 1.0)
    position = clipped * (len(stops) - 1)
    low = np.floor(position).astype(np.int16)
    high = np.minimum(low + 1, len(stops) - 1)
    fraction = (position - low)[..., None]
    stop_array = np.asarray(stops, dtype=np.float32)
    rgb = stop_array[low] * (1.0 - fraction) + stop_array[high] * fraction
    rgb = np.clip(rgb, 0, 255).astype(np.uint8)
    rgb[~valid] = (155, 155, 155)
    return rgb


def _panels(out_path, corr, phase, unwrap, los, residue_mask, bad_pixels):
    limit = 600
    step = max(1, int(math.ceil(max(corr.shape) / limit)))
    corr_small = corr[::step, ::step]
    phase_small = phase[::step, ::step]
    unwrap_small = unwrap[::step, ::step]
    los_small = los[::step, ::step]
    residue_pixels = np.zeros(corr.shape, dtype=bool)
    residue_pixels[:-1, :-1] = residue_mask
    bad_small = bad_pixels[::step, ::step] | residue_pixels[::step, ::step]
    height, width = corr_small.shape
    panel_height = height + 76
    panel_width = 600

    coherence_colors = [(68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98), (253, 231, 37)]
    coherence_rgb = _palette(corr_small, coherence_colors, np.isfinite(corr_small))

    phase_valid = np.isfinite(phase_small)
    phase_filled = np.where(phase_valid, phase_small, 0.0)
    hue = ((phase_filled + math.pi) / (2.0 * math.pi) * 255.0).astype(np.uint8)
    hsv = np.zeros((height, width, 3), dtype=np.uint8)
    hsv[:, :, 0] = hue
    hsv[:, :, 1] = 255
    hsv[:, :, 2] = 255
    phase_rgb = np.asarray(Image.fromarray(hsv, mode="HSV").convert("RGB"))
    phase_rgb[~phase_valid] = (155, 155, 155)

    unwrap_rgb = np.empty((height, width, 3), dtype=np.uint8)
    unwrap_rgb[:] = (190, 190, 190)
    unwrap_rgb[~np.isfinite(unwrap_small)] = (120, 120, 120)
    unwrap_rgb[bad_small] = (220, 35, 35)

    los_valid = np.isfinite(los_small) & np.isfinite(corr_small) & (corr_small >= 0.5)
    full_los_valid = np.isfinite(los) & np.isfinite(corr) & (corr >= 0.5)
    scale = _finite_percentile(np.abs(los[full_los_valid]), 99) or 1.0
    normalized = np.clip(los_small / scale, -1.0, 1.0)
    los_rgb = np.empty((height, width, 3), dtype=np.uint8)
    negative = los_valid & (normalized < 0)
    positive = los_valid & (normalized >= 0)
    fraction = np.abs(normalized)
    los_rgb[:] = (155, 155, 155)
    los_rgb[negative] = np.stack(
        (
            (255.0 * (1.0 - fraction[negative])).astype(np.uint8),
            (255.0 * (1.0 - fraction[negative])).astype(np.uint8),
            np.full(np.count_nonzero(negative), 255, dtype=np.uint8),
        ),
        axis=-1,
    )
    los_rgb[positive] = np.stack(
        (
            np.full(np.count_nonzero(positive), 255, dtype=np.uint8),
            (255.0 * (1.0 - fraction[positive])).astype(np.uint8),
            (255.0 * (1.0 - fraction[positive])).astype(np.uint8),
        ),
        axis=-1,
    )
    los_rgb[~los_valid] = (155, 155, 155)

    panels = (
        ("Coherence corr", coherence_rgb, "0        0.5        1"),
        ("Wrapped phasefilt", phase_rgb, "-pi                 0                +pi"),
        ("Unwrap jumps / residues", unwrap_rgb, "Gray: continuous/invalid   Red: anomaly"),
        (f"LOS displacement (corr>=0.5, +/-{scale:.2f} mm)", los_rgb, f"-{scale:.2f}        0        +{scale:.2f} mm"),
    )
    canvas = Image.new("RGB", (panel_width * 4, panel_height), "white")
    draw = ImageDraw.Draw(canvas)
    for index, (title, rgb, legend) in enumerate(panels):
        x0 = index * panel_width
        draw.text((x0 + 8, 5), title, fill="black")
        image = Image.fromarray(rgb, mode="RGB").resize((width, height), Image.NEAREST)
        canvas.paste(image, (x0, 28))
        draw.text((x0 + 8, 31 + height), legend, fill="black")
    canvas.save(out_path)


def _markdown_value(check):
    value = check["value"]
    if check["id"] == "coherence" and value:
        return (
            f"均值 {value['mean']:.4f}；覆盖率 {value['frac_ge_corr_thresh']:.3%}；"
            f"≥0.3 {value['frac_ge_0p3']:.3%}；≥0.5 {value['frac_ge_0p5']:.3%}；"
            f"直方图 {value['histogram_10_bins_0_to_1']}"
        )
    if check["id"] == "congruence" and value:
        number = value["max_abs_wrapped_difference_rad"]
        return f"最大包裹差 {number:.6g} rad" if number is not None else "无有效像元"
    if check["id"] == "coreg" and value:
        if "verify_valid" in value:
            return (
                f"estimate gate={value['estimate_gate_passed']}；verify gate={value['verify_gate_passed']}；"
                f"有效点 {value['verify_valid']}/{value['verify_total']}；"
                f"σ={value['robust_sigma_range_px']}/{value['robust_sigma_azimuth_px']} px；"
                f"≤0.3px={value['frac_within_0p3']}"
            )
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _write_markdown(report, path):
    lines = [
        "# DInSAR 产品质量检查",
        "",
        "| 检查项 | 数值 | 阈值 | 结论 | 说明 |",
        "|---|---|---|---|---|",
    ]
    for check in report["checks"]:
        value = _markdown_value(check).replace("|", "\\|").replace("\n", " ")
        threshold = json.dumps(check["threshold"], ensure_ascii=False, separators=(",", ":"))
        threshold = threshold.replace("|", "\\|")
        note = check["note_zh"].replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {check['name_zh']} | {value} | {threshold} | {check['status']} | {note} |")
    lines.extend(("", f"**总体结论：{report['overall']}**", ""))
    path.write_text("\n".join(lines), encoding="utf-8")


def _run(args):
    intf_dir = Path(args.intf_dir)
    project_dir = Path(args.project_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tiles = tuple(int(item) for item in args.snaphu_tiles.split())
    if len(tiles) != 2 or any(item < 1 for item in tiles):
        raise ValueError("--snaphu-tiles 必须是两个正整数: 行数 列数")
    if args.downsample < 1:
        raise ValueError("--downsample 必须为正整数")

    grids = {
        "corr": _read_grid(intf_dir / "corr.grd"),
        "phase": _read_grid(intf_dir / "phasefilt.grd"),
        "unwrap": _read_grid(intf_dir / "unwrap.grd"),
        "detrended": _read_grid(intf_dir / "unwrap_detrended.grd"),
        "los": _read_grid(intf_dir / "los_displacement.grd"),
    }
    shape = grids["corr"].shape
    for name, array in grids.items():
        if array.shape != shape:
            raise ValueError(f"{name} 网格尺寸 {array.shape} 与 corr {shape} 不一致")
    wavelength = float(args.wavelength)

    checks = []
    checks.append(_coreg_check(project_dir))
    checks.append(_coherence_check(grids["corr"], args.corr_thresh))
    residue_value, residue_mask, valid_loops = _residue_metrics(grids["phase"], grids["corr"])
    checks.append(_residue_check(residue_value))
    checks.append(_congruence_check(grids["phase"], grids["unwrap"]))
    discontinuity_check, bad_pixels = _unwrap_discontinuity_check(grids["unwrap"], grids["corr"])
    checks.append(discontinuity_check)
    checks.append(_tile_seam_check(grids["unwrap"], grids["corr"], tiles, args.snaphu_overlap))
    edge_check = _edge_strip_check(residue_mask, valid_loops, grids["phase"])
    checks.append(edge_check)
    product_mask_check = _product_mask_check(intf_dir, shape, edge_check)
    if product_mask_check is not None:
        checks.append(product_mask_check)
    ramp = _read_grid(project_dir / "debug" / "orbital_ramp.grd")
    checks.append(_detrend_check(ramp, grids["unwrap"], grids["detrended"], grids["corr"], wavelength))
    del ramp
    checks.append(_los_check(grids["los"], grids["corr"]))
    dem = _read_grid(project_dir / "topo" / "dem.grd")
    checks.append(_dem_check(dem))
    del dem
    checks.append(_geotiff_check(project_dir))
    checks.append(_snaphu_config_check(args.snaphu_conf, args.master_prm, wavelength))

    severity = {"PASS": 0, "WARN": 1, "FAIL": 2}
    overall = max(
        (check["status"] for check in checks if check["status"] in severity),
        key=severity.get,
        default="PASS",
    )
    report = {
        "metadata": {
            "intf_dir": str(intf_dir),
            "project_dir": str(project_dir),
            "master_prm": str(args.master_prm),
            "wavelength_m": wavelength,
            "downsample": args.downsample,
            "snaphu_tiles_rows_columns": list(tiles),
            "snaphu_overlap_pixels": args.snaphu_overlap,
            "grid_rows_columns": [int(shape[0]), int(shape[1])],
            "corr_thresh": args.corr_thresh,
        },
        "thresholds": THRESHOLDS,
        "checks": checks,
        "overall": overall,
    }
    json_path = out_dir / "qa_report.json"
    markdown_path = out_dir / "qa_report.md"
    panels_path = out_dir / "qa_panels.png"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_markdown(report, markdown_path)
    _panels(
        panels_path,
        grids["corr"],
        grids["phase"],
        grids["unwrap"],
        grids["los"],
        residue_mask,
        bad_pixels,
    )
    print(f"QA overall: {overall}")
    print(f"JSON: {json_path}")
    print(f"Markdown: {markdown_path}")
    print(f"Panels: {panels_path}")
    return overall


def main(argv=None):
    parser = argparse.ArgumentParser(description="DInSAR product quality checks")
    parser.add_argument("--intf-dir", required=True)
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--master-prm", required=True)
    parser.add_argument("--wavelength", required=True, type=float)
    parser.add_argument("--downsample", required=True, type=int)
    parser.add_argument("--snaphu-tiles", required=True, help="row and column tile counts, e.g. '2 2'")
    parser.add_argument("--snaphu-overlap", required=True, type=int)
    parser.add_argument("--snaphu-conf")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--corr-thresh", type=float, default=THRESHOLDS["coherence"]["default_corr_thresh"])
    args = parser.parse_args(argv)
    try:
        overall = _run(args)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"QA 执行失败: {exc}", file=sys.stderr)
        return 1
    return 2 if args.strict and overall == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
