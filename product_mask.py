#!/usr/bin/env python3
"""Build a validity mask from resampling overlap and coherence."""

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import dinsar_qa


SHIFT_KEYS = ("rshift", "ashift", "stretch_r", "a_stretch_r", "stretch_a", "a_stretch_a")
SUBPIXEL_KEYS = ("sub_int_r", "sub_int_a")


def _read_prm(path, required_keys):
    path = Path(path)
    values = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.replace("=", " ").split()
        if len(fields) < 2:
            continue
        key = fields[0]
        if key not in required_keys and key not in SUBPIXEL_KEYS:
            continue
        try:
            value = float(fields[1])
        except ValueError:
            continue
        if not math.isfinite(value):
            raise ValueError(f"PRM {path} key {key} is not finite")
        values[key] = value
    missing = [key for key in required_keys if key not in values]
    if missing:
        raise ValueError(f"missing required PRM key(s) in {path}: {', '.join(missing)}")
    defaults = [key for key in SUBPIXEL_KEYS if key not in values]
    for key in defaults:
        values[key] = 0.0
    return values, defaults


def _grid_info(path):
    result = subprocess.run(
        ["gmt", "grdinfo", "-C", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    fields = result.stdout.strip().split()
    if len(fields) < 11:
        raise ValueError(f"unable to parse GMT grid geometry: {path}")
    metadata = fields[-12:] if len(fields) >= 13 else fields[-11:] if len(fields) == 12 else fields[-10:]
    x_min, x_max, y_min, y_max = map(float, metadata[:4])
    x_inc, y_inc = map(float, metadata[6:8])
    nx, ny = int(float(metadata[8])), int(float(metadata[9]))
    registration = int(float(metadata[10])) if len(metadata) >= 11 else 0
    if nx < 1 or ny < 1 or x_inc <= 0 or y_inc <= 0:
        raise ValueError(f"invalid GMT grid geometry: {path}")
    return {
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "x_inc": x_inc,
        "y_inc": y_inc,
        "nx": nx,
        "ny": ny,
        "registration": registration,
    }


def _read_grid(path, info):
    raw = subprocess.run(
        ["gmt", "grd2xyz", str(path), "-ZTLf"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    values = np.frombuffer(raw, dtype=np.float32)
    expected = info["nx"] * info["ny"]
    if values.size != expected:
        raise ValueError(f"grid data size mismatch for {path}: {values.size} != {expected}")
    return values.reshape((info["ny"], info["nx"]))


def _write_grid(path, values, info):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    region = "{}/{}/{}/{}".format(info["x_min"], info["x_max"], info["y_min"], info["y_max"])
    increment = "{}/{}".format(info["x_inc"], info["y_inc"])
    with tempfile.TemporaryDirectory(prefix="product-mask-") as temp_dir:
        raw_path = Path(temp_dir) / "mask.bin"
        values.astype(np.float32, copy=False).tofile(raw_path)
        command = [
            "gmt",
            "xyz2grd",
            str(raw_path),
            f"-R{region}",
            f"-I{increment}",
            "-ZTLf",
            f"-G{path}",
        ]
        if info["registration"] == 1:
            command.append("-r")
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def _edge_widths(overlap_ok):
    invalid = ~overlap_ok
    profiles = {
        "top": np.mean(invalid, axis=1),
        "bottom": np.mean(invalid, axis=1)[::-1],
        "left": np.mean(invalid, axis=0),
        "right": np.mean(invalid, axis=0)[::-1],
    }
    widths = {}
    for edge, profile in profiles.items():
        width = 0
        for fraction in profile:
            if fraction < 0.9:
                break
            width += 1
        widths[edge] = width
    return widths


def build_mask(
    grid_path,
    slave_prm,
    slave_size_prm,
    corr_thresh,
    margin_px,
    out_mask,
    summary_json=None,
    master_prm=None,
    phase_grid=None,
    edge_noise="on",
):
    if not math.isfinite(corr_thresh) or not 0.0 <= corr_thresh < 1.0:
        raise ValueError("--corr-thresh must be in [0, 1)")
    if margin_px < 0:
        raise ValueError("--margin-px must be non-negative")
    if edge_noise not in ("on", "off"):
        raise ValueError("--edge-noise must be on or off")

    transform, subpixel_defaults = _read_prm(slave_prm, SHIFT_KEYS)
    slave_size, _ = _read_prm(slave_size_prm, ("num_rng_bins", "num_lines"))
    slave_nx = int(slave_size["num_rng_bins"])
    slave_ny = int(slave_size["num_lines"])
    if slave_nx < 1 or slave_ny < 1:
        raise ValueError(f"invalid slave dimensions in {slave_size_prm}")
    if master_prm:
        master_size, _ = _read_prm(master_prm, ("num_lines",))
        master_y_origin = int(master_size["num_lines"])
    else:
        master_y_origin = None

    grid = _grid_info(grid_path)
    corr = _read_grid(grid_path, grid)
    if grid["registration"] == 1:
        x = grid["x_min"] + (np.arange(grid["nx"], dtype=np.float64) + 0.5) * grid["x_inc"]
        y_grid = grid["y_max"] - (np.arange(grid["ny"], dtype=np.float64) + 0.5) * grid["y_inc"]
    else:
        x = grid["x_min"] + np.arange(grid["nx"], dtype=np.float64) * grid["x_inc"]
        y_grid = grid["y_max"] - np.arange(grid["ny"], dtype=np.float64) * grid["y_inc"]

    # GMT grids are north-up; resamp.c uses zero-based, top-down master line indices.
    y_origin = master_y_origin if master_y_origin is not None else grid["y_max"]
    y = y_origin - y_grid
    overlap_ok = np.empty(corr.shape, dtype=bool)
    for row, y_master in enumerate(y):
        xs = (
            x
            + transform["rshift"]
            + transform["sub_int_r"]
            + transform["stretch_r"] * x
            + transform["a_stretch_r"] * y_master
        )
        ys = (
            y_master
            + transform["ashift"]
            + transform["sub_int_a"]
            + transform["stretch_a"] * x
            + transform["a_stretch_a"] * y_master
        )
        overlap_ok[row, :] = (
            (xs >= margin_px)
            & (xs <= slave_nx - 1 - margin_px)
            & (ys >= margin_px)
            & (ys <= slave_ny - 1 - margin_px)
        )
    corr_ok = np.ones(corr.shape, dtype=bool) if corr_thresh == 0.0 else np.isfinite(corr) & (corr >= corr_thresh)
    edge_noise_widths = {edge: 0 for edge in ("top", "bottom", "left", "right")}
    edge_noise_extents = {edge: 0 for edge in ("top", "bottom", "left", "right")}
    edge_noise_pixels = np.zeros(corr.shape, dtype=bool)
    if edge_noise == "on":
        if not phase_grid:
            raise ValueError("--phase-grid is required when --edge-noise is on")
        phase = dinsar_qa._read_grid(phase_grid)
        if phase.shape != corr.shape:
            raise ValueError(f"phase grid shape {phase.shape} does not match coherence grid {corr.shape}")
        _, residue_mask, valid_loops = dinsar_qa._residue_metrics(phase, corr)
        edge_check = dinsar_qa._edge_strip_check(residue_mask, valid_loops, phase)
        edge_widths = edge_check["value"]["width_pixels"]
        edge_starts = edge_check["value"]["noisy_band_start_offset"]
        ny, nx = corr.shape
        for edge in edge_noise_widths:
            width = int(edge_widths[edge]["noisy_width"])
            start = edge_starts[edge]
            edge_noise_widths[edge] = width
            if width < 1 or start is None:
                continue
            extent = int(start) + width
            edge_noise_extents[edge] = extent
            if edge == "top":
                edge_noise_pixels[:extent, :] = True
            elif edge == "bottom":
                edge_noise_pixels[ny - extent :, :] = True
            elif edge == "left":
                edge_noise_pixels[:, :extent] = True
            else:
                edge_noise_pixels[:, nx - extent :] = True

    base_valid = overlap_ok & corr_ok
    newly_masked = base_valid & edge_noise_pixels
    valid = base_valid & ~edge_noise_pixels
    mask = np.where(valid, 1.0, np.nan).astype(np.float32)
    _write_grid(out_mask, mask, grid)

    total = int(mask.size)
    overlap_masked = int(np.count_nonzero(~overlap_ok))
    corr_masked = int(np.count_nonzero(overlap_ok & ~corr_ok))
    edge_noise_masked = int(np.count_nonzero(newly_masked))
    edge_noise_fraction = edge_noise_masked / total
    summary = {
        "grid": {
            **grid,
            "coordinate_convention": "x uses GMT grid coordinates; y is converted from north-up GMT coordinates to top-down resamp master line coordinates",
        },
        "slave_prm": str(Path(slave_prm)),
        "slave_size_prm": str(Path(slave_size_prm)),
        "master_prm": str(Path(master_prm)) if master_prm else None,
        "master_y_origin": y_origin,
        "slave_size": {"num_rng_bins": slave_nx, "num_lines": slave_ny},
        "corr_thresh": corr_thresh,
        "margin_px": margin_px,
        "subpixel_defaults": subpixel_defaults,
        "fractions": {
            "masked_fraction_by_overlap": overlap_masked / total,
            "masked_fraction_by_corr": corr_masked / total,
            "masked_fraction_total": float(np.count_nonzero(~valid) / total),
            "corr_below_threshold_fraction": float(np.count_nonzero(~corr_ok) / total),
        },
        "overlap_invalid_widths_px": _edge_widths(overlap_ok),
        "edge_noise_enabled": edge_noise == "on",
        "edge_noise_widths_px": edge_noise_widths,
        "edge_noise_masked_extent_px": edge_noise_extents,
        "edge_noise_mask_fraction": float(edge_noise_fraction),
        "edge_noise_masked_pixels": edge_noise_masked,
        "phase_grid": str(Path(phase_grid)) if phase_grid else None,
        "edge_width_invalid_line_fraction_min": 0.9,
        "output_mask": str(Path(out_mask)),
    }
    if summary_json:
        summary_path = Path(summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", required=True, help="input coherence GMT grid")
    parser.add_argument("--slave-prm", required=True, help="PRM containing resampling shifts")
    parser.add_argument("--slave-size-prm", required=True, help="original slave PRM with full dimensions")
    parser.add_argument("--master-prm", help="master PRM; supplies the top-down line origin for cropped grids")
    parser.add_argument("--corr-thresh", type=float, default=0.3)
    parser.add_argument("--margin-px", type=int, default=4)
    parser.add_argument("--phase-grid", help="filtered phase grid used by the QA edge-residual detector")
    parser.add_argument(
        "--edge-noise", choices=("on", "off"), default="on",
        help="mask edge residual-noise bands detected by dinsar_qa (default: on)",
    )
    parser.add_argument("--out-mask", required=True)
    parser.add_argument("--summary-json")
    args = parser.parse_args(argv)
    try:
        summary = build_mask(
            args.grid,
            args.slave_prm,
            args.slave_size_prm,
            args.corr_thresh,
            args.margin_px,
            args.out_mask,
            args.summary_json,
            args.master_prm,
            args.phase_grid,
            args.edge_noise,
        )
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"product_mask.py: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary["fractions"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
