#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""低内存指标/门禁回归：直接加载模块，不导入整景合成测试，不读取 SLC。"""

import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "gmtsar-master/gmtsar/csh/ms_register.py"
SPEC = importlib.util.spec_from_file_location("ms_register_metrics_target", SCRIPT)
ms = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ms)


def measurements(residuals, model=None):
    residuals = np.asarray(residuals, dtype=float).reshape(-1, 2)
    rows = np.zeros((len(residuals), 7))
    rows[:, 0] = np.arange(len(rows)) * 3
    rows[:, 1] = np.arange(len(rows)) * 2
    rows[:, 2:4] = residuals
    if model is not None:
        rows[:, 2:4] += np.column_stack(model.predict(rows[:, 0], rows[:, 1]))
    rows[:, 4:6] = (0.8, 10)
    return rows


def check_rows(rows, n_candidates, model=None):
    model = model or ms.Affine.translation(0, 0)
    points = [(i, i) for i in range(n_candidates)]
    master, slave = object(), object()
    with mock.patch.object(ms, "make_lattice", return_value=points) as lattice, \
            mock.patch.object(ms, "measure_points", return_value=rows) as measure:
        result = ms.check_model(master, slave, model)
    lattice.assert_called_once_with(master, slave, model, 128, 8, (14, 14), tex=None, phase=(0.25, 0.75))
    measure.assert_called_once_with(master, slave, model, points, 128, 8, prep="sqrt", ovs=ms.CHECK_OVS)
    return result


class CheckMetrics(unittest.TestCase):
    def assert_strict_json(self, report):
        clean = ms._clean(report)
        self.assertEqual(json.loads(json.dumps(clean, allow_nan=False)), clean)
        # Exercise the real report writer without creating a file.
        with mock.patch("builtins.open", mock.mock_open()) as opened:
            ms.write_report("unused.json", report)
        text = "".join(call.args[0] for call in opened().write.call_args_list)
        self.assertNotIn("NaN", text)
        self.assertNotIn("Infinity", text)
        self.assertEqual(json.loads(text), clean)

    def test_bias_and_outliers_are_not_removed_from_rmse(self):
        model = ms.Affine([10, 0.125, -0.25], [-4, -0.5, 0.125])
        residuals = np.array([[1, -2], [1, -2], [1, -2], [5, -6]], dtype=float)
        rows = measurements(residuals, model)
        result = check_rows(rows, 10, model)
        # All legacy fields, including centered max_abs and the valid denominator, stay unchanged.
        expected = dict(n_points=10, n_valid=4, valid_frac=0.4, median_r=1.0, median_a=-2.0,
                        rms_r=0.0, rms_a=0.0, max_abs_r=4.0, max_abs_a=4.0, frac_within_0p3=0.0)
        for key, value in expected.items():
            self.assertEqual(result[key], value, key)
        np.testing.assert_array_equal(result["residuals"], np.column_stack([rows[:, :2], residuals, rows[:, 4]]))
        for axis, rmse in (("r", np.sqrt(7)), ("a", np.sqrt(12))):
            self.assertAlmostEqual(result["rmse_" + axis], rmse)
            self.assertEqual(result["robust_sigma_" + axis], result["rms_" + axis])
            self.assertEqual(result["mad_" + axis], 0.0)
        self.assertNotAlmostEqual(result["rmse_r"], float(np.std(residuals[:, 0])))
        self.assertIn("not RMSE", result["rms_definition"])
        self.assertIn("no mean subtraction", result["rmse_definition"])
        self.assertEqual(result["metric_n_points"], 4)
        self.assert_strict_json(result)

    def test_inclusive_radial_boundary_and_unequal_denominators(self):
        residuals = np.array([[0.3, 0], [-0.3, 0], [0, 0.3], [0, -0.3], [0.18, 0.24], [0, 0],
                              [np.nextafter(0.3, np.inf), 0], [0.25, 0.25]])
        rows = measurements(np.vstack([residuals, [0, 0], [0, 0]]))
        rows[0, 4] = 0.08  # ncc equality is valid
        rows[-2, 4] = np.nextafter(0.08, 0)  # low NCC is not valid
        rows[-1, 6] = 0.5  # edge equality is not valid
        result = check_rows(rows, 12)  # 12 candidates, 10 measured, 8 valid, 6 within
        self.assertEqual(result["n_valid"], 8)
        self.assertEqual(result["valid_frac"], 8 / 12)
        self.assertEqual(result["n_within_0p3"], 6)
        self.assertEqual(result["frac_within_0p3"], 6 / 8)
        self.assertEqual(result["frac_within_0p3_all_candidates"], 6 / 12)
        for i, axis in enumerate(("r", "a")):
            values = residuals[:, i]
            mad = float(np.median(np.abs(values - np.median(values))))
            self.assertEqual(result["mad_" + axis], mad)
            self.assertEqual(result["robust_sigma_" + axis], 1.4826 * mad)
            self.assertEqual(result["rms_" + axis], 1.4826 * mad)
            self.assertEqual(result["rmse_" + axis], float(np.sqrt(np.mean(values ** 2))))
        self.assert_strict_json(result)

    def test_no_candidates_no_measurements_and_all_rejected(self):
        rejected = measurements([[0, 0], [1, 1]])
        rejected[:, 4] = 0
        for rows, n_candidates in ((measurements([]), 0), (measurements([]), 10), (rejected, 10)):
            with self.subTest(n_candidates=n_candidates, n_measured=len(rows)):
                result = check_rows(rows, n_candidates)
                self.assertEqual(result["n_valid"], 0)
                self.assertEqual(result["valid_frac"], 0.0)
                self.assertEqual(result["n_within_0p3"], 0)
                self.assertEqual(result["frac_within_0p3_all_candidates"], 0.0 if n_candidates else None)
                for axis in ("r", "a"):
                    for name in ("mad_", "robust_sigma_", "rmse_"):
                        self.assertIsNone(result[name + axis])
                    self.assertNotIn("rms_" + axis, result)
                self.assertNotIn("frac_within_0p3", result)
                self.assert_strict_json(result)

    def test_low_valid_counts_preserve_legacy_absence_and_gate_failure(self):
        for n in (1, 2, 3, 4, 7):
            with self.subTest(n_valid=n):
                result = check_rows(measurements([[0.1, -0.05]] * n), 10)
                self.assertEqual(result["n_valid"], n)
                self.assertEqual(result["valid_frac"], n / 10)
                self.assertEqual(result["n_within_0p3"], n)
                self.assertEqual(result["frac_within_0p3_all_candidates"], n / 10)
                for axis, value in (("r", 0.1), ("a", 0.05)):
                    if n < 4:
                        self.assertIsNone(result["rmse_" + axis])
                        self.assertIsNone(result["robust_sigma_" + axis])
                        self.assertNotIn("rms_" + axis, result)
                        self.assertNotIn("median_" + axis, result)
                    else:
                        self.assertAlmostEqual(result["rmse_" + axis], value)
                        self.assertEqual(result["rms_" + axis], 0.0)
                report, _, master, points = gate_fixture()
                passed, reasons = ms.evaluate_gate(report, master, points, result, ms.GATE_DEFAULTS)
                self.assertFalse(passed)
                self.assertEqual(reasons, ["独立点复核有效点不足 (%s)" % n])
                self.assert_strict_json(result)

    def test_added_overflow_metrics_serialize_as_null(self):
        result = ms._residual_metrics(np.full(4, np.finfo(float).max), np.zeros(4))
        self.assertIsNone(result["rmse_r"])
        self.assertEqual(result["rmse_a"], 0.0)  # measured zero, not missing data
        self.assert_strict_json(result)


class FitMetrics(unittest.TestCase):
    def test_real_fit_reports_original_final_lstsq_inliers_only(self):
        y, x = np.mgrid[0:6, 0:6]
        x, y = x.ravel().astype(float), y.ravel().astype(float)
        rng = np.random.default_rng(42)
        dr = 3 + 0.03 * x - 0.02 * y + rng.normal(0, 0.04, len(x))
        da = -2 + 0.01 * x + 0.04 * y + rng.normal(0, 0.02, len(x))
        dr[[2, 17, 31]] += 5
        da[[2, 17, 31]] -= 7
        real_lstsq = np.linalg.lstsq
        solves = []

        def record_lstsq(a, b, **kwargs):
            result = real_lstsq(a, b, **kwargs)
            solves.append((a.copy(), b.copy(), result[0].copy()))
            return result

        # Spy on real tiny solves, without mocking the numerical result or inlier selection.
        with mock.patch.object(ms.np.linalg, "lstsq", side_effect=record_lstsq):
            fit = ms.fit_affine_robust(x, y, dr, da, 0.5, np.random.default_rng(0), n_iter=256)
        self.assertIsNotNone(fit)
        model, inliers, info = fit
        a, b, coef = solves[-1]
        residuals = b - a @ coef
        self.assertLess(len(residuals), len(x))
        self.assertGreaterEqual(len(residuals), 20)
        self.assertEqual(info["metric_n_points"], len(residuals))
        self.assertIn("before returned mask", info["metric_population"])
        self.assertIn("not median-centered MAD", info["sigma_definition"])
        self.assertEqual(info["rms_definition"], info["rmse_definition"])
        for i, axis in enumerate(("r", "a")):
            v = residuals[:, i]
            rmse = float(np.sqrt(np.mean(v ** 2)))
            mad = float(np.median(np.abs(v - np.median(v))))
            self.assertEqual(info["rms_" + axis], rmse)
            self.assertEqual(info["rmse_" + axis], rmse)
            self.assertEqual(info["mad_" + axis], mad)
            self.assertEqual(info["robust_sigma_" + axis], 1.4826 * mad)
            self.assertNotAlmostEqual(rmse, 1.4826 * mad, places=5)
        pr, pa = model.predict(x, y)
        np.testing.assert_array_equal(inliers, np.hypot(dr - pr, da - pa) <= info["t_in"])
        self.assertGreater(float(np.sqrt(np.mean((dr - pr) ** 2))), 10 * info["rmse_r"])
        rows = np.column_stack([x, y, dr, da, np.ones((len(x), 3))])
        stage = ms._stage_stats("S128", dict(W=128, R=8), list(zip(x, y)), rows,
                                np.ones(len(x), bool), inliers, info, 1.234)
        for key, value in info.items():
            self.assertEqual(stage[key], value, key)
        self.assertEqual(stage["n_inliers"], int(inliers.sum()))
        self.assertEqual(stage["time_s"], 1.23)
        json.dumps(ms._clean(stage), allow_nan=False)
        stage["accepted"] = True
        with mock.patch("builtins.open", mock.mock_open()) as opened:
            ms.write_text_log("unused.log", dict(stages=[stage]), "master", "slave")
        text = opened().write.call_args.args[0]
        self.assertIn("rmse=(", text)
        self.assertIn("robust_sigma=(", text)
        self.assertNotIn(" rms=(", text)

    def test_insufficient_or_degenerate_fit_keeps_none_and_null_stats(self):
        for n in (0, 5, 6):
            with self.subTest(n=n):
                v = np.zeros(n)
                self.assertIsNone(ms.fit_affine_robust(v, v, v, v, 0.5, np.random.default_rng(0), n_iter=64))
                stats = ms._stage_stats("S128", dict(W=128, R=8), [], measurements([]),
                                        np.zeros(0, bool), None, None, 0)
                for axis in ("r", "a"):
                    for key in ("rms_", "rmse_", "mad_", "robust_sigma_"):
                        self.assertIsNone(stats[key + axis])
                self.assertEqual(stats["metric_n_points"], 0)
                json.dumps(stats, allow_nan=False)


def gate_fixture(coverage=5):
    # Every existing numeric threshold is exactly on its passing boundary.
    stage = dict(name="S256", accepted=True, n_inliers=20, inlier_frac=0.30, rms_r=0.35, rms_a=0.35)
    check = dict(n_valid=8, valid_frac=0.30, median_r=0.12, median_a=-0.12, rms_r=0.35, rms_a=0.35)
    cells = [(5, 5), (15, 5), (25, 5), (5, 15), (15, 15)]
    points = np.zeros((coverage + 1, 8))
    points[:coverage, :2] = cells[:coverage]
    points[:coverage, 7] = 1
    points[-1, :2] = (25, 25)  # excluded outlier must not increase coverage
    return dict(stages=[stage]), check, SimpleNamespace(nx=30, ny=30), (dict(W=256), points)


class GateCompatibility(unittest.TestCase):
    def assert_gate_unchanged(self, report, check, master, points, passed):
        legacy = ms.evaluate_gate(copy.deepcopy(report), master, points, copy.deepcopy(check), ms.GATE_DEFAULTS)
        self.assertEqual(legacy[0], passed, legacy[1])
        self.assertEqual(len(legacy[1]), 0 if passed else 1)
        # Neither excellent nor terrible new diagnostics may override legacy decisions/reasons.
        for value in (0.0, 999.0):
            decorated, chk = copy.deepcopy(report), copy.deepcopy(check)
            extra = dict(rmse_r=value, rmse_a=value, robust_sigma_r=value, robust_sigma_a=value,
                         mad_r=value, mad_a=value, n_within_0p3=0, frac_within_0p3_all_candidates=0.0,
                         rms_definition="metadata only")
            decorated["stages"][-1].update(extra)
            chk.update(extra)
            self.assertEqual(ms.evaluate_gate(decorated, master, points, chk, ms.GATE_DEFAULTS), legacy)

    def test_all_threshold_equalities_pass(self):
        self.assertEqual(ms.GATE_DEFAULTS, dict(min_stage="S256", min_inliers=20, min_inlier_frac=0.30,
                                              max_rms=0.35, min_coverage=5, check_median=0.12,
                                              check_rms=0.35, check_valid_frac=0.30))
        report, check, master, points = gate_fixture()
        self.assert_gate_unchanged(report, check, master, points, True)
        ms.evaluate_gate(report, master, points, check, ms.GATE_DEFAULTS)
        self.assertEqual(report["coverage_cells"], 5)

    def test_each_threshold_just_outside_fails(self):
        cases = [("stage", "name", "S512"), ("stage", "n_inliers", 19),
                 ("stage", "inlier_frac", np.nextafter(0.30, 0)), ("check", "n_valid", 7),
                 ("check", "valid_frac", np.nextafter(0.30, 0)), ("coverage", None, 4)]
        for axis in ("r", "a"):
            cases += [("stage", "rms_" + axis, np.nextafter(0.35, np.inf)),
                      ("check", "rms_" + axis, np.nextafter(0.35, np.inf)),
                      ("check", "median_" + axis, np.nextafter(0.12, np.inf)),
                      ("check", "median_" + axis, np.nextafter(-0.12, -np.inf))]
        for target, key, value in cases:
            with self.subTest(target=target, key=key, value=value):
                report, check, master, points = gate_fixture(4 if target == "coverage" else 5)
                if target != "coverage":
                    (report["stages"][-1] if target == "stage" else check)[key] = value
                self.assert_gate_unchanged(report, check, master, points, False)

    def test_only_accepted_stages_count_and_optional_checks_remain_optional(self):
        report, check, master, points = gate_fixture()
        report["stages"].append(dict(name="S64", accepted=False))
        self.assert_gate_unchanged(report, check, master, points, True)
        self.assertEqual(ms.evaluate_gate(report, master, None, None, ms.GATE_DEFAULTS), (True, []))
        report["stages"][0]["accepted"] = False
        self.assertEqual(ms.evaluate_gate(report, master, points, check, ms.GATE_DEFAULTS),
                         (False, ["没有任何一级通过检验"]))

    def test_estimate_exit_status_and_check_log_keep_legacy_semantics(self):
        for passed in (True, False):
            with self.subTest(passed=passed):
                report, check, master, _ = gate_fixture()
                report["time_total_s"] = 0.0
                check.update(n_points=10, rmse_r=999.0, rmse_a=999.0)
                if not passed:
                    check["rms_a"] = 0.36
                    check.update(rmse_r=0.0, rmse_a=0.0)
                args = SimpleNamespace(master="master", slave="slave", scales=None, ovs=None, prep="sqrt",
                                       global_target=1500, seed=0, output="output", write_on_fail=False,
                                       points=None, report="report", log="log")
                with mock.patch.object(ms, "SLC", return_value=master), \
                        mock.patch.object(ms, "estimate", return_value=(report, ms.Affine.translation(0, 0), None)), \
                        mock.patch.object(ms, "check_model", return_value=check), \
                        mock.patch.object(ms, "write_prm") as write_prm, \
                        mock.patch.object(ms, "write_report"), mock.patch.object(ms, "write_text_log"), \
                        mock.patch.object(ms, "log") as log:
                    # Absolute output and PRM paths avoid touching any real inputs.
                    args.output = "/unused/output.PRM"
                    master.prm_path = "/unused/slave.PRM"
                    rc, result = ms.run_estimate(args)
                self.assertEqual(rc, ms.EXIT_OK if passed else ms.EXIT_GATE)
                self.assertEqual(result["status"], "passed" if passed else "failed")
                self.assertEqual(result["gate"]["passed"], passed)
                self.assertEqual(write_prm.call_count, int(passed))
                check_log = next(call.args[0] for call in log.call_args_list if "模型预测" in call.args[0])
                self.assertIn("robust_sigma=(", check_log)
                self.assertNotIn(" rms=(", check_log)


if __name__ == "__main__":
    unittest.main()
