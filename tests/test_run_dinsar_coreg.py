#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_dinsar.sh 的 --coreg 开关 (把脚本当库 source 进来, 直接调用 parse_args / run_alignment)。

  1. 命令行: 默认 legacy; multiscale / ms / legacy / xcorr 可识别; 无效值报错退出; --help 里有说明;
  2. 不加 --coreg: 仍然走 legacy (调用 xcorr 和 resamp, 不碰 ms_register), 行为与加开关之前一致;
  3. --coreg multiscale (真实 GMTSAR resamp, 找不到则跳过): 产出与 legacy 同名的 ${从影像}RS.PRM / RS.SLC,
     RS.PRM 的 SLC_file 是从影像原来的文件名 (run_interferogram 会把 RS.SLC 链接成这个名字),
     输入的 PRM/SLC 不被改动, 偏移模型准确;
  4. 配准没通过 (不相关的影像): 流程中止、不退回 legacy、不留下 RS 文件。

运行: python3 -m unittest tests.test_run_dinsar_coreg -v    (需要 bash、numpy、scipy)
"""

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_ms_register as base  # noqa: E402

ms = base.ms
REPO = os.path.abspath(os.path.join(HERE, ".."))
SCRIPT = os.path.join(REPO, "run_dinsar.sh")
RESAMP = shutil.which("resamp") or os.path.expanduser("~/gmtsar/gmtsar-master/bin/resamp")
HAVE_RESAMP = os.path.isfile(RESAMP) and os.access(RESAMP, os.X_OK)
MASTER, SLAVE = "20300101", "20300113"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_bash(snippet, path_prefix=None, timeout=900):
    env = dict(os.environ)
    if path_prefix:
        env["PATH"] = path_prefix + os.pathsep + env.get("PATH", "")
    return subprocess.run(["bash", "-c", snippet], cwd=REPO, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)


def align_snippet(project, mode):
    return ('source ./run_dinsar.sh; LOG_FILE="%s/align.log"; PROJECT_DIR="%s"; MASTER_DATE=%s; SLAVE_DATE=%s; '
            'PROCESS_MODE=FAST; SATELLITE_TYPE=DJ1; COREG_MODE=%s; run_alignment'
            % (project, project, MASTER, SLAVE, mode))


@unittest.skipUnless(shutil.which("bash"), "需要 bash")
class CommandLine(unittest.TestCase):
    def test_default_and_parsing(self):
        with tempfile.TemporaryDirectory() as d:
            r = run_bash('source ./run_dinsar.sh; LOG_FILE=/dev/null; echo "default=$COREG_MODE"; '
                         'parse_args --coreg MultiScale "%s"; echo "ms=$COREG_MODE"; '
                         'parse_args --coreg xcorr "%s"; echo "legacy=$COREG_MODE"; '
                         'parse_args --coreg ms "%s"; echo "short=$COREG_MODE"' % (d, d, d))
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertEqual([ln for ln in r.stdout.splitlines() if "=" in ln],
                         ["default=legacy", "ms=multiscale", "legacy=legacy", "short=multiscale"])

    def test_invalid_value_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            r = run_bash('source ./run_dinsar.sh; LOG_FILE=/dev/null; parse_args --coreg bogus "%s"; echo reached' % d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("无效 --coreg", r.stdout)
        self.assertNotIn("reached", r.stdout)

    def test_help_mentions_flag(self):
        r = subprocess.run(["bash", SCRIPT, "--help"], cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, timeout=60)
        self.assertEqual(r.returncode, 0)
        self.assertIn("--coreg legacy|multiscale", r.stdout)

    def test_multiscale_binary_preflight_and_python_override(self):
        with tempfile.TemporaryDirectory() as script_dir:
            r = run_bash(
                'source ./run_dinsar.sh; LOG_FILE=/dev/null; SCRIPT_DIR="%s"; PATH=/nonexistent; '
                'COREG_MODE=multiscale; MS_REGISTER_BIN=/missing/ms_register; MS_REGISTER_ARGS=""; '
                'error_exit() { echo "preflight_error: $*"; return 42; }; '
                'preflight_ms_register; echo "missing_rc=$?"; '
                'MS_REGISTER_ARGS="--engine python"; preflight_ms_register; echo "python_rc=$?"'
                % script_dir
            )
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertIn("missing_rc=42", r.stdout)
        self.assertIn("make -C gmtsar-master/gmtsar ms_register", r.stdout)
        self.assertIn('MS_REGISTER_ARGS="--engine python"', r.stdout)
        self.assertIn("python_rc=0", r.stdout)

    def test_sourcing_does_not_run_main(self):
        r = run_bash("source ./run_dinsar.sh; echo sourced")
        self.assertEqual(r.stdout.strip().splitlines()[-1], "sourced")
        self.assertNotIn("开始处理", r.stdout)


@unittest.skipUnless(shutil.which("bash") and base.HAVE_SCIPY, "需要 bash 与 scipy")
class Alignment(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="run_dinsar_coreg_")
        ms._QUIET = True
        cls.pair = os.path.join(cls.tmp, "pair")
        os.makedirs(cls.pair)
        base.make_pair(cls.pair, 0.9)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _project(self, name, unrelated=False):
        """SLC/ 目录里放 ${主}.PRM/.SLC 和 ${从}.PRM/.SLC (run_alignment 约定的布局); 返回项目目录。"""
        project = os.path.join(self.tmp, name)
        slc = os.path.join(project, "SLC")
        os.makedirs(slc)
        shutil.copy(os.path.join(self.pair, "m.SLC"), os.path.join(slc, MASTER + ".SLC"))
        if unrelated:
            rng = np.random.default_rng(5)
            noise = rng.integers(-3000, 3000, size=2 * base.NX * base.NY, dtype=np.int16)
            noise.tofile(os.path.join(slc, SLAVE + ".SLC"))
        else:
            shutil.copy(os.path.join(self.pair, "s.SLC"), os.path.join(slc, SLAVE + ".SLC"))
        base.write_prm(os.path.join(slc, MASTER + ".PRM"), MASTER + ".SLC")
        base.write_prm(os.path.join(slc, SLAVE + ".PRM"), SLAVE + ".SLC")
        return project

    def _stubs(self, name, resamp=True):
        """假的 xcorr (和 resamp): 只记录调用, resamp 桩还会创建输出文件。
        multiscale 测试要用真实 resamp, 所以那里 resamp=False, 免得桩排在 PATH 前面被 ms_register 当成 resamp。"""
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        log = os.path.join(d, "calls.log")
        scripts = [("xcorr", 'echo "xcorr $*" >> "%s"\n' % log)]
        if resamp:
            scripts.append(("resamp", 'echo "resamp $*" >> "%s"\ntouch "$3" "$4"\n' % log))
        for exe, body in scripts:
            p = os.path.join(d, exe)
            with open(p, "w") as f:
                f.write("#!/bin/sh\n" + body)
            os.chmod(p, 0o755)
        return d, log

    def test_default_path_is_still_legacy(self):
        project = self._project("p_legacy")
        stub_dir, calls = self._stubs("stubs_legacy")
        r = run_bash(align_snippet(project, "legacy"), path_prefix=stub_dir)
        self.assertEqual(r.returncode, 0, r.stdout)
        with open(calls) as f:
            lines = f.read().splitlines()
        self.assertEqual(len(lines), 2, lines)
        self.assertTrue(lines[0].startswith("xcorr %s.PRM %s.PRM" % (MASTER, SLAVE)), lines[0])
        self.assertEqual(lines[1], "resamp %s.PRM %s.PRM %sRS.PRM %sRS.SLC 4" % (MASTER, SLAVE, SLAVE, SLAVE))
        self.assertFalse(os.path.exists(os.path.join(project, "SLC", "coreg_report.json")))

    @unittest.skipUnless(HAVE_RESAMP, "需要 GMTSAR resamp")
    def test_multiscale_produces_legacy_named_outputs(self):
        project = self._project("p_ms")
        slc = os.path.join(project, "SLC")
        names = [MASTER + ".SLC", SLAVE + ".SLC", MASTER + ".PRM", SLAVE + ".PRM"]
        before = {n: md5(os.path.join(slc, n)) for n in names}
        stub_dir, calls = self._stubs("stubs_ms", resamp=False)    # 若误走 legacy, xcorr 桩会留下记录
        r = run_bash(align_snippet(project, "multiscale"), path_prefix=stub_dir + os.pathsep + os.path.dirname(RESAMP))
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertFalse(os.path.exists(calls), "multiscale 模式不应调用 xcorr")
        for n in (SLAVE + "RS.SLC", SLAVE + "RS.PRM", SLAVE + "_est.PRM", "coreg_report.json"):
            self.assertTrue(os.path.isfile(os.path.join(slc, n)), n)
        self.assertEqual(before, {n: md5(os.path.join(slc, n)) for n in names}, "输入文件被改动")
        self.assertEqual(os.path.getsize(os.path.join(slc, SLAVE + "RS.SLC")), 4 * base.NX * base.NY)
        prm = ms.read_prm(os.path.join(slc, SLAVE + "RS.PRM"))
        self.assertEqual(prm["SLC_file"], SLAVE + ".SLC")
        self.assertLess(base.model_error(base.model_from_prm(os.path.join(slc, SLAVE + "_est.PRM")), base.TRUTH), 0.1)
        with open(os.path.join(slc, "coreg_report.json")) as f:
            rep = json.load(f)
        self.assertEqual(rep["status"], "passed")
        self.assertTrue(rep["verify"]["gate"]["passed"])
        # run_interferogram 的做法: 在干涉目录里把 RS.SLC 链接成 ${从}.SLC, PRM 的 SLC_file 随即能解析到重采样后的数据
        intf = os.path.join(project, "intf")
        os.makedirs(intf)
        os.symlink(os.path.join(slc, SLAVE + "RS.SLC"), os.path.join(intf, SLAVE + ".SLC"))
        self.assertEqual(os.path.getsize(os.path.join(intf, prm["SLC_file"])), 4 * base.NX * base.NY)

    def test_multiscale_gate_failure_aborts_without_fallback(self):
        project = self._project("p_bad", unrelated=True)
        slc = os.path.join(project, "SLC")
        stub_dir, calls = self._stubs("stubs_bad")                 # 失败路径里 resamp 不会被调用, 桩只用来抓"退回 legacy"
        r = run_bash(align_snippet(project, "multiscale"), path_prefix=stub_dir)
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("多尺度配准未通过 (退出码 3)", r.stdout)      # 3 = 质量门限未通过, 排除环境问题冒充失败
        self.assertFalse(os.path.exists(calls), "未通过时不应退回 legacy (xcorr/resamp)")
        self.assertFalse(os.path.exists(os.path.join(slc, SLAVE + "RS.SLC")))
        self.assertFalse(os.path.exists(os.path.join(slc, SLAVE + "RS.PRM")))
        with open(os.path.join(slc, "coreg_report.json")) as f:
            self.assertEqual(json.load(f)["status"], "failed")


if __name__ == "__main__":
    unittest.main()
