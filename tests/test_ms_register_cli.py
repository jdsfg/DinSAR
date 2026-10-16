#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ms_register register 子命令 (估计 -> resamp 重采样 -> 独立复核) 的集成测试。

  1. 端到端 (真实 GMTSAR resamp, 找不到则跳过): 输出齐全、模型恢复准确、输入文件逐字节不变、
     配准后的 PRM 指向配准后的 SLC、对输出独立 verify 通过、临时目录已清理;
  2. resamp 失败: 不留下半成品, 报告写成 failed;
  3. resamp 抛出意外异常: 磁盘上不能残留估计阶段写下的 status=passed 报告;
  4. 输出路径与输入相同: 拒绝执行, 输入不被改动;
  5. 质量门限未通过: 清除上一次遗留的同名输出, 不产生新输出。

运行: python3 -m unittest discover -s tests -v    (需要 numpy 与 scipy)
"""

import hashlib
import json
import os
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_ms_register as base  # noqa: E402

ms = base.ms
RESAMP = shutil.which("resamp") or os.path.expanduser("~/gmtsar/gmtsar-master/bin/resamp")
HAVE_RESAMP = os.path.isfile(RESAMP) and os.access(RESAMP, os.X_OK)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@unittest.skipUnless(base.HAVE_SCIPY, "合成数据生成需要 scipy")
class RegisterCommand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ms_register_cli_")
        ms._QUIET = True
        cls.good = os.path.join(cls.tmp, "good")
        os.makedirs(cls.good)
        base.make_pair(cls.good, 0.9)
        cls.bad = os.path.join(cls.tmp, "unrelated")
        os.makedirs(cls.bad)
        base.make_pair(cls.bad, 0.9, unrelated=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _exe(self, name, body):
        p = os.path.join(self.tmp, name)
        with open(p, "w") as f:
            f.write(body)
        os.chmod(p, 0o755)
        return p

    def _register(self, d, out, *extra):
        return ms.main(["register", os.path.join(d, "m.PRM"), os.path.join(d, "s.PRM"), "--outdir", out, "-q"]
                       + list(extra))

    @staticmethod
    def _report(out):
        with open(os.path.join(out, "coreg_report.json")) as f:
            return json.load(f)

    def test_native_verify_uses_requested_options(self):
        check = {"n_valid": 8}
        gate = {"passed": True, "reasons": [], "thresholds": {}}

        def fake_run(cmd, **kwargs):
            report_path = cmd[cmd.index("--report") + 1]
            with open(report_path, "w") as f:
                json.dump({"check": check, "gate": gate}, f)
            return SimpleNamespace(returncode=ms.EXIT_OK, stdout="")

        with mock.patch.object(ms, "_find_c_engine", return_value="/fake/ms_register"), \
                mock.patch.object(ms.subprocess, "run", side_effect=fake_run) as run:
            result = ms.verify_registered_c("/master.PRM", "/out.PRM", "/out.SLC",
                                            grid=(10, 12), prep="log", threads=3)

        cmd = run.call_args[0][0]
        self.assertEqual(cmd[:4], ["/fake/ms_register", "verify", "/master.PRM", "/out.PRM"])
        for option, value in (("--slave-slc", "/out.SLC"), ("--grid", "10x12"),
                              ("--prep", "log"), ("--threads", "3")):
            self.assertEqual(cmd[cmd.index(option) + 1], value)
        self.assertEqual(result, (check, gate))

    @unittest.skipUnless(HAVE_RESAMP, "需要 GMTSAR resamp")
    def test_end_to_end_with_real_resamp(self):
        d, out = self.good, os.path.join(self.tmp, "out_ok")
        before = {n: md5(os.path.join(d, n)) for n in ("m.SLC", "s.SLC", "m.PRM", "s.PRM")}
        self.assertEqual(self._register(d, out, "--resamp", RESAMP), ms.EXIT_OK)
        for n in ("s_reg.SLC", "s_reg.PRM", "s_est.PRM", "coreg_report.json", "coreg_offsets.log", "coreg_points.dat"):
            self.assertTrue(os.path.isfile(os.path.join(out, n)), n)
        self.assertEqual(before, {n: md5(os.path.join(d, n)) for n in before}, "输入文件被改动")
        self.assertEqual(os.path.getsize(os.path.join(out, "s_reg.SLC")), 4 * base.NX * base.NY)
        prm = ms.read_prm(os.path.join(out, "s_reg.PRM"))
        self.assertEqual(prm["SLC_file"], "s_reg.SLC")
        self.assertLess(base.model_error(base.model_from_prm(os.path.join(out, "s_est.PRM")), base.TRUTH), 0.1)
        rep = self._report(out)
        self.assertEqual(rep["status"], "passed")
        self.assertEqual(rep["mode"], "register")
        self.assertEqual(rep["resample"]["interp"], 4)
        self.assertTrue(rep["verify"]["gate"]["passed"])
        self.assertEqual([n for n in os.listdir(out) if n.startswith(".ms_register_")], [], "临时目录没有清理")
        # 对输出再做一次独立复核
        self.assertEqual(ms.main(["verify", os.path.join(d, "m.PRM"), os.path.join(out, "s_reg.PRM"), "-q"]), ms.EXIT_OK)

    def test_resamp_failure_leaves_no_product(self):
        out = os.path.join(self.tmp, "out_fail")
        exe = self._exe("fail_resamp.sh", "#!/bin/sh\necho boom >&2\nexit 1\n")
        self.assertEqual(self._register(self.good, out, "--resamp", exe), ms.EXIT_INTERNAL)
        self.assertFalse(os.path.exists(os.path.join(out, "s_reg.SLC")))
        self.assertFalse(os.path.exists(os.path.join(out, "s_reg.PRM")))
        rep = self._report(out)
        self.assertEqual(rep["status"], "failed")
        self.assertIn("重采样失败", rep["fail_reason"])
        self.assertEqual([n for n in os.listdir(out) if n.startswith(".ms_register_")], [])

    def test_unexpected_error_rewrites_stale_passed_report(self):
        """解释器不存在 -> subprocess 抛 OSError (不是 ResampError); 估计阶段写下的 passed 报告必须被改写"""
        out = os.path.join(self.tmp, "out_oserror")
        exe = self._exe("bad_interp.sh", "#!/nonexistent/interpreter\n")
        self.assertEqual(self._register(self.good, out, "--resamp", exe), ms.EXIT_INTERNAL)
        rep = self._report(out)
        self.assertEqual(rep["status"], "failed")
        self.assertIn("意外错误", rep["fail_reason"])
        self.assertFalse(os.path.exists(os.path.join(out, "s_reg.SLC")))

    def test_refuses_to_overwrite_inputs(self):
        d = self.good
        before = md5(os.path.join(d, "s.SLC"))
        rc = ms.main(["register", os.path.join(d, "m.PRM"), os.path.join(d, "s.PRM"), "--outdir",
                      os.path.join(self.tmp, "out_x"), "--out-slc", os.path.join(d, "s.SLC"), "--resamp", "/bin/true", "-q"])
        self.assertEqual(rc, ms.EXIT_INPUT)
        self.assertEqual(before, md5(os.path.join(d, "s.SLC")))

    def test_gate_failure_clears_stale_outputs(self):
        out = os.path.join(self.tmp, "out_gate")
        os.makedirs(out)
        for n in ("s_reg.SLC", "s_reg.PRM", "s_est.PRM"):      # 上一次运行遗留的结果
            with open(os.path.join(out, n), "w") as f:
                f.write("stale")
        exe = self._exe("never_called.sh", "#!/bin/sh\nexit 0\n")
        self.assertEqual(self._register(self.bad, out, "--resamp", exe), ms.EXIT_GATE)
        for n in ("s_reg.SLC", "s_reg.PRM", "s_est.PRM"):
            self.assertFalse(os.path.exists(os.path.join(out, n)), n + " 没有被清除")
        self.assertEqual(self._report(out)["status"], "failed")


if __name__ == "__main__":
    unittest.main()
