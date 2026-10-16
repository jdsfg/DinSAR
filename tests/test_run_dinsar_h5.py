#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_dinsar.sh 对新版电建一号 (XML + HDF5) 的接入 (把脚本当库 source 进来, 直接调用各函数)。

  1. 数据发现: 同名 .h5 优先于 8 bit 预览 .tiff, INPUT_FORMAT=h5; 没给 --coreg 时自动用 multiscale, 给了 --coreg legacy 仍尊重;
     旧版 XML+TIFF 不受影响 (INPUT_FORMAT=tiff, 配准仍默认 legacy); 主辅格式不一致时报错; 认得 Dianjian 卫星名;
  2. add_orbit_params_to_prm: 成对调用 calc_dop_orb (辅影像沿用主影像的 earth_radius), 不重复写 fd1;
  3. prepare_dem: GeoTIFF -> topo/dem.grd; 没有 DEM 返回 0, 转换失败返回 1;
  4. run_interferogram: 新版数据先转 DEM, 再生成 4 倍粗网格的 dem_ra (dec=auto), 生成失败直接中止;
     旧版数据仍是整景分辨率 (dec=1), 失败只警告, 且这一步不转 GeoTIFF DEM; TOPO_RA_DEC 可覆盖。
  5. 去倾斜与位移图: 新版数据默认去平面 (DETREND_ORDER=3) 且色标取 98 分位, 旧版数据仍是只减均值 + 最大绝对值;
     环境变量 DETREND_ORDER / DISP_STRETCH_PCT 明确指定的值优先; grdtrend -N3 真的去掉倾斜, -N1 不去;
     色标取分位数后, 个别离群值不再把位移图冲淡。

运行: python3 -m unittest tests.test_run_dinsar_h5 -v    (需要 bash; 第 3~5 部分还需要 gmt、gdalwarp、osgeo, 第 5 部分还需要 numpy、PIL)
"""

import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
HAVE_BASH = shutil.which("bash") is not None
HAVE_GEO = all(shutil.which(c) for c in ("gmt", "gdalwarp"))
try:
    from osgeo import gdal, osr
    HAVE_GDAL = True
except ImportError:
    HAVE_GDAL = False
try:
    import numpy as np
    from PIL import Image
    HAVE_IMAGING = True
except ImportError:
    HAVE_IMAGING = False

MASTER, SLAVE = "20260904", "20260912"
XML_H5 = ("<root><satellite_name>Dianjian-01</satellite_name>"
          "<acquisition_start_utc>%s</acquisition_start_utc></root>")
XML_OLD = "<root><SatelliteID>DJ1</SatelliteID><startTime>%s</startTime></root>"


def run_bash(snippet, path_prefix=None, env_extra=None, timeout=300):
    env = dict(os.environ)
    if path_prefix:
        env["PATH"] = path_prefix + os.pathsep + env.get("PATH", "")
    env.update(env_extra or {})
    return subprocess.run(["bash", "-c", snippet], cwd=REPO, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)


def write(path, text="", mode=None):
    with open(path, "w") as f:
        f.write(text)
    if mode:
        os.chmod(path, mode)


def write_stub(directory, name, body):
    write(os.path.join(directory, name), "#!/bin/sh\n" + body, 0o755)


def make_raw(project, kind):
    """raw/ 里放两景: kind = h5 (XML + h5 + 预览 tiff) / old (XML + tiff) / mixed (主影像有 h5, 辅影像只有 tiff)。"""
    raw = os.path.join(project, "raw")
    os.makedirs(raw)
    for date, iso in ((MASTER, "2026-09-04T04:23:44"), (SLAVE, "2026-09-12T04:23:43")):
        stem = "scene_" + date
        if kind == "old":
            write(os.path.join(raw, stem + ".xml"), XML_OLD % (iso + "Z"))
        else:
            write(os.path.join(raw, stem + ".XML"), XML_H5 % iso)
        write(os.path.join(raw, stem + ".tiff"))
        if kind == "h5" or (kind == "mixed" and date == MASTER):
            write(os.path.join(raw, stem + ".h5"))
    return raw


@unittest.skipUnless(HAVE_BASH, "需要 bash")
class Discovery(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="run_dinsar_h5_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _discover(self, kind, extra_args="", env_extra=None):
        project = os.path.join(self.tmp, kind + extra_args.replace(" ", "_").replace("-", "")
                               + "_".join(sorted((env_extra or {}).values())))
        os.makedirs(project)
        make_raw(project, kind)
        snippet = ('source ./run_dinsar.sh; LOG_FILE=/dev/null; %s; PROJECT_DIR="%s"; discover_data; rc=$?; '
                   'echo "rc=$rc fmt=$INPUT_FORMAT coreg=$COREG_MODE sat=$SATELLITE_TYPE '
                   'detrend=$DETREND_ORDER stretch=$DISP_STRETCH_PCT '
                   'm=$(basename "$MASTER_TIFF") s=$(basename "$SLAVE_TIFF") master=$MASTER_DATE slave=$SLAVE_DATE"'
                   % ("parse_args %s \"%s\"" % (extra_args, project) if extra_args else ":", project))
        r = run_bash(snippet, env_extra=env_extra)
        line = [ln for ln in r.stdout.splitlines() if ln.startswith("rc=")]
        self.assertTrue(line, r.stdout)
        return dict(kv.split("=", 1) for kv in line[-1].split()), r.stdout

    def test_h5_is_preferred_over_the_preview_tiff(self):
        res, out = self._discover("h5")
        self.assertEqual(res["rc"], "0", out)
        self.assertEqual((res["fmt"], res["sat"]), ("h5", "DJ1"))
        self.assertEqual((res["m"], res["s"]), ("scene_%s.h5" % MASTER, "scene_%s.h5" % SLAVE))
        self.assertEqual((res["master"], res["slave"]), (MASTER, SLAVE))

    def test_h5_defaults_to_multiscale_but_respects_an_explicit_coreg(self):
        self.assertEqual(self._discover("h5")[0]["coreg"], "multiscale")
        self.assertEqual(self._discover("h5", "--coreg legacy")[0]["coreg"], "legacy")
        self.assertEqual(self._discover("h5", "--coreg multiscale")[0]["coreg"], "multiscale")

    def test_old_xml_tiff_is_unchanged(self):
        res, out = self._discover("old")
        self.assertEqual(res["rc"], "0", out)
        self.assertEqual((res["fmt"], res["coreg"], res["sat"]), ("tiff", "legacy", "DJ1"))
        self.assertEqual(res["m"], "scene_%s.tiff" % MASTER)
        self.assertEqual((res["detrend"], res["stretch"]), ("1", "100"), "旧版数据的去趋势和色标不能变")

    def test_h5_defaults_to_plane_detrend_and_percentile_stretch(self):
        res, out = self._discover("h5")
        self.assertEqual((res["detrend"], res["stretch"]), ("3", "98"), out)
        self.assertIn("去趋势: 平面", out)

    def test_environment_variables_win_over_the_new_data_defaults(self):
        res, _ = self._discover("h5", env_extra={"DETREND_ORDER": "6", "DISP_STRETCH_PCT": "95"})
        self.assertEqual((res["detrend"], res["stretch"]), ("6", "95"))
        res, _ = self._discover("h5", env_extra={"DETREND_ORDER": "1", "DISP_STRETCH_PCT": "100"})
        self.assertEqual((res["detrend"], res["stretch"]), ("1", "100"), "明确要旧行为 (只减均值) 时不能被新版默认覆盖")
        res, _ = self._discover("old", env_extra={"DETREND_ORDER": "3"})
        self.assertEqual(res["detrend"], "3", "旧版数据也可以用环境变量改去趋势阶数")

    def test_mixed_formats_are_rejected(self):
        res, out = self._discover("mixed")
        self.assertEqual(res["rc"], "1", out)
        self.assertIn("格式不一致", out)

    def test_dianjian_satellite_name_is_recognised(self):
        r = run_bash('source ./run_dinsar.sh; normalize_satellite_type "Dianjian-01"; normalize_satellite_type "ICEYE-X"')
        self.assertEqual(r.stdout.split(), ["DJ1", "BC3"], r.stdout)


@unittest.skipUnless(HAVE_BASH, "需要 bash")
class OrbitParameters(unittest.TestCase):
    def test_pairwise_calc_dop_orb_and_no_duplicate_fd1(self):
        with tempfile.TemporaryDirectory(prefix="run_dinsar_orb_") as d:
            slc = os.path.join(d, "SLC")
            stubs = os.path.join(d, "stubs")
            os.makedirs(slc)
            os.makedirs(stubs)
            calls = os.path.join(d, "calls.log")
            write_stub(stubs, "calc_dop_orb",
                       'echo "$@" >> "%s"\n'
                       'RE=$3; [ "$3" = 0 ] && RE=6371000\n'
                       'printf "SC_vel = 7000.0\\nearth_radius = %%s\\nSC_height = 530000.0\\nfd1 = 99.0\\n" "$RE" > "$2"\n'
                       % calls)
            for date in (MASTER, SLAVE):
                write(os.path.join(slc, date + ".PRM"), "num_rng_bins = 8\nfd1 = 0.0\n")
            r = run_bash('source ./run_dinsar.sh; LOG_FILE=/dev/null; MASTER_DATE=%s; SLAVE_DATE=%s; cd "%s"; '
                         'add_orbit_params_to_prm; echo rc=$?' % (MASTER, SLAVE, slc), path_prefix=stubs)
            self.assertIn("rc=0", r.stdout, r.stdout)
            with open(calls) as f:
                self.assertEqual(f.read().splitlines(),
                                 ["%s.PRM %s.orb.tmp 0 0" % (MASTER, MASTER),
                                  "%s.PRM %s.orb.tmp 6371000 0" % (SLAVE, SLAVE)])
            for date in (MASTER, SLAVE):
                with open(os.path.join(slc, date + ".PRM")) as f:
                    prm = f.read().splitlines()
                self.assertIn("earth_radius = 6371000", prm)
                self.assertEqual(sum(1 for ln in prm if ln.startswith("SC_vel")), 1)
                self.assertEqual([ln for ln in prm if ln.startswith("fd1")], ["fd1 = 0.0"])
                self.assertFalse(os.path.exists(os.path.join(slc, date + ".orb.tmp")))


def make_dem_tiff(path):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, 11, 11, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((109.0, 0.001, 0, 40.011, 0, -0.001))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).Fill(1500.0)
    ds.FlushCache()


@unittest.skipUnless(HAVE_BASH and HAVE_GEO and HAVE_GDAL, "需要 bash、gmt、gdalwarp、osgeo")
class DemAndInterferogramStage(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="run_dinsar_dem_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _project(self, name, dem="tif"):
        project = os.path.join(self.tmp, name)
        for sub in ("SLC", "topo", "logs"):
            os.makedirs(os.path.join(project, sub))
        if dem == "tif":
            make_dem_tiff(os.path.join(project, "topo", "dem.tif"))
        elif dem == "corrupt":
            write(os.path.join(project, "topo", "dem.tif"), "not a tiff")
        return project

    def _stubs(self, name, make_dem_ra=True):
        """假的 fast_dem2topo_ra.sh (记录参数, 可选地生成 dem_ra.grd)、phasediff 和 SAT_baseline。"""
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        calls = os.path.join(d, "calls.log")
        write(os.path.join(d, "fast_dem2topo_ra.sh"),
              '#!/bin/bash\necho "dem2topo prm=$1 dem=$2 dec=$4" >> "%s"\n%s'
              % (calls, 'gmt grdmath -R0/4/0/4 -I1 1 = dem_ra.grd\n' if make_dem_ra else ""), 0o755)
        write_stub(d, "phasediff",
                   'echo "phasediff $*" >> "%s"\n'
                   'gmt grdmath -R0/4/0/4 -I1 1 = real.grd\ngmt grdmath -R0/4/0/4 -I1 1 = imag.grd\n' % calls)
        write_stub(d, "SAT_baseline", "exit 0\n")
        return d, calls

    def _interferogram(self, project, stubs, input_format, env_extra=None):
        snippet = ('source ./run_dinsar.sh; LOG_FILE="%s/logs/t.log"; PROJECT_DIR="%s"; SCRIPT_DIR="%s"; '
                   'MASTER_DATE=%s; SLAVE_DATE=%s; INPUT_FORMAT=%s; DEM_PATH=""; '
                   'run_interferogram; echo "rc=$?"' % (project, project, stubs, MASTER, SLAVE, input_format))
        return run_bash(snippet, path_prefix=stubs, env_extra=env_extra)

    def _calls(self, path):
        if not os.path.exists(path):
            return []
        with open(path) as f:
            return f.read().splitlines()

    def test_prepare_dem_converts_geotiff_and_reports_failure(self):
        project = self._project("prep")
        r = run_bash('source ./run_dinsar.sh; LOG_FILE="%s/logs/t.log"; PROJECT_DIR="%s"; DEM_PATH=""; '
                     'prepare_dem; echo "rc=$?"; gmt grdinfo -C "%s/topo/dem.grd" | cut -f2-5' % (project, project, project))
        self.assertIn("rc=0", r.stdout, r.stdout)
        w, e, s, n = [float(x) for x in r.stdout.strip().splitlines()[-1].split()]
        self.assertAlmostEqual(w, 109.0, delta=0.001)
        self.assertAlmostEqual(n, 40.011, delta=0.001)
        self.assertGreater(e, w)
        self.assertGreater(n, s)
        self.assertFalse(os.path.exists(os.path.join(project, "topo", "dem_tmp.nc")))

        none = self._project("prep_none", dem=None)
        r = run_bash('source ./run_dinsar.sh; LOG_FILE="%s/logs/t.log"; PROJECT_DIR="%s"; DEM_PATH=""; '
                     'prepare_dem; echo "rc=$?"' % (none, none))
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertFalse(os.path.exists(os.path.join(none, "topo", "dem.grd")))

        bad = self._project("prep_bad", dem="corrupt")
        r = run_bash('source ./run_dinsar.sh; LOG_FILE="%s/logs/t.log"; PROJECT_DIR="%s"; DEM_PATH=""; '
                     'prepare_dem; echo "rc=$?"' % (bad, bad))
        self.assertIn("rc=1", r.stdout, r.stdout)
        self.assertFalse(os.path.exists(os.path.join(bad, "topo", "dem.grd")))

    def test_h5_converts_dem_first_and_uses_a_coarse_topo_grid(self):
        project = self._project("h5_ok")
        stubs, calls = self._stubs("stubs_h5_ok")
        r = self._interferogram(project, stubs, "h5")
        self.assertIn("rc=0", r.stdout, r.stdout)
        lines = self._calls(calls)
        self.assertEqual(len(lines), 2, lines)
        self.assertEqual(lines[0], "dem2topo prm=%s.PRM dem=%s/topo/dem.grd dec=auto" % (MASTER, project))
        self.assertTrue(os.path.isfile(os.path.join(project, "topo", "dem.grd")), "DEM 应在去地形相位之前已转成 dem.grd")
        self.assertIn("-topo dem_ra.grd", lines[1])

    def test_h5_aborts_when_the_topo_grid_cannot_be_made(self):
        project = self._project("h5_fail")
        stubs, calls = self._stubs("stubs_h5_fail", make_dem_ra=False)
        r = self._interferogram(project, stubs, "h5")
        self.assertNotIn("rc=0", r.stdout)
        self.assertIn("雷达坐标 DEM (dem_ra.grd) 生成失败", r.stdout)
        self.assertEqual([ln for ln in self._calls(calls) if ln.startswith("phasediff")], [], "不应在没有地形的情况下继续")

    def test_h5_aborts_when_the_dem_cannot_be_converted(self):
        project = self._project("h5_bad_dem", dem="corrupt")
        stubs, calls = self._stubs("stubs_h5_bad_dem")
        r = self._interferogram(project, stubs, "h5")
        self.assertNotIn("rc=0", r.stdout)
        self.assertIn("DEM 转换失败", r.stdout)
        self.assertEqual(self._calls(calls), [])

    def test_h5_without_any_dem_still_runs_but_says_so(self):
        project = self._project("h5_nodem", dem=None)
        stubs, calls = self._stubs("stubs_h5_nodem")
        r = self._interferogram(project, stubs, "h5")
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertIn("含地形相位", r.stdout)
        lines = self._calls(calls)
        self.assertEqual(len(lines), 1, lines)
        self.assertNotIn("-topo", lines[0])

    def test_legacy_keeps_full_resolution_and_only_warns(self):
        project = self._project("old_ok", dem=None)
        os.rename(os.path.join(project, "topo"), os.path.join(project, "topo_unused"))
        os.makedirs(os.path.join(project, "topo"))
        subprocess.run(["gmt", "grdmath", "-R109/109.01/40/40.01", "-I0.001", "1500", "=",
                        os.path.join(project, "topo", "dem.grd")], check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
        stubs, calls = self._stubs("stubs_old_ok")
        r = self._interferogram(project, stubs, "tiff")
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertEqual(self._calls(calls)[0], "dem2topo prm=%s.PRM dem=%s/topo/dem.grd dec=1" % (MASTER, project))

        stubs2, calls2 = self._stubs("stubs_old_fail", make_dem_ra=False)
        project2 = self._project("old_fail", dem=None)
        shutil.copy(os.path.join(project, "topo", "dem.grd"), os.path.join(project2, "topo", "dem.grd"))
        r = self._interferogram(project2, stubs2, "tiff")
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertIn("干涉图将含地形相位", r.stdout)
        self.assertNotIn("-topo", self._calls(calls2)[1])

    def test_legacy_does_not_convert_a_geotiff_dem_at_this_stage(self):
        project = self._project("old_tif")          # topo/ 里只有 dem.tif, 旧流程在这一步不转它
        stubs, calls = self._stubs("stubs_old_tif")
        r = self._interferogram(project, stubs, "tiff")
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertFalse(os.path.exists(os.path.join(project, "topo", "dem.grd")))
        self.assertEqual([ln for ln in self._calls(calls) if ln.startswith("dem2topo")], [])

    def test_topo_ra_dec_environment_variable_overrides(self):
        project = self._project("h5_env")
        stubs, calls = self._stubs("stubs_h5_env")
        r = self._interferogram(project, stubs, "h5", env_extra={"TOPO_RA_DEC": "2"})
        self.assertIn("rc=0", r.stdout, r.stdout)
        self.assertTrue(self._calls(calls)[0].endswith("dec=2"), self._calls(calls))


@unittest.skipUnless(HAVE_BASH and HAVE_GEO and HAVE_GDAL and HAVE_IMAGING, "需要 bash、gmt、osgeo、numpy、PIL")
class DetrendAndDisplacementImage(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="run_dinsar_detrend_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.project = os.path.join(self.tmp, "p")
        self.intf = os.path.join(self.project, "intf", "pair")
        for sub in ("debug", "logs", "results"):
            os.makedirs(os.path.join(self.project, sub))
        os.makedirs(self.intf)
        self.prefix = ('source ./run_dinsar.sh; LOG_FILE="%s/logs/t.log"; PROJECT_DIR="%s"; INTF_DIR="%s"; '
                       'MASTER_DATE=%s; SLAVE_DATE=%s; SATELLITE_TYPE=DJ1; ' % (self.project, self.project, self.intf,
                                                                              MASTER, SLAVE))

    def _gmt_grid_from_array(self, name, array):
        """把 numpy 数组 (第 0 行 = 方位向第 0 行) 写成 GMT 网格, 与流水线里 -ZTLf 的约定一致。"""
        raw = os.path.join(self.tmp, name + ".f32")
        array.astype("<f4").tofile(raw)
        ny, nx = array.shape
        subprocess.run(["gmt", "xyz2grd", raw, "-R0/%d/0/%d" % (nx, ny), "-I1", "-r", "-ZTLf",
                        "-G" + os.path.join(self.intf, name)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def _range_of(self, name):
        out = subprocess.run(["gmt", "grdinfo", "-C", os.path.join(self.intf, name)], check=True,
                             stdout=subprocess.PIPE, universal_newlines=True).stdout.split()
        return float(out[5]), float(out[6])

    def test_detrend_n3_removes_a_plane_and_n1_does_not(self):
        yy, xx = np.mgrid[0:40, 0:60]
        unw = 0.5 * xx + 2.0 * yy + 3.0 + 0.01 * np.sin(3.0 * xx)
        for order, flat in (("3", True), ("1", False)):
            for stale in ("unwrap_detrended.grd", "unwrap.grd"):      # run_detrend 见到已有结果会跳过
                if os.path.exists(os.path.join(self.intf, stale)):
                    os.remove(os.path.join(self.intf, stale))
            self._gmt_grid_from_array("unwrap.grd", unw)
            r = run_bash(self.prefix + 'DETREND_ORDER=%s; run_detrend; echo "rc=$?"' % order)
            self.assertIn("rc=0", r.stdout, r.stdout)
            lo, hi = self._range_of("unwrap_detrended.grd")
            if flat:
                self.assertLess(hi - lo, 0.1, "平面应被去掉, 只剩 0.01 的小波动: %s" % r.stdout)
            else:
                self.assertGreater(hi - lo, 80.0, "只减均值时倾斜原样保留 (峰峰值约 109)")
                self.assertAlmostEqual((hi + lo) / 2.0, 0.0, delta=3.0)
            self.assertIn("grdtrend -N%s" % order, r.stdout)

    def _render(self, pct):
        """跑 generate_outputs (无 DEM 时只渲染 PNG), 返回位移图; pct=None 表示不设 DISP_STRETCH_PCT (旧默认)。"""
        phase = np.zeros((100, 100))
        disp = np.full((100, 100), 10.0)
        disp[:, 50:] = -10.0
        disp[10:12, :] = 20.0                     # 2% 的像元 |位移| = 20
        disp[0, 0], disp[0, 1], disp[1, 0] = 300.0, -300.0, 300.0     # 3 个离群值
        self._gmt_grid_from_array("phasefilt.grd", phase)
        self._gmt_grid_from_array("los_displacement.grd", disp)
        env = "DISP_STRETCH_PCT=%s; " % pct if pct is not None else ""
        r = run_bash(self.prefix + 'FILT_PHASE=phasefilt.grd; DEM_PATH=""; ' + env + 'generate_outputs; echo "rc=$?"')
        self.assertIn("rc=0", r.stdout, r.stdout)
        png = os.path.join(self.project, "results", "displacement.png")
        self.assertTrue(os.path.isfile(png), r.stdout)
        return np.array(Image.open(png).convert("RGB")), r.stdout

    def _log(self):
        with open(os.path.join(self.project, "logs", "t.log")) as f:
            return f.read()

    def test_percentile_colour_scale_keeps_a_few_outliers_from_washing_out_the_image(self):
        img, out = self._render(98)
        self.assertIn("Displacement colour scale: +/-20.0 mm", self._log())
        r, g, b = img[50, 10]                     # +10 mm
        self.assertEqual(r, 255)
        self.assertLess(g, 160, "色标 ±20 mm 时 +10 mm 应是明显的红色, 实际 RGB=%s" % ((r, g, b),))
        r, g, b = img[50, 80]                     # -10 mm
        self.assertEqual(b, 255)
        self.assertLess(g, 160)

    def test_default_colour_scale_is_still_the_maximum_absolute_value(self):
        img, _ = self._render(None)
        self.assertIn("Displacement colour scale: +/-300.0 mm", self._log())
        self.assertGreater(img[50, 10][1], 235, "色标 ±300 mm 时 +10 mm 几乎是白色 (旧行为)")


if __name__ == "__main__":
    unittest.main()
