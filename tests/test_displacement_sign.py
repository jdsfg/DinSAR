import math
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
from osgeo import gdal


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "run_dinsar.sh"
GEOCODE = REPO / "gmtsar-master/gmtsar/csh/geocode.csh"
LAMBDA = 0.0311


@unittest.skipUnless(all(shutil.which(c) for c in ("bash", "gmt")), "需要 bash 和 GMT")
class DisplacementSign(unittest.TestCase):
    """钉死相位到 LOS 位移的符号, 防止被误"更正"成相反方向。

    约定来源: GMTSAR 自带 geocode.csh 注释 "constant is negative to make
    LOS = -1 * range change", 常数为 1000/(4*pi)。正相位 -> 负位移(LOS 减小)。
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="disp_sign_")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def pipeline_expression(self):
        text = SCRIPT.read_text()
        match = re.search(r'gmt grdmath "\$input_phase"[^\n]*los_displacement\.grd', text)
        self.assertIsNotNone(match, "run_dinsar.sh 中未找到位移换算表达式")
        return match.group(0)

    def convert(self, phase):
        source = self.work / "unwrap.grd"
        output = self.work / "out.grd"
        subprocess.run(["gmt", "grdmath", "-R0/4/0/4", "-I1", str(phase), "=", str(source)],
                       check=True, capture_output=True, timeout=30)
        command = (self.pipeline_expression()
                   .replace('"$input_phase"', str(source))
                   .replace("$WAVELENGTH", repr(LAMBDA))
                   .replace("los_displacement.grd", str(output)))
        subprocess.run(["bash", "-c", command], check=True, capture_output=True, timeout=30)
        return gdal.Open(str(output)).ReadAsArray()

    def test_positive_phase_gives_negative_los(self):
        values = self.convert(2.0)
        expected = -2.0 * LAMBDA * 1000.0 / (4.0 * math.pi)
        np.testing.assert_allclose(values, expected, rtol=0, atol=1e-5)
        self.assertLess(float(values[0, 0]), 0.0)

    def test_conversion_is_linear_in_phase(self):
        np.testing.assert_allclose(self.convert(-3.0), -3.0 * self.convert(1.0), rtol=0, atol=1e-5)

    def test_constant_matches_gmtsar_geocode_convention(self):
        doc = GEOCODE.read_text()
        self.assertIn("LOS = -1 * range change", doc)
        constant = abs(float(re.search(r"-\d+\.\d+ MUL", doc).group(0).split()[0]))
        self.assertAlmostEqual(constant, 1000.0 / (4.0 * math.pi), delta=0.01)
        self.assertAlmostEqual(self.convert(1.0)[0, 0], -constant * LAMBDA, delta=0.02 * LAMBDA)


if __name__ == "__main__":
    unittest.main()