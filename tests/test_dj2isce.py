"""Unit tests for gmtsar-master/gmtsar/csh/dj2isce.py

The module imports osgeo (GDAL) at the top level which is a heavy C
dependency not needed for testing pure-Python helpers. We mock it
before importing the module under test.
"""

import os
import sys
import types
import xml.etree.ElementTree as ET
import numpy as np
import pytest

# Stub osgeo.gdal so dj2isce can be imported without GDAL installed
_osgeo = types.ModuleType("osgeo")
_gdal = types.ModuleType("osgeo.gdal")
_osgeo.gdal = _gdal
_gdal.GA_ReadOnly = 0
_gdal.Open = lambda *a, **k: None
sys.modules.setdefault("osgeo", _osgeo)
sys.modules.setdefault("osgeo.gdal", _gdal)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'gmtsar-master', 'gmtsar', 'csh'))
import dj2isce


# ---------- find_text ----------

class TestFindText:

    def _root(self, xml_str):
        return ET.fromstring(xml_str)

    def test_finds_first_matching_tag(self):
        root = self._root("<r><PRF>1000</PRF></r>")
        assert dj2isce.find_text(root, ["PRF"]) == "1000"

    def test_falls_back_to_second_tag(self):
        root = self._root("<r><chirpBandwidth>50e6</chirpBandwidth></r>")
        assert dj2isce.find_text(root, ["rangeBandwidth", "chirpBandwidth"]) == "50e6"

    def test_returns_none_when_not_found(self):
        root = self._root("<r><other>1</other></r>")
        assert dj2isce.find_text(root, ["PRF"]) is None

    def test_strips_whitespace(self):
        root = self._root("<r><PRF>  1500  </PRF></r>")
        assert dj2isce.find_text(root, ["PRF"]) == "1500"

    def test_empty_text_returns_none(self):
        root = self._root("<r><PRF></PRF></r>")
        assert dj2isce.find_text(root, ["PRF"]) is None


# ---------- estimate_pulse_duration ----------

class TestEstimatePulseDuration:

    def test_with_chirp_duration(self):
        meta = {"chirpDuration": "2.5e-5"}
        assert dj2isce.estimate_pulse_duration(meta) == pytest.approx(2.5e-5)

    def test_fallback_to_bandwidth(self):
        meta = {"rangeBandwidth": "50000000"}
        assert dj2isce.estimate_pulse_duration(meta) == pytest.approx(1.0 / 50000000)

    def test_raises_when_no_info(self):
        with pytest.raises(RuntimeError, match="Cannot estimate pulseDuration"):
            dj2isce.estimate_pulse_duration({})


# ---------- estimate_chirp_slope ----------

class TestEstimateChirpSlope:

    def test_with_direct_chirp_slope(self):
        meta = {"chirpSlope": "5e12"}
        assert dj2isce.estimate_chirp_slope(meta, 1e-5) == pytest.approx(5e12)

    def test_small_chirp_slope_uses_fallback(self):
        meta = {"chirpSlope": "100", "chirp_rate": "-1", "rangeBandwidth": "50e6"}
        pulse_dur = 2e-5
        expected = -1 * 50e6 / pulse_dur
        assert dj2isce.estimate_chirp_slope(meta, pulse_dur) == pytest.approx(expected)

    def test_raises_when_no_info(self):
        with pytest.raises(RuntimeError, match="Cannot estimate chirpSlope"):
            dj2isce.estimate_chirp_slope({}, 1e-5)


# ---------- parse_orbit ----------

class TestParseOrbit:

    def _orbit_xml(self, n_vectors):
        parts = ["<root>"]
        for i in range(n_vectors):
            parts.append(f"""
            <stateVector>
                <time>2020-01-01T00:00:{i:02d}.000Z</time>
                <x>{1000 + i}</x><y>{2000 + i}</y><z>{3000 + i}</z>
                <vx>{10 + i}</vx><vy>{20 + i}</vy><vz>{30 + i}</vz>
            </stateVector>""")
        parts.append("</root>")
        return ET.fromstring("".join(parts))

    def test_parses_enough_vectors(self):
        root = self._orbit_xml(6)
        svs = dj2isce.parse_orbit(root)
        assert len(svs) == 6
        assert svs[0]["x"] == "1000"

    def test_raises_with_too_few_vectors(self):
        root = self._orbit_xml(3)
        with pytest.raises(RuntimeError, match="Orbit vectors < 5"):
            dj2isce.parse_orbit(root)

    def test_skips_incomplete_vectors(self):
        xml = """<root>
            <stateVector><time>T</time><x>1</x></stateVector>
            <stateVector><time>T</time><x>1</x><y>2</y><z>3</z><vx>4</vx><vy>5</vy><vz>6</vz></stateVector>
            <stateVector><time>T</time><x>1</x><y>2</y><z>3</z><vx>4</vx><vy>5</vy><vz>6</vz></stateVector>
            <stateVector><time>T</time><x>1</x><y>2</y><z>3</z><vx>4</vx><vy>5</vy><vz>6</vz></stateVector>
            <stateVector><time>T</time><x>1</x><y>2</y><z>3</z><vx>4</vx><vy>5</vy><vz>6</vz></stateVector>
            <stateVector><time>T</time><x>1</x><y>2</y><z>3</z><vx>4</vx><vy>5</vy><vz>6</vz></stateVector>
        </root>"""
        root = ET.fromstring(xml)
        svs = dj2isce.parse_orbit(root)
        assert len(svs) == 5

    def test_alternate_tag_names(self):
        xml = """<root>
            <OrbitStateVector><UTC>T1</UTC><PositionX>1</PositionX><PositionY>2</PositionY><PositionZ>3</PositionZ><VelocityX>4</VelocityX><VelocityY>5</VelocityY><VelocityZ>6</VelocityZ></OrbitStateVector>
            <OrbitStateVector><UTC>T2</UTC><PositionX>1</PositionX><PositionY>2</PositionY><PositionZ>3</PositionZ><VelocityX>4</VelocityX><VelocityY>5</VelocityY><VelocityZ>6</VelocityZ></OrbitStateVector>
            <OrbitStateVector><UTC>T3</UTC><PositionX>1</PositionX><PositionY>2</PositionY><PositionZ>3</PositionZ><VelocityX>4</VelocityX><VelocityY>5</VelocityY><VelocityZ>6</VelocityZ></OrbitStateVector>
            <OrbitStateVector><UTC>T4</UTC><PositionX>1</PositionX><PositionY>2</PositionY><PositionZ>3</PositionZ><VelocityX>4</VelocityX><VelocityY>5</VelocityY><VelocityZ>6</VelocityZ></OrbitStateVector>
            <OrbitStateVector><UTC>T5</UTC><PositionX>1</PositionX><PositionY>2</PositionY><PositionZ>3</PositionZ><VelocityX>4</VelocityX><VelocityY>5</VelocityY><VelocityZ>6</VelocityZ></OrbitStateVector>
        </root>"""
        root = ET.fromstring(xml)
        svs = dj2isce.parse_orbit(root)
        assert len(svs) == 5
