"""Unit tests for gmtsar-master/gmtsar/csh/fetchOrbit.py"""

import os
import sys
import datetime
from unittest.mock import MagicMock, patch
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'gmtsar-master', 'gmtsar', 'csh'))
import fetchOrbit


# ---------- FileToTimeStamp ----------

class TestFileToTimeStamp:

    def test_standard_safe_name(self):
        name = "S1A_IW_SLC__1SDV_20200101T120000_20200101T120030_030000_037000_ABCD.SAFE"
        tstamp, sat, sstamp = fetchOrbit.FileToTimeStamp(name)
        assert sat == "S1A"
        assert isinstance(tstamp, datetime.datetime)
        assert tstamp.year == 2020

    def test_with_directory_prefix(self):
        name = "/data/orbits/S1B_IW_SLC__1SDV_20210315T060000_20210315T060030_025000_030000_1234.SAFE"
        tstamp, sat, sstamp = fetchOrbit.FileToTimeStamp(name)
        assert sat == "S1B"

    def test_short_date_format(self):
        name = "S1A_something_20200615_extra.SAFE"
        tstamp, sat, _ = fetchOrbit.FileToTimeStamp(name)
        assert tstamp == datetime.datetime(2020, 6, 15)


# ---------- MyHTMLParser ----------

class TestMyHTMLParser:

    def test_parses_eof_data(self):
        parser = fetchOrbit.MyHTMLParser("https://example.com")
        parser.feed('<html><body>S1A_OPER_AUX_POEORB_20200102.EOF</body></html>')
        assert len(parser.fileList) == 1
        assert parser.fileList[0][1] == "S1A_OPER_AUX_POEORB_20200102.EOF"

    def test_ignores_non_eof_data(self):
        parser = fetchOrbit.MyHTMLParser("https://example.com")
        parser.feed('<html><body>some random text</body></html>')
        assert len(parser.fileList) == 0

    def test_handles_quicklook_href(self):
        parser = fetchOrbit.MyHTMLParser("https://example.com")
        # URL with Products('Quicklook') gets the Quicklook part stripped
        html = '<a href="https://scihub.copernicus.eu/gnss/odata/v1/Products(\'abc\')/Products(\'Quicklook\')/rest">S1A_TEST.EOF</a>'
        parser.feed(html)
        assert len(parser.fileList) == 1
        assert "Quicklook" not in parser.fileList[0][0]


# ---------- fileToRange ----------

class TestFileToRange:

    def test_standard_orbit_filename(self):
        fname = "S1A_OPER_AUX_POEORB_OPOD_20200120T121212_V20200101T000000_20200102T235959.EOF"
        start, stop, mission = fetchOrbit.fileToRange(fname)
        assert mission == "S1A"
        assert start.year == 2020
        assert start.month == 1
        assert start.day == 1
        assert stop.day == 2

    def test_with_path(self):
        fname = "/data/S1B_OPER_AUX_RESORB_OPOD_20210301T121212_V20210301T000000_20210302T235959.EOF"
        start, stop, mission = fetchOrbit.fileToRange(fname)
        assert mission == "S1B"
        assert start.month == 3
        assert stop.month == 3


# ---------- download_file ----------

class TestDownloadFile:

    def test_successful_download(self, tmp_path):
        out = tmp_path / "orbit.eof"
        mock_response = MagicMock()
        mock_response.raise_for_status.return_value = None
        mock_response.iter_content.return_value = [b"orbit data"]

        session = MagicMock()
        session.get.return_value = mock_response

        result = fetchOrbit.download_file("https://example.com/orbit.eof", str(out), session)
        assert result is True
        assert out.read_bytes() == b"orbit data"

    def test_failed_download(self, tmp_path):
        out = tmp_path / "orbit.eof"
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = Exception("404")

        session = MagicMock()
        session.get.return_value = mock_response

        result = fetchOrbit.download_file("https://example.com/bad", str(out), session)
        assert result is False
