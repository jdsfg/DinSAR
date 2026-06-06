"""Unit tests for gmtsar-master/gmtsar/python/utils/gmtsar_lib.py"""

import os
import sys
import tempfile
import pytest

# Add module path so we can import gmtsar_lib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'gmtsar-master', 'gmtsar', 'python', 'utils'))
import gmtsar_lib


# ---------- check_file_report ----------

class TestCheckFileReport:

    def test_existing_file(self, tmp_path):
        f = tmp_path / "exists.txt"
        f.write_text("hello")
        assert gmtsar_lib.check_file_report(str(f)) is True

    def test_missing_file(self, tmp_path):
        assert gmtsar_lib.check_file_report(str(tmp_path / "nope.txt")) is False


# ---------- intFloatOrString ----------

class TestIntFloatOrString:

    def test_integer(self):
        assert gmtsar_lib.intFloatOrString("42") == 42

    def test_negative_integer_string(self):
        # "-5" is not purely digits, so isdigit() returns False;
        # float conversion succeeds.
        assert gmtsar_lib.intFloatOrString("-5") == -5.0

    def test_float(self):
        assert gmtsar_lib.intFloatOrString("3.14") == pytest.approx(3.14)

    def test_non_numeric(self):
        assert gmtsar_lib.intFloatOrString("abc") == ""

    def test_zero(self):
        assert gmtsar_lib.intFloatOrString("0") == 0


# ---------- grep_value ----------

class TestGrepValue:

    def test_grep_integer(self, tmp_path):
        f = tmp_path / "data.txt"
        f.write_text("key1 100\nkey2 200\n")
        val = gmtsar_lib.grep_value(str(f), "key2", 2)
        assert val == 200

    def test_grep_float(self, tmp_path):
        f = tmp_path / "data.txt"
        f.write_text("wavelength 5.6\n")
        val = gmtsar_lib.grep_value(str(f), "wavelength", 2)
        assert val == pytest.approx(5.6)


# ---------- replace_strings ----------

class TestReplaceStrings:

    def test_basic_replacement(self, tmp_path):
        f = tmp_path / "cfg.txt"
        f.write_text("DEFOMAX_CYCLE 0\nOTHER 1\n")
        gmtsar_lib.replace_strings(str(f), "DEFOMAX_CYCLE", "DEFOMAX_CYCLE 40")
        lines = f.read_text().splitlines()
        assert lines[0] == "DEFOMAX_CYCLE 40"
        assert lines[1] == "OTHER 1"

    def test_no_match_leaves_file_unchanged(self, tmp_path):
        f = tmp_path / "cfg.txt"
        original = "alpha 1\nbeta 2\n"
        f.write_text(original)
        gmtsar_lib.replace_strings(str(f), "gamma", "gamma 3")
        assert f.read_text() == original


# ---------- append_new_line ----------

class TestAppendNewLine:

    def test_append_to_empty_file(self, tmp_path):
        f = tmp_path / "out.txt"
        f.write_text("")
        gmtsar_lib.append_new_line(str(f), "first line")
        assert f.read_text() == "first line"

    def test_append_to_nonempty_file(self, tmp_path):
        f = tmp_path / "out.txt"
        f.write_text("line1")
        gmtsar_lib.append_new_line(str(f), "line2")
        content = f.read_text()
        assert "line1" in content
        assert "line2" in content
        assert "\n" in content


# ---------- file_shuttle ----------

class TestFileShuttle:

    def test_copy(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "dst.txt"
        src.write_text("data")
        gmtsar_lib.file_shuttle(str(src), str(dst), "cp")
        assert dst.read_text() == "data"
        assert src.exists()

    def test_move(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "dst.txt"
        src.write_text("data")
        gmtsar_lib.file_shuttle(str(src), str(dst), "mv")
        assert dst.read_text() == "data"
        assert not src.exists()

    def test_link(self, tmp_path):
        src = tmp_path / "src.txt"
        dst = tmp_path / "dst.txt"
        src.write_text("data")
        gmtsar_lib.file_shuttle(str(src), str(dst), "link")
        assert dst.read_text() == "data"


# ---------- delete ----------

class TestDelete:

    def test_delete_file(self, tmp_path):
        f = tmp_path / "todelete.txt"
        f.write_text("bye")
        gmtsar_lib.delete(str(f))
        assert not f.exists()

    def test_delete_nonexistent_no_error(self, tmp_path):
        # Should not raise even if file does not exist
        gmtsar_lib.delete(str(tmp_path / "ghost.txt"))


# ---------- assign_arg ----------

class TestAssignArg:

    def test_found_int(self):
        args = ["script.py", "--threshold", "42", "--output", "out.grd"]
        assert gmtsar_lib.assign_arg(args, "--threshold") == 42

    def test_found_float(self):
        args = ["script.py", "--corr", "0.12"]
        assert gmtsar_lib.assign_arg(args, "--corr") == pytest.approx(0.12)

    def test_not_found(self):
        args = ["script.py", "--other", "val"]
        assert gmtsar_lib.assign_arg(args, "--missing") == 0


# ---------- renameMasterAlignedForS1tops ----------

class TestRenameMasterAlignedForS1tops:

    def test_basic_rename(self):
        # Typical S1_TOPS filename format:
        # F1_xxx_xxx_xxx_20200101_123456_xxx...
        master0 = "F1_xxx_xxx_xxx_20200101T120000_123456_rest"
        aligned0 = "F2_xxx_xxx_xxx_20210601T060000_654321_rest"
        master, aligned = gmtsar_lib.renameMasterAlignedForS1tops(master0, aligned0)
        assert master.startswith("S1_")
        assert aligned.startswith("S1_")


# ---------- catch_output_cmd ----------

class TestCatchOutputCmd:

    def test_simple_command(self):
        result = gmtsar_lib.catch_output_cmd(["echo", "hello world"])
        assert result == "hello world"

    def test_split(self):
        result = gmtsar_lib.catch_output_cmd(["echo", "a b c"], choose_split=True)
        assert result == ["a", "b", "c"]

    def test_split_with_id(self):
        # split_id=2 means index 1 (split_id-1)
        result = gmtsar_lib.catch_output_cmd(["echo", "a b c"], choose_split=True, split_id=2)
        assert result == "b"

    def test_split_with_digit_id(self):
        # digit_id selects a character from the chosen split element
        # split_id=1 -> "hello", digit_id=1 -> "hello"[0] = "h"
        result = gmtsar_lib.catch_output_cmd(
            ["echo", "hello world"], choose_split=True, split_id=1, digit_id=1
        )
        assert result == "h"


# ---------- run ----------

class TestRun:

    def test_run_echoes_and_executes(self, tmp_path):
        f = tmp_path / "marker.txt"
        gmtsar_lib.run(f"touch {f}")
        assert f.exists()
