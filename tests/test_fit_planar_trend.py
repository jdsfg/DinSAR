"""Unit tests for gmtsar-master/gmtsar/csh/fit_planar_trend.py"""

import os
import sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'gmtsar-master', 'gmtsar', 'csh'))
from fit_planar_trend import get_trend


class TestGetTrend:

    def test_flat_plane(self):
        # z = 5 everywhere; x and y must be linearly independent
        rng = np.random.RandomState(0)
        x = rng.uniform(0, 10, 20)
        y = rng.uniform(0, 10, 20)
        z = np.full_like(x, 5.0)
        p = get_trend(x, y, z)
        assert p[0] == pytest.approx(5.0, abs=1e-8)
        assert p[1] == pytest.approx(0.0, abs=1e-8)
        assert p[2] == pytest.approx(0.0, abs=1e-8)

    def test_linear_x_slope(self):
        # z = 2 + 3*x; y independent of x
        rng = np.random.RandomState(1)
        x = rng.uniform(0, 10, 30)
        y = rng.uniform(0, 10, 30)
        z = 2.0 + 3.0 * x
        p = get_trend(x, y, z)
        assert p[0] == pytest.approx(2.0, abs=1e-6)
        assert p[1] == pytest.approx(3.0, abs=1e-6)
        assert p[2] == pytest.approx(0.0, abs=1e-6)

    def test_linear_y_slope(self):
        # z = 1 + 4*y
        rng = np.random.RandomState(2)
        x = rng.uniform(0, 10, 30)
        y = rng.uniform(0, 10, 30)
        z = 1.0 + 4.0 * y
        p = get_trend(x, y, z)
        assert p[0] == pytest.approx(1.0, abs=1e-6)
        assert p[1] == pytest.approx(0.0, abs=1e-6)
        assert p[2] == pytest.approx(4.0, abs=1e-6)

    def test_full_planar_trend(self):
        # z = 10 + 2*x - 3*y
        rng = np.random.RandomState(42)
        x = rng.uniform(0, 100, 50)
        y = rng.uniform(0, 100, 50)
        z = 10.0 + 2.0 * x - 3.0 * y
        p = get_trend(x, y, z)
        assert p[0] == pytest.approx(10.0, abs=1e-6)
        assert p[1] == pytest.approx(2.0, abs=1e-6)
        assert p[2] == pytest.approx(-3.0, abs=1e-6)

    def test_returns_three_coefficients(self):
        rng = np.random.RandomState(3)
        x = rng.uniform(0, 5, 10)
        y = rng.uniform(0, 5, 10)
        z = 7.0 + 1.5 * x - 2.0 * y
        p = get_trend(x, y, z)
        assert len(p) == 3
