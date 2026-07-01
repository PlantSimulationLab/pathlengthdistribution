"""Analytical checks for the path-length distribution of a sphere.

For a sphere of radius R the path-length PDF is triangular,
    p(r) = r / (2 R^2),   0 <= r <= 2R,
with CDF F(r) = (r/2R)^2, mean 4R/3, and E[r^2] = 2R^2 (Bailey et al. 2020).
"""

import numpy as np
import pytest

from conftest import R, SCALE, NRAYS


def _sphere_paths(pld, zenith=0.0, azimuth=np.pi / 2, nrays=NRAYS):
    pl, Ap = pld.pathlengths('ellipsoid', SCALE, SCALE, SCALE, zenith, azimuth,
                             nrays)
    return pl, Ap


def test_mean_path_length(pld_module):
    pl, _ = _sphere_paths(pld_module)
    assert pl.mean() == pytest.approx(4.0 * R / 3.0, rel=1e-2)


def test_max_path_length(pld_module):
    pl, _ = _sphere_paths(pld_module)
    assert pl.max() == pytest.approx(2.0 * R, abs=0.05)


def test_second_moment(pld_module):
    pl, _ = _sphere_paths(pld_module)
    # E[r^2] = int_0^{2R} r^2 * r/(2R^2) dr = 2 R^2
    assert np.mean(pl ** 2) == pytest.approx(2.0 * R ** 2, rel=2e-2)


def test_triangular_cdf_ks(pld_module):
    pl, _ = _sphere_paths(pld_module)

    def cdf(r):
        return np.clip((r / (2.0 * R)) ** 2, 0.0, 1.0)

    try:
        from scipy import stats
        D = stats.kstest(pl, cdf).statistic
    except ImportError:
        # Manual KS statistic if scipy is unavailable.
        s = np.sort(pl)
        n = s.size
        emp = np.arange(1, n + 1) / n
        D = np.max(np.abs(emp - cdf(s)))
    assert D < 0.03


def test_projected_area_is_pi_r_squared(pld_module):
    # At theta=0 the returned area is the exact silhouette pi R^2.
    _, Ap = _sphere_paths(pld_module, zenith=0.0)
    assert Ap == pytest.approx(np.pi * R ** 2, rel=1e-2)


def test_determinism(pld_module):
    a, _ = _sphere_paths(pld_module, nrays=10000)
    b, _ = _sphere_paths(pld_module, nrays=10000)
    assert np.array_equal(a, b)
