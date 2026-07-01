"""Per-crown interception probability P_leaf (Bailey et al. 2020)."""

import numpy as np
import pytest

from conftest import R, SCALE


def _P(pld, G, a, mult=1.0, absorp=1.0, nrays=4000):
    return pld.crown_interception(G, a, 'ellipsoid', SCALE, SCALE, SCALE, 0.0,
                                  np.radians(90.0), nrays, path_multiplier=mult,
                                  absorptivity=absorp)


def test_beer_thin_limit(pld_module):
    # Thin canopy: P_leaf -> G*a*<r> = G*a*(4R/3).
    G, a = 0.5, 1e-3
    P = _P(pld_module, G, a)
    assert P == pytest.approx(G * a * 4.0 * R / 3.0, rel=3e-2)


def test_bounded_probability(pld_module):
    for a in [1e-3, 0.1, 1.0, 10.0]:
        P = _P(pld_module, 0.5, a)
        assert 0.0 <= P <= 1.0


def test_monotonic_in_lad(pld_module):
    vals = [_P(pld_module, 0.5, a) for a in [0.01, 0.1, 0.5, 1.0, 3.0]]
    assert all(b > a for a, b in zip(vals, vals[1:]))


def test_thick_limit_approaches_one(pld_module):
    assert _P(pld_module, 0.5, 100.0) > 0.999


def test_multiplier_two_ge_one(pld_module):
    P1 = _P(pld_module, 0.5, 0.2, mult=1.0)
    P2 = _P(pld_module, 0.5, 0.2, mult=2.0)
    assert P2 >= P1


def test_absorptivity_reduces_interception(pld_module):
    P_full = _P(pld_module, 0.5, 0.3, absorp=1.0)
    P_scat = _P(pld_module, 0.5, 0.3, absorp=0.87)  # zeta = 1 - 0.09 - 0.04
    assert P_scat < P_full


def test_crownabsorptionfraction_is_probability(pld_module):
    # The historically buggy /Ap division is fixed: result is a probability,
    # and equals crown_interception with unit multiplier/absorptivity.
    fa = pld_module.crownabsorptionfraction(0.5, 0.1, 'ellipsoid', SCALE, SCALE,
                                            SCALE, 0.0, np.radians(90.0), 4000)
    P = _P(pld_module, 0.5, 0.1)
    assert 0.0 <= fa <= 1.0
    assert fa == pytest.approx(P, rel=1e-9)
