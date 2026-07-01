"""Crown sunlit-fraction (fixed from the former copy-paste stub)."""

import numpy as np
import pytest

from conftest import SCALE


def _fsun(pld, G, a, nrays=3000):
    return pld.crown_sunlit_fraction(G, a, 'ellipsoid', SCALE, SCALE, SCALE, 0.0,
                                     np.radians(90.0), nrays)


def test_thin_limit(pld_module):
    assert _fsun(pld_module, 0.5, 1e-4) > 0.99


def test_thick_limit(pld_module):
    assert _fsun(pld_module, 0.5, 50.0) < 0.01


def test_decreasing_in_lad(pld_module):
    vals = [_fsun(pld_module, 0.5, a) for a in [0.05, 0.2, 0.5, 1.0, 3.0]]
    assert all(b < a for a, b in zip(vals, vals[1:]))


def test_bounded(pld_module):
    for a in [1e-3, 0.1, 1.0, 10.0]:
        assert 0.0 <= _fsun(pld_module, 0.5, a) <= 1.0


def test_closed_form_two_ways(pld_module):
    # f_sun = sum_i (1-exp(-Ga r))/(Ga) / sum_i r, recomputed from raw paths.
    G, a = 0.5, 0.4
    pl, _ = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, 0.0,
                                   np.radians(90.0), 3000)
    Ga = G * a
    expected = np.sum((1 - np.exp(-Ga * pl)) / Ga) / np.sum(pl)
    got = _fsun(pld_module, G, a)
    assert got == pytest.approx(expected, rel=1e-9)


def test_not_equal_to_interception(pld_module):
    # Regression guard: catch re-introduction of the copy-paste bug where
    # crownsunlitfraction duplicated the absorption function.
    G, a = 0.5, 0.5
    fsun = pld_module.crownsunlitfraction(G, a, 'ellipsoid', SCALE, SCALE, SCALE,
                                          0.0, np.radians(90.0), 3000)
    Pabs = pld_module.crownabsorptionfraction(G, a, 'ellipsoid', SCALE, SCALE,
                                              SCALE, 0.0, np.radians(90.0), 3000)
    assert abs(fsun - Pabs) > 0.05
