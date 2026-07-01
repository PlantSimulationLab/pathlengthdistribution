"""Shape aliases, unit handling, and input validation."""

import numpy as np
import pytest

from conftest import SCALE, PLY_CONE


def test_sphere_alias_equals_ellipsoid(pld_module):
    a, Aa = pld_module.pathlengths('sphere', SCALE, SCALE, SCALE, np.radians(20.0),
                                   np.radians(30.0), 4000)
    b, Ab = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE,
                                   np.radians(20.0), np.radians(30.0), 4000)
    assert np.array_equal(a, b)
    assert Aa == Ab


def test_degrees_equals_radians(pld_module):
    a, _ = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, 45.0, 90.0,
                                  4000, degrees=True)
    b, _ = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE,
                                  np.radians(45.0), np.radians(90.0), 4000)
    assert np.array_equal(a, b)


def test_case_insensitive_shape(pld_module):
    a, _ = pld_module.pathlengths('ELLIPSOID', SCALE, SCALE, SCALE, 0.1, 0.0, 1000)
    b, _ = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, 0.1, 0.0, 1000)
    assert np.array_equal(a, b)


def test_invalid_shape_raises(pld_module):
    with pytest.raises(ValueError):
        pld_module.pathlengths('banana', SCALE, SCALE, SCALE, 0.1, 0.0, 100)


def test_cone_alias_needs_ply(pld_module):
    with pytest.raises(Exception):
        pld_module.pathlengths('cone', SCALE, SCALE, SCALE, 0.1, 0.0, 100)
    # ...and works with one.
    pl, _ = pld_module.pathlengths('cone', SCALE, SCALE, SCALE, 0.1, 0.0, 1000,
                                   plyfile=PLY_CONE)
    assert pl.size >= 0
