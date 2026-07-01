"""Regression: numba-accelerated kernels reproduce the pure-Python reference.

The reference implementation is a frozen copy of the original module
(tests/_reference_impl.py).  Because the ray grid is a deterministic lattice,
both implementations trace the same rays in the same order.  Analytic shapes
should agree to floating-point noise; the polymesh path differs only by tiny
reassociation error, so a modest tolerance is used.
"""

import numpy as np
import pytest

import _reference_impl as ref
from conftest import SCALE, PLY_CONE


ANALYTIC = ['ellipsoid', 'cylinder', 'prism']


@pytest.mark.parametrize('shape', ANALYTIC)
def test_analytic_paths_match_reference(pld_module, shape):
    theta, azi = np.radians(35.0), np.radians(70.0)
    plr, Apr = ref.pathlengths(shape, 10.0, 7.0, 12.0, theta, azi, 2000)
    pln, Apn = pld_module.pathlengths(shape, 10.0, 7.0, 12.0, theta, azi, 2000)
    assert pln.size == plr.size
    assert Apn == pytest.approx(Apr, rel=1e-9)
    np.testing.assert_allclose(np.sort(pln), np.sort(plr), atol=1e-6)


def test_polymesh_paths_match_reference(pld_module):
    theta, azi = np.radians(30.0), np.radians(90.0)
    plr, Apr = ref.pathlengths('polymesh', SCALE, SCALE, SCALE, theta, azi, 500,
                               plyfile=PLY_CONE)
    pln, Apn = pld_module.pathlengths('polymesh', SCALE, SCALE, SCALE, theta, azi,
                                      500, plyfile=PLY_CONE)
    assert pln.size == plr.size
    assert Apn == pytest.approx(Apr, rel=1e-6)
    # float reassociation between pure-Python and numba: ~1e-5 worst case.
    np.testing.assert_allclose(np.sort(pln), np.sort(plr), atol=1e-4)
