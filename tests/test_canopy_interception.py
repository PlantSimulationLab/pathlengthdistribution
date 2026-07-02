"""Canopy-level binomial interception P_c (Bailey et al. 2020)."""

import numpy as np
import pytest

from conftest import SCALE


def _Pc(pld, G, a, sr, sp, theta_deg=20.0, azi_deg=0.0, phi=None, nrays=3000):
    return pld.canopy_interception(G, a, 'ellipsoid', SCALE, SCALE, SCALE,
                                   np.radians(theta_deg), np.radians(azi_deg),
                                   nrays, sr=sr, sp=sp, phi=phi)


def test_bounded(pld_module):
    for sr in [1.0, 4.0, 12.0]:
        Pc = _Pc(pld_module, 0.5, 1.0, sr, 2.0)
        assert 0.0 <= Pc <= 1.0


def test_dense_canopy_approaches_one(pld_module):
    Pc = _Pc(pld_module, 0.5, 5.0, sr=1.2, sp=1.2, theta_deg=5.0)
    assert Pc > 0.9


def test_monotonic_decreasing_in_row_spacing(pld_module):
    # Spacings must exceed the crown (SCALE=10) in both directions, else the
    # canopy is fully closed (Pc saturates at 1.0) and monotonicity is masked.
    vals = [_Pc(pld_module, 0.5, 1.0, sr, 14.0) for sr in [12.0, 16.0, 24.0, 32.0]]
    assert all(b < a for a, b in zip(vals, vals[1:]))


def test_azimuth_symmetry_when_sr_equals_sp(pld_module):
    # With sr == sp the effective spacing s is azimuth-independent.
    a = _Pc(pld_module, 0.5, 1.0, 3.0, 3.0, azi_deg=0.0, phi=0.0)
    b = _Pc(pld_module, 0.5, 1.0, 3.0, 3.0, azi_deg=0.0, phi=np.radians(45.0))
    assert a == pytest.approx(b, rel=1e-9)


def test_sparse_limit_first_order(pld_module):
    # Very sparse: P_c ~ (S0/(sr sp)) * P_leaf to first order (Nc*(S0/s^2)*P_leaf,
    # with s^2/(sr sp) factor cancelling to S0*Nc/(sr sp) = S(theta)/(sr sp)*P_leaf).
    G, a, sr, sp = 0.5, 0.05, 30.0, 30.0
    theta = np.radians(20.0)
    Pc = pld_module.canopy_interception(G, a, 'ellipsoid', SCALE, SCALE, SCALE,
                                        theta, 0.0, 4000, sr=sr, sp=sp)
    Pleaf = pld_module.crown_interception(G, a, 'ellipsoid', SCALE, SCALE, SCALE,
                                          theta, 0.0, 4000)
    Stheta = pld_module.silhouette_area('ellipsoid', SCALE, SCALE, SCALE, theta)
    approx = Stheta / (sr * sp) * Pleaf
    assert Pc == pytest.approx(approx, rel=5e-2)
