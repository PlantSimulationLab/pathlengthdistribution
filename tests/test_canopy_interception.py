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


# --- closure-corrected ceiling (replaces the Eq. 13 s^2/(sr sp) prefactor) ---

def test_nadir_is_exact_cover_fraction(pld_module):
    # At theta=0 (Nc=1) Pc = S(0) P_leaf/(sr sp) exactly, independent of ceiling.
    sr, sp = 6.0, 4.0
    S0 = pld_module.silhouette_area('ellipsoid', 4.0, 6.0, 6.0, 0.0)
    Pl = pld_module.crown_interception(0.5, 1.0, 'ellipsoid', 4.0, 6.0, 6.0,
                                       0.0, 0.0, 6400)
    Pc = pld_module.canopy_interception(0.5, 1.0, 'ellipsoid', 4.0, 6.0, 6.0,
                                        0.0, 0.0, 6400, sr=sr, sp=sp)
    assert Pc == pytest.approx(S0 / (sr * sp) * Pl, rel=1e-9)


def test_closed_canopy_rises_to_full_at_grazing(pld_module):
    # Touching crowns filling the cell must climb with zenith and approach ~1
    # at grazing for BOTH azimuths.  The published s^2/(sr sp) prefactor pinned
    # azimuth 0 at sp/sr = 2/3.
    for azi in [0.0, 90.0]:
        vals = [pld_module.canopy_interception(
            0.5, 1.0, 'ellipsoid', 4.0, 6.0, 6.0, np.radians(z),
            np.radians(azi), 6400, sr=6.0, sp=4.0)
            for z in [10.0, 40.0, 70.0, 85.0]]
        assert all(b > a for a, b in zip(vals, vals[1:]))
        assert vals[-1] > 0.97


def test_grazing_ceiling_is_geometric(pld_module):
    # As theta->90 a dense canopy saturates at C = min(1, w_perp/s_perp).
    # crown 2x2x8, sr=8, sp=2, azimuth 0 -> w_perp=Dy=2, s_perp=sr=8 -> C=0.25.
    Pc = pld_module.canopy_interception(0.5, 5.0, 'ellipsoid', 2.0, 2.0, 8.0,
                                        np.radians(88.0), 0.0, 6400,
                                        sr=8.0, sp=2.0)
    assert Pc == pytest.approx(0.25, abs=0.02)


def test_row_anisotropy_preserved(pld_module):
    # Row-induced anisotropy must survive: a tall crown with very different
    # sr, sp intercepts very differently along vs across rows.
    th = np.radians(60.0)
    a0 = pld_module.canopy_interception(0.5, 1.0, 'ellipsoid', 2.0, 2.0, 8.0,
                                        th, 0.0, 6400, sr=8.0, sp=2.0)
    a90 = pld_module.canopy_interception(0.5, 1.0, 'ellipsoid', 2.0, 2.0, 8.0,
                                         th, np.radians(90.0), 6400,
                                         sr=8.0, sp=2.0)
    assert abs(a0 - a90) > 0.2
