"""Three-mode scattering absorbed-fraction model (Ponce de Leon et al. 2025)."""

import numpy as np
import pytest

from conftest import SCALE


def _abs(pld, LAD=1.0, G=0.5, rho_l=0.09, tau_l=0.04, rho_s=0.18, sr=6.0, sp=4.0,
         theta_deg=30.0, nrays=2000):
    return pld.absorbed_fraction(LAD, G, 'ellipsoid', SCALE, SCALE, SCALE,
                                 np.radians(theta_deg), 0.0, nrays, sr=sr, sp=sp,
                                 rho_l=rho_l, tau_l=tau_l, rho_s=rho_s)


def test_absorbed_bounded(pld_module):
    Q = _abs(pld_module)
    assert 0.0 <= Q <= 1.0


def test_ground_reflection_adds_energy(pld_module):
    Q0 = _abs(pld_module, rho_s=0.0)
    Q1 = _abs(pld_module, rho_s=0.3)
    assert Q1 > Q0


def test_monotonic_in_lad(pld_module):
    vals = [_abs(pld_module, LAD=a) for a in [0.2, 0.5, 1.0, 2.0]]
    assert all(b >= a for a, b in zip(vals, vals[1:]))


def test_higher_absorptivity_absorbs_more(pld_module):
    # Lower leaf reflectance/transmittance (higher zeta) -> more absorbed.
    # Use a thin canopy (LAD=0.1) so interception is not saturated at the
    # geometric ceiling s^2/(sr sp); otherwise zeta cannot matter.
    Q_low = _abs(pld_module, LAD=0.1, rho_l=0.4, tau_l=0.4, rho_s=0.0)   # zeta=0.2
    Q_high = _abs(pld_module, LAD=0.1, rho_l=0.05, tau_l=0.05, rho_s=0.0)  # zeta=0.9
    assert Q_high > Q_low


def test_no_scatter_matches_mode2_canopy(pld_module):
    # With zeta=1 (rho_l=tau_l=0) and rho_s=0, absorbed = P_c2 = canopy
    # interception with the doubled path (mode 2), by the paper's formulation.
    theta = np.radians(30.0)
    Q = pld_module.absorbed_fraction(1.0, 0.5, 'ellipsoid', SCALE, SCALE, SCALE,
                                     theta, 0.0, 2000, sr=6.0, sp=4.0,
                                     rho_l=0.0, tau_l=0.0, rho_s=0.0)
    Pc2 = pld_module.canopy_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                         SCALE, theta, 0.0, 2000, sr=6.0, sp=4.0,
                                         path_multiplier=2.0, absorptivity=1.0)
    assert Q == pytest.approx(Pc2, rel=1e-9)


def test_Q0_scales_linearly(pld_module):
    Qf = _abs(pld_module)
    Q = pld_module.absorbed_fraction(1.0, 0.5, 'ellipsoid', SCALE, SCALE, SCALE,
                                     np.radians(30.0), 0.0, 2000, sr=6.0, sp=4.0,
                                     rho_l=0.09, tau_l=0.04, rho_s=0.18, Q0=500.0)
    assert Q == pytest.approx(500.0 * Qf, rel=1e-9)
