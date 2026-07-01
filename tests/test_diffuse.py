"""Hemispherical (diffuse) interception integration (Bailey et al. 2020)."""

import numpy as np
import pytest

from conftest import SCALE


def test_constant_P_identity(pld_module, monkeypatch):
    # 2 * int_0^{pi/2} P0 cos sin dtheta = P0.  Patch canopy_interception -> 0.4.
    monkeypatch.setattr(pld_module, 'canopy_interception',
                        lambda *a, **k: 0.4)
    Pd = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                         SCALE, 100, sr=6.0, sp=4.0, n_zenith=16)
    assert Pd == pytest.approx(0.4, rel=5e-3)


def test_constant_P_identity_2d(pld_module, monkeypatch):
    monkeypatch.setattr(pld_module, 'canopy_interception',
                        lambda *a, **k: 0.4)
    Pd = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                         SCALE, 100, sr=6.0, sp=4.0, n_zenith=12,
                                         n_azimuth=8)
    assert Pd == pytest.approx(0.4, rel=5e-3)


def test_bracketed_by_directional_values(pld_module):
    # P_diff must lie between the min and max directional interception.
    Pd = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                         SCALE, 2000, sr=6.0, sp=4.0, n_zenith=8)
    dir_vals = [
        pld_module.canopy_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE, SCALE,
                                       th, 0.0, 2000, sr=6.0, sp=4.0)
        for th in np.radians([5.0, 30.0, 60.0, 85.0])
    ]
    assert min(dir_vals) - 1e-6 <= Pd <= max(dir_vals) + 1e-6


def test_quadrature_convergence(pld_module):
    kw = dict(sr=6.0, sp=4.0)
    coarse = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                             SCALE, 2000, n_zenith=6, **kw)
    fine = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                           SCALE, 2000, n_zenith=24, **kw)
    assert coarse == pytest.approx(fine, rel=1e-2)


def test_crown_level_diffuse_bounded(pld_module):
    Pd = pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                         SCALE, 2000, level='crown', n_zenith=8)
    assert 0.0 <= Pd <= 1.0


def test_canopy_level_requires_spacing(pld_module):
    with pytest.raises(ValueError):
        pld_module.diffuse_interception(0.5, 1.0, 'ellipsoid', SCALE, SCALE,
                                        SCALE, 100, level='canopy')
