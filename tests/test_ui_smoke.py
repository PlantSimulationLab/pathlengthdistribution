"""Headless smoke tests for the GUI's compute layer (`pathlength_ui.run_model`).

These do not open a window: importing ``pathlength_ui`` only wires up the
compute helpers (the Tk/matplotlib app is built lazily inside ``main``), so the
module imports cleanly with no display. The tests confirm ``run_model`` returns
the expected keys and that its numbers match the analytic sphere anchors used
elsewhere in the suite.
"""

import numpy as np
import pytest

from conftest import R, SCALE, NRAYS_FAST, PLY_SPHERE

import pathlength_ui as ui


def _base_params(**overrides):
    params = {
        'shape': 'sphere',
        'scale_x': SCALE, 'scale_y': SCALE, 'scale_z': SCALE,
        'zenith': 0.0, 'azimuth': 90.0,
        'nrays': NRAYS_FAST,
        'plyfile': '',
        'bins': 15, 'normalize': True,
        'Gtheta': 0.5, 'LAD': 1.0,
        # Canopy geometry is always supplied.
        'sr': 6.0, 'sp': 4.0, 'phi': None,
        'include_scatter': False,
        'include_diffuse': False,
    }
    params.update(overrides)
    return params


def test_import_is_side_effect_free():
    # Importing must not construct any Tk objects; the compute helpers exist.
    assert hasattr(ui, 'run_model')
    assert callable(ui.run_model)


def test_distribution_matches_sphere_analytics():
    res = ui.run_model(_base_params())
    assert res['n_intersecting'] > 0
    assert res['mean_path'] == pytest.approx(4.0 * R / 3.0, rel=2e-2)
    assert res['max_path'] == pytest.approx(2.0 * R, abs=0.1)
    # At theta=0 the projected area is the exact silhouette pi R^2.
    assert res['projected_area'] == pytest.approx(np.pi * R ** 2, rel=2e-2)


def test_interception_components_always_present_and_bounded():
    # Crown and canopy interception are the headline outputs and are computed
    # on every run, no toggle required.
    res = ui.run_model(_base_params(zenith=30.0))
    for key in ('S_theta', 'S_zero', 'P_leaf', 'P_canopy', 'sunlit_fraction'):
        assert key in res
    assert 0.0 <= res['P_leaf'] <= 1.0
    assert 0.0 <= res['P_canopy'] <= 1.0
    assert 0.0 <= res['sunlit_fraction'] <= 1.0
    assert res['S_theta'] > 0.0 and res['S_zero'] > 0.0


def test_scatter_absorbed_bounded():
    res = ui.run_model(_base_params(
        include_scatter=True, zenith=30.0, azimuth=0.0,
        rho_l=0.09, tau_l=0.04, rho_s=0.18, Q0=1.0))
    assert 'absorbed' in res
    assert 0.0 <= res['absorbed'] <= 1.0
    assert res['zeta'] == pytest.approx(1.0 - 0.09 - 0.04)


def test_diffuse_canopy_level():
    res = ui.run_model(_base_params(
        include_diffuse=True, diffuse_level='canopy', n_zenith=8))
    assert 0.0 <= res['P_diffuse'] <= 1.0
    assert res['diffuse_level'] == 'canopy'


def test_polymesh_requires_plyfile():
    with pytest.raises(Exception):
        ui.run_model(_base_params(shape='polymesh', plyfile=''))


def test_polymesh_with_ply_runs():
    res = ui.run_model(_base_params(shape='polymesh', plyfile=PLY_SPHERE,
                                    zenith=20.0, azimuth=0.0))
    assert res['n_intersecting'] > 0
    assert res['max_path'] == pytest.approx(2.0 * R, abs=0.5)
