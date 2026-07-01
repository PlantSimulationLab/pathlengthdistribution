"""Crown silhouette / projected areas.

pathlengths() returns the beam-normal silhouette S(theta).  At theta=0 this is
the exact silhouette; at oblique angles the periodic boundaries inflate it (a ray
can cross the tiled crown multiple times), so the *isolated* S(theta) needed by
the canopy model comes from silhouette_area(), which is validated here against
the analytic formulas (ellipsoid, cylinder) and by a non-periodic mesh trace.
"""

import numpy as np
import pytest

from conftest import R, SCALE, NRAYS_FAST, PLY_SPHERE


def test_pathlengths_area_at_zenith0_is_pi_r2(pld_module):
    _, Ap = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, 0.0,
                                   np.radians(45.0), NRAYS_FAST)
    assert Ap == pytest.approx(np.pi * R ** 2, rel=1.5e-2)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0, 60.0])
def test_sphere_silhouette_zenith_independent(pld_module, deg):
    # Analytic silhouette_area for a sphere is pi R^2 at all zenith angles.
    S = pld_module.silhouette_area('ellipsoid', SCALE, SCALE, SCALE,
                                   np.radians(deg))
    assert S == pytest.approx(np.pi * R ** 2, rel=1e-9)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0])
def test_ellipsoid_silhouette_formula(pld_module, deg):
    # Norman & Welles (1983): S = pi R sqrt(R^2 cos^2 + b^2 sin^2), b = vertical semi-axis.
    theta = np.radians(deg)
    Rh, b = 5.0, 8.0
    S = pld_module.silhouette_area('ellipsoid', 2 * Rh, 2 * Rh, 2 * b, theta)
    expected = np.pi * Rh * np.sqrt(Rh ** 2 * np.cos(theta) ** 2 + b ** 2 * np.sin(theta) ** 2)
    assert S == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0, 60.0])
def test_ellipsoid_silhouette_vs_mesh(pld_module, deg):
    # Independent ground truth: non-periodic trace of a scaled sphere mesh.
    theta = np.radians(deg)
    Rh, b = 5.0, 8.0
    S_mesh = pld_module.silhouette_area('polymesh', 2 * Rh, 2 * Rh, 2 * b, theta,
                                        nrays=NRAYS_FAST, plyfile=PLY_SPHERE)
    S_analytic = pld_module.silhouette_area('ellipsoid', 2 * Rh, 2 * Rh, 2 * b, theta)
    assert S_mesh == pytest.approx(S_analytic, rel=3e-2)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0])
def test_cylinder_silhouette_formula(pld_module, deg):
    theta = np.radians(deg)
    Rc, H = 5.0, 10.0
    S = pld_module.silhouette_area('cylinder', 2 * Rc, 2 * Rc, H, theta)
    expected = np.pi * Rc ** 2 + 2.0 * Rc * H * np.tan(theta)
    assert S == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0, 60.0])
def test_mesh_sphere_silhouette_matches_analytic(pld_module, deg):
    # Non-periodic mesh silhouette of a sphere should be ~pi R^2 at all angles.
    theta = np.radians(deg)
    S = pld_module.silhouette_area('polymesh', SCALE, SCALE, SCALE, theta,
                                   nrays=NRAYS_FAST, plyfile=PLY_SPHERE)
    assert S == pytest.approx(np.pi * R ** 2, rel=4e-2)
