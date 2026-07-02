"""Crown silhouette / projected areas.

pathlengths() returns the beam-normal projected area (pi R^2 at theta=0).  The
canopy model instead needs the isolated horizontal ground shadow S(theta) =
S_beam_normal / cos(theta) (Ponce de Leon et al. 2026, Table 2); that comes from
silhouette_area(), validated here against the analytic formulas (ellipsoid,
cylinder) and by a non-periodic mesh trace.
"""

import numpy as np
import pytest

from conftest import R, SCALE, NRAYS_FAST, PLY_SPHERE


def test_pathlengths_area_at_zenith0_is_pi_r2(pld_module):
    _, Ap = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, 0.0,
                                   np.radians(45.0), NRAYS_FAST)
    assert Ap == pytest.approx(np.pi * R ** 2, rel=1.5e-2)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0, 60.0])
def test_sphere_silhouette_is_horizontal_shadow(pld_module, deg):
    # silhouette_area returns the horizontal ground shadow (paper Table 2):
    # for a sphere that is pi R^2 / cos(theta), NOT the zenith-independent
    # beam-normal pi R^2.
    theta = np.radians(deg)
    S = pld_module.silhouette_area('ellipsoid', SCALE, SCALE, SCALE, theta)
    assert S == pytest.approx(np.pi * R ** 2 / np.cos(theta), rel=1e-9)


@pytest.mark.parametrize('deg', [0.0, 30.0, 45.0])
def test_ellipsoid_silhouette_formula(pld_module, deg):
    # Horizontal ground shadow = beam-normal silhouette / cos(theta):
    # S = pi R sqrt(R^2 cos^2 + b^2 sin^2) / cos, b = vertical semi-axis.
    theta = np.radians(deg)
    Rh, b = 5.0, 8.0
    S = pld_module.silhouette_area('ellipsoid', 2 * Rh, 2 * Rh, 2 * b, theta)
    expected = (np.pi * Rh * np.sqrt(Rh ** 2 * np.cos(theta) ** 2 +
                b ** 2 * np.sin(theta) ** 2) / np.cos(theta))
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
    # Non-periodic mesh horizontal shadow of a sphere ~ pi R^2 / cos(theta).
    theta = np.radians(deg)
    S = pld_module.silhouette_area('polymesh', SCALE, SCALE, SCALE, theta,
                                   nrays=NRAYS_FAST, plyfile=PLY_SPHERE)
    assert S == pytest.approx(np.pi * R ** 2 / np.cos(theta), rel=4e-2)


# --- azimuth dependence (triaxial ellipsoid, elliptical cylinder, prism) -----

def test_sphere_silhouette_azimuth_independent(pld_module):
    # A horizontally circular crown casts an azimuth-independent shadow.
    theta = np.radians(40.0)
    vals = [pld_module.silhouette_area('ellipsoid', SCALE, SCALE, SCALE, theta,
                                       np.radians(a))
            for a in [0.0, 30.0, 60.0, 90.0]]
    assert max(vals) - min(vals) < 1e-9


@pytest.mark.parametrize('azi', [0.0, 30.0, 45.0, 90.0])
@pytest.mark.parametrize('deg', [20.0, 45.0])
def test_triaxial_ellipsoid_silhouette_azimuth(pld_module, deg, azi):
    # Ground shadow of a triaxial ellipsoid = beam-normal projection / cos:
    #   pi sqrt(Rx^2 Ry^2 dz^2 + Ry^2 Rz^2 dx^2 + Rx^2 Rz^2 dy^2) / dz
    theta, phi = np.radians(deg), np.radians(azi)
    Rx, Ry, Rz = 4.0, 6.0, 8.0
    S = pld_module.silhouette_area('ellipsoid', 2 * Rx, 2 * Ry, 2 * Rz, theta, phi)
    dx, dy, dz = (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi),
                  np.cos(theta))
    perp = np.pi * np.sqrt(Rx ** 2 * Ry ** 2 * dz ** 2 +
                           Ry ** 2 * Rz ** 2 * dx ** 2 +
                           Rx ** 2 * Rz ** 2 * dy ** 2)
    assert S == pytest.approx(perp / dz, rel=1e-9)


def test_triaxial_ellipsoid_azimuth_actually_varies(pld_module):
    # Regression guard: for Rx != Ry the shadow must change with azimuth
    # (would be constant if silhouette_area ignored the beam azimuth).
    theta = np.radians(45.0)
    S0 = pld_module.silhouette_area('ellipsoid', 8.0, 12.0, 16.0, theta, 0.0)
    S90 = pld_module.silhouette_area('ellipsoid', 8.0, 12.0, 16.0, theta,
                                     np.radians(90.0))
    assert abs(S0 - S90) > 0.05 * S0


@pytest.mark.parametrize('azi', [0.0, 45.0, 90.0])
def test_triaxial_ellipsoid_silhouette_vs_mesh_azimuth(pld_module, azi):
    # Independent ground truth: non-periodic trace of a scaled sphere mesh at
    # the same beam azimuth.
    theta, phi = np.radians(45.0), np.radians(azi)
    Rx, Ry, Rz = 5.0, 8.0, 6.0
    S_mesh = pld_module.silhouette_area('polymesh', 2 * Rx, 2 * Ry, 2 * Rz,
                                        theta, phi, nrays=NRAYS_FAST,
                                        plyfile=PLY_SPHERE)
    S_analytic = pld_module.silhouette_area('ellipsoid', 2 * Rx, 2 * Ry, 2 * Rz,
                                            theta, phi)
    assert S_mesh == pytest.approx(S_analytic, rel=4e-2)


@pytest.mark.parametrize('azi', [0.0, 30.0, 45.0, 90.0])
def test_elliptical_cylinder_silhouette_azimuth(pld_module, azi):
    # pi Rx Ry + w(phi) H tan(theta), w(phi) = 2 sqrt(Rx^2 sin^2 + Ry^2 cos^2).
    theta, phi = np.radians(35.0), np.radians(azi)
    Rx, Ry, H = 4.0, 7.0, 9.0
    S = pld_module.silhouette_area('cylinder', 2 * Rx, 2 * Ry, H, theta, phi)
    width = 2.0 * np.sqrt(Rx ** 2 * np.sin(phi) ** 2 + Ry ** 2 * np.cos(phi) ** 2)
    expected = np.pi * Rx * Ry + width * H * np.tan(theta)
    assert S == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize('azi', [0.0, 30.0, 90.0])
def test_prism_silhouette_azimuth_matches_analytic(pld_module, azi):
    # Box ground shadow = Minkowski sum of the Lx x Ly footprint with the beam
    # shift (Sx, Sy) = H tan(theta) (cos phi, sin phi):
    #   Lx Ly + Ly |Sx| + Lx |Sy|
    theta, phi = np.radians(40.0), np.radians(azi)
    Lx, Ly, H = 6.0, 10.0, 8.0
    S = pld_module.silhouette_area('prism', Lx, Ly, H, theta, phi, nrays=90000)
    t = np.tan(theta)
    expected = (Lx * Ly + Ly * H * t * abs(np.cos(phi)) +
                Lx * H * t * abs(np.sin(phi)))
    assert S == pytest.approx(expected, rel=3e-2)
