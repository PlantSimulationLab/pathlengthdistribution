"""Triangular-mesh (polymesh) checks and error handling."""

import numpy as np
import pytest

from conftest import R, SCALE, NRAYS, PLY_SPHERE, PLY_CONE


def test_mesh_sphere_matches_analytic_sphere(pld_module):
    theta, azi = np.radians(30.0), np.radians(45.0)
    plm, Apm = pld_module.pathlengths('polymesh', SCALE, SCALE, SCALE, theta, azi,
                                      NRAYS, plyfile=PLY_SPHERE)
    pla, Apa = pld_module.pathlengths('ellipsoid', SCALE, SCALE, SCALE, theta, azi,
                                      NRAYS)
    # Coarse mesh (960 triangles) -> loose tolerance, but means & areas agree.
    assert plm.mean() == pytest.approx(pla.mean(), rel=4e-2)
    assert Apm == pytest.approx(Apa, rel=4e-2)


def test_mesh_sphere_mean_is_4R_over_3(pld_module):
    pl, _ = pld_module.pathlengths('polymesh', SCALE, SCALE, SCALE, 0.0,
                                   np.radians(90.0), NRAYS, plyfile=PLY_SPHERE)
    assert pl.mean() == pytest.approx(4.0 * R / 3.0, rel=3e-2)


def test_cone_alias_routes_to_polymesh(pld_module):
    # 'cone' is an alias for polymesh and needs a PLY file.
    pl, Ap = pld_module.pathlengths('cone', SCALE, SCALE, SCALE, np.radians(20.0),
                                    0.0, 4000, plyfile=PLY_CONE)
    assert pl.size > 0 and Ap > 0.0


def test_missing_plyfile_raises(pld_module):
    with pytest.raises(Exception):
        pld_module.pathlengths('polymesh', SCALE, SCALE, SCALE, 0.0, 0.0, 100)


def test_nonexistent_plyfile_raises(pld_module):
    with pytest.raises(Exception):
        pld_module.pathlengths('polymesh', SCALE, SCALE, SCALE, 0.0, 0.0, 100,
                               plyfile='does_not_exist.ply')
