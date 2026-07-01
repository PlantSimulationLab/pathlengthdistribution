"""Shared fixtures and constants for the pathlengthdistribution test suite.

Determinism
-----------
``pathlengths`` launches a fixed ``N x N`` lattice of rays (``N = ceil(sqrt(nrays))``)
with no random sampling, so every result is exactly reproducible.  Tests pick
``nrays`` as perfect squares so ``N`` is exact, and several tests assert
``np.array_equal`` across repeated calls to lock in this invariant.
"""

import os
import sys

import numpy as np
import pytest

# Make the top-level module importable when pytest is run from anywhere.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import pathlengthdistribution as pld  # noqa: E402


# Canonical test crown: an ellipsoid with equal scales == sphere of radius R.
R = 5.0
SCALE = 2.0 * R          # 10.0  (full axis length; semi-axis = SCALE/2 = R)
NRAYS = 40000            # perfect square -> N = 200
NRAYS_FAST = 10000       # N = 100, for cheaper tests

PLY_SPHERE = os.path.join(REPO_ROOT, 'PLY', 'sphere.ply')
PLY_CONE = os.path.join(REPO_ROOT, 'PLY', 'cone.ply')


@pytest.fixture(scope='session')
def pld_module():
    return pld


@pytest.fixture(scope='session', autouse=True)
def _warmup():
    """Trigger numba JIT compilation once so per-test timing/behavior is clean."""
    for shape in ('ellipsoid', 'cylinder', 'prism'):
        pld.pathlengths(shape, SCALE, SCALE, SCALE, 0.1, 0.0, 16)
    pld.pathlengths('polymesh', SCALE, SCALE, SCALE, 0.1, 0.0, 16,
                    plyfile=PLY_CONE)
    yield
