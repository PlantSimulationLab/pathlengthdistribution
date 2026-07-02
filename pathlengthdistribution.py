"""Radiative path-length distributions and geometric radiation interception.

This module computes the distribution of radiative path lengths through 3D plant
crown shapes (sphere/ellipsoid, cylinder, rectangular prism, or an arbitrary
triangular mesh) and, from those distributions, the geometric radiation
interception and absorption quantities described by:

  * Bailey, Ponce de Leon & Krayenhoff (2020), Geosci. Model Dev. 13, 4789-4808
    -- path-length distributions, the canopy-level binomial interception model,
    and hemispherical diffuse integration.
  * Ponce de Leon et al. (2025), Agric. For. Meteorol. 373, 110706
    -- three-mode scattering (leaf reflectance/transmittance + ground reflection).
  * Ponce de Leon et al. (2026), J. Geophys. Res. Biogeosciences
    -- ellipsoidal crowns with diffuse radiation and scattering.

Conventions
-----------
* ``ray_zenith`` and ``ray_azimuth`` are in **radians** by default.  Pass
  ``degrees=True`` to the public functions to supply them in degrees instead.
* Shapes: ``'ellipsoid'`` (equal scales give a sphere), ``'cylinder'``,
  ``'prism'`` (rectangular prism), and ``'polymesh'`` (triangular mesh from a
  PLY file).  Aliases ``'sphere'`` -> ``'ellipsoid'`` and ``'cone'`` ->
  ``'polymesh'`` are accepted.
* The ray grid is a deterministic ``N x N`` lattice (``N = ceil(sqrt(nrays))``);
  there is no random sampling, so results are fully reproducible.

Performance
-----------
The ray-marching and ray/triangle intersection kernels are JIT-compiled with
numba.  The pure-Python reference implementations are retained (suffixed
``_ref``) for regression testing and readability.
"""

import os
import numpy as np
from numpy import sqrt, sin, cos, exp, pi, ceil

try:
    from numba import njit
    _HAVE_NUMBA = True
except ImportError:  # pragma: no cover - numba is a listed dependency
    _HAVE_NUMBA = False

    def njit(*args, **kwargs):
        """Fallback no-op decorator if numba is unavailable."""
        def _decorate(func):
            return func
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return _decorate

from plyfile import PlyData, PlyElement


# ---------------------------------------------------------------------------
# Shape-name / unit helpers
# ---------------------------------------------------------------------------

_SHAPE_ALIASES = {
    'ellipsoid': 'ellipsoid',
    'sphere': 'ellipsoid',
    'cylinder': 'cylinder',
    'prism': 'prism',
    'box': 'prism',
    'rectangular_prism': 'prism',
    'polymesh': 'polymesh',
    'mesh': 'polymesh',
    'cone': 'polymesh',
}


def _normalize_shape(shape):
    """Map a user-facing shape name/alias to an internal shape key.

    Accepted (case-insensitive): 'sphere'/'ellipsoid', 'cylinder',
    'prism'/'box'/'rectangular_prism', 'polymesh'/'mesh'/'cone'.

    ``'sphere'`` is an ellipsoid with equal scales; ``'cone'`` is handled via a
    triangular mesh (there is no analytic cone primitive), so it maps to
    'polymesh' and requires a ``plyfile``.
    """
    key = str(shape).strip().lower()
    if key not in _SHAPE_ALIASES:
        raise ValueError(
            "Invalid shape '{}'. Options: prism, ellipsoid (alias sphere), "
            "cylinder, polymesh (alias cone).".format(shape))
    return _SHAPE_ALIASES[key]


def _to_radians(angle, degrees):
    """Return ``angle`` in radians, converting from degrees if requested."""
    return np.radians(angle) if degrees else float(angle)


# ===========================================================================
# Numba-accelerated numerical kernels
# ===========================================================================
#
# These are faithful ports of the pure-Python reference functions further down
# in this file.  They operate only on plain float64/int64 scalars and arrays
# (numba cannot consume the structured record arrays returned by plyfile), so
# mesh geometry is pre-extracted into a contiguous (Nfaces, 3, 3) float64 array
# once in ``pathlengths`` before these kernels are called.
#
# A BVH / spatial acceleration structure is intentionally not used: brute-force
# ray/triangle testing under numba is fast enough for meshes up to ~1e5 faces.
# For much larger meshes a BVH would be the natural next step.


@njit(cache=True, fastmath=False)
def _intersect_bbox(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):
    """Ray/axis-aligned-box slab intersection (Suffern 2007, Listing 19.1).

    Returns ``(dr, xe, ye, ze)`` where ``dr`` is the segment length inside the
    box and ``(xe, ye, ze)`` is the exit point.
    """
    x0 = -0.5 * sizex
    x1 = 0.5 * sizex
    y0 = -0.5 * sizey
    y1 = 0.5 * sizey
    z0 = -1e-6
    z1 = sizez

    if dx == 0.0:
        a = 1e6
    else:
        a = 1.0 / dx
    if a >= 0.0:
        tx_min = (x0 - ox) * a
        tx_max = (x1 - ox) * a
    else:
        tx_min = (x1 - ox) * a
        tx_max = (x0 - ox) * a

    if dy == 0.0:
        b = 1e6
    else:
        b = 1.0 / dy
    if b >= 0.0:
        ty_min = (y0 - oy) * b
        ty_max = (y1 - oy) * b
    else:
        ty_min = (y1 - oy) * b
        ty_max = (y0 - oy) * b

    if dz == 0.0:
        c = 1e6
    else:
        c = 1.0 / dz
    if c >= 0.0:
        tz_min = (z0 - oz) * c
        tz_max = (z1 - oz) * c
    else:
        tz_min = (z1 - oz) * c
        tz_max = (z0 - oz) * c

    if tx_min > ty_min:
        t0 = tx_min
    else:
        t0 = ty_min
    if tz_min > t0:
        t0 = tz_min

    if tx_max < ty_max:
        t1 = tx_max
    else:
        t1 = ty_max
    if tz_max < t1:
        t1 = tz_max

    if t0 < t1 and t1 > 1e-6:
        if t0 > 1e-6:
            dr = t1 - t0
        else:
            dr = t1
    else:
        dr = 0.0

    xe = ox + t1 * dx
    ye = oy + t1 * dy
    ze = oz + t1 * dz

    return dr, xe, ye, ze


@njit(cache=True, fastmath=False)
def _intersect_ellipsoid(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):
    """Path length through an ellipsoid centered at height ``sizez/2``."""
    tempx = ox / sizex
    tempy = oy / sizey
    tempz = (oz - 0.5) / sizez

    ddx = dx / sizex
    ddy = dy / sizey
    ddz = dz / sizez

    a = ddx * ddx + ddy * ddy + ddz * ddz
    b = 2.0 * (tempx * ddx + tempy * ddy + tempz * ddz)
    c = (tempx * tempx + tempy * tempy + tempz * tempz) - 0.25
    disc = b * b - 4.0 * a * c

    if disc < 0.0:
        return 0.0
    e = sqrt(disc)
    denom = 2.0 * a
    t_small = (-b - e) / denom
    t_big = (-b + e) / denom
    if t_small > 1e-6 or t_big > 1e-6:
        return abs(t_big - t_small)
    return 0.0


@njit(cache=True, fastmath=False)
def _intersect_cylinder(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):
    """Path length through an upright cylinder (Suffern 2007, 19.5.3)."""
    tempx = ox / sizex
    tempy = oy / sizey
    tempz = oz / sizez

    ddx = dx / sizex
    ddy = dy / sizey
    ddz = dz / sizez

    a = ddx * ddx + ddy * ddy
    b = 2.0 * (tempx * ddx + tempy * ddy)
    c = (tempx * tempx + tempy * tempy) - 0.25

    disc = b * b - 4.0 * a * c
    if disc >= 0.0 and a != 0.0:
        t0 = (-b + sqrt(disc)) / (2.0 * a)
        t1 = (-b - sqrt(disc)) / (2.0 * a)
    else:
        t0 = -1.0
        t1 = -1.0

    # Candidate intersection parameters (side x2, top, bottom); -1 = no hit.
    c0 = -1.0
    c1 = -1.0
    c2 = -1.0
    c3 = -1.0

    z0 = tempz + t0 * ddz
    if t0 > 1e-6 and 0.0 <= z0 <= 1.0:
        c0 = t0

    z1 = tempz + t1 * ddz
    if t1 > 1e-6 and 0.0 <= z1 <= 1.0:
        c1 = t1

    if ddz != 0.0:
        t_top = (1.0 - tempz) / ddz
        hx = tempx + t_top * ddx
        hy = tempy + t_top * ddy
        if t_top > 0.0 and hx * hx + hy * hy <= 0.25:
            c2 = t_top

        t_bot = (0.0 - tempz) / ddz
        hx = tempx + t_bot * ddx
        hy = tempy + t_bot * ddy
        if t_bot > 0.0 and hx * hx + hy * hy <= 0.25:
            c3 = t_bot

    # Smallest non-negative entry and overall largest, matching the reference
    # (Tin = min over non-negative candidates, Tout = max over all four).
    tin = 1e30
    have_in = False
    for cval in (c0, c1, c2, c3):
        if cval >= 0.0 and cval < tin:
            tin = cval
            have_in = True
    if not have_in:
        return 0.0

    tout = c0
    if c1 > tout:
        tout = c1
    if c2 > tout:
        tout = c2
    if c3 > tout:
        tout = c3

    if 0.0 < tin < tout and tout > 0.0:
        return tout - tin
    return 0.0


@njit(cache=True, fastmath=False)
def _intersect_triangle(ox, oy, oz, dx, dy, dz, v0x, v0y, v0z,
                        v1x, v1y, v1z, v2x, v2y, v2z):
    """Ray/triangle intersection parameter ``t`` (0 if no hit)."""
    kEpsilon = 1e-6

    a = v0x - v1x
    b = v0x - v2x
    c = dx
    d = v0x - ox

    e = v0y - v1y
    f = v0y - v2y
    g = dy
    h = v0y - oy

    i = v0z - v1z
    j = v0z - v2z
    k = dz
    l = v0z - oz

    m = f * k - g * j
    n = h * k - g * l
    p = f * l - h * j

    q = g * i - e * k
    s = e * j - f * i

    denom = a * m + b * q + c * s
    if denom == 0.0:
        inv_denom = 1e8
    else:
        inv_denom = 1.0 / denom

    e1 = d * m - b * n - c * p
    beta = e1 * inv_denom
    if beta < 0.0:
        return 0.0

    r = e * l - h * i
    e2 = a * n + d * q + c * r
    gamma = e2 * inv_denom
    if gamma < 0.0:
        return 0.0

    if beta + gamma > 1.0:
        return 0.0

    e3 = a * p - b * r + d * s
    t = e3 * inv_denom
    if t < kEpsilon:
        return 0.0

    return t


@njit(cache=True, fastmath=False)
def _intersect_polymesh(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez, faces):
    """Path length through a closed triangular mesh (first two hits)."""
    ox = ox / sizex
    oy = oy / sizey
    oz = oz / sizez
    dx = dx / sizex
    dy = dy / sizey
    dz = dz / sizez

    t0 = 0.0
    t1 = 0.0
    for kf in range(faces.shape[0]):
        t = _intersect_triangle(
            ox, oy, oz, dx, dy, dz,
            faces[kf, 0, 0], faces[kf, 0, 1], faces[kf, 0, 2],
            faces[kf, 1, 0], faces[kf, 1, 1], faces[kf, 1, 2],
            faces[kf, 2, 0], faces[kf, 2, 1], faces[kf, 2, 2])
        if t > 0.0:
            if t0 > 0.0:
                t1 = t
            else:
                t0 = t

    if t0 > 0.0 and t1 > 0.0:
        return abs(t0 - t1)
    return 0.0


@njit(cache=True, fastmath=False)
def _march_kernel(shape_code, faces, N, dx, dy, dz,
                  bbox_sizex, bbox_sizey, z_min, z_max,
                  scale_x, scale_y, scale_z, kEpsilon):
    """March the full N x N ray grid and return recorded path-length segments.

    ``shape_code``: 0=prism, 1=ellipsoid, 2=cylinder, 3=polymesh.  Mirrors the
    reference ``pathlengths`` marcher exactly, including the periodic wall
    cycling and the append pattern (one ``dr`` per march step plus one trailing
    ``dr`` per ray).
    """
    sx = bbox_sizex / N
    sy = bbox_sizey / N

    # Upper bound on recorded segments; grown-safe preallocation.  Each ray
    # records at least one trailing value plus one per z-slab crossing.
    out = np.empty(N * N * 64, dtype=np.float64)
    count = 0

    for j in range(N):
        for i in range(N):
            ox = -0.5 * bbox_sizex + (i + 0.5) * sx
            oy = -0.5 * bbox_sizey + (j + 0.5) * sy
            oz = z_min - kEpsilon

            ze = 0.0
            dr = 0.0
            while ze <= z_max:
                if shape_code == 0:
                    dr, _, _, _ = _intersect_bbox(
                        ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                elif shape_code == 1:
                    dr = _intersect_ellipsoid(
                        ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                elif shape_code == 2:
                    dr = _intersect_cylinder(
                        ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                else:
                    dr = _intersect_polymesh(
                        ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z, faces)

                _, xe, ye, ze = _intersect_bbox(
                    ox, oy, oz, dx, dy, dz, bbox_sizex, bbox_sizey, 1e6)

                if ze <= z_max:
                    if count >= out.shape[0]:
                        tmp = np.empty(out.shape[0] * 2, dtype=np.float64)
                        tmp[:count] = out[:count]
                        out = tmp
                    out[count] = dr
                    count += 1

                    ox = xe
                    oy = ye
                    oz = ze

                    if abs(ox - 0.5 * bbox_sizex) < kEpsilon:
                        ox = ox - bbox_sizex + kEpsilon
                    elif abs(ox + 0.5 * bbox_sizex) < kEpsilon:
                        ox = ox + bbox_sizex - kEpsilon

                    if abs(oy - 0.5 * bbox_sizey) < kEpsilon:
                        oy = oy - bbox_sizey + kEpsilon
                    elif abs(oy + 0.5 * bbox_sizey) < kEpsilon:
                        oy = oy + bbox_sizey - kEpsilon

            if count >= out.shape[0]:
                tmp = np.empty(out.shape[0] * 2, dtype=np.float64)
                tmp[:count] = out[:count]
                out = tmp
            out[count] = dr
            count += 1

    return out[:count]


@njit(cache=True, fastmath=False)
def _silhouette_shadow_kernel(shape_code, faces, N, dx, dy, dz,
                              bbox_sizex, bbox_sizey, z_launch,
                              scale_x, scale_y, scale_z):
    """Fraction of an N x N launch grid whose ray hits the crown (no periodicity).

    Rays are launched over a bounding box (which the caller enlarges so the
    tilted crown shadow fits), and each ray is counted at most once.  The single
    crown's horizontal shadow area is ``fraction * bbox_sizex * bbox_sizey`` and
    the beam-normal silhouette is that divided by ``cos(theta)``.
    """
    sx = bbox_sizex / N
    sy = bbox_sizey / N
    hits = 0
    for j in range(N):
        for i in range(N):
            ox = -0.5 * bbox_sizex + (i + 0.5) * sx
            oy = -0.5 * bbox_sizey + (j + 0.5) * sy
            oz = z_launch
            if shape_code == 0:
                dr, _, _, _ = _intersect_bbox(
                    ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
            elif shape_code == 1:
                dr = _intersect_ellipsoid(
                    ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
            elif shape_code == 2:
                dr = _intersect_cylinder(
                    ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
            else:
                dr = _intersect_polymesh(
                    ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z, faces)
            if dr > 1e-6:
                hits += 1
    return hits / (N * N)


# ===========================================================================
# Pure-Python reference implementations (retained for testing / clarity)
# ===========================================================================

def intersectBBox(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):

    # Intersection code below is adapted from Suffern (2007) Listing 19.1

    x0 = -0.5*sizex
    x1 = 0.5*sizex
    y0 = -0.5 * sizey
    y1 = 0.5 * sizey
    z0 = -1e-6
    z1 = sizez

    if dx == 0:
        a = 1e6
    else:
        a = 1.0 / dx
    if a >= 0:
        tx_min = (x0 - ox) * a
        tx_max = (x1 - ox) * a
    else:
        tx_min = (x1 - ox) * a
        tx_max = (x0 - ox) * a

    if dy == 0:
        b = 1e6
    else:
        b = 1.0 / dy
    if b >= 0:
        ty_min = (y0 - oy) * b
        ty_max = (y1 - oy) * b
    else:
        ty_min = (y1 - oy) * b
        ty_max = (y0 - oy) * b

    if dz == 0:
        c = 1e6
    else:
        c = 1.0 / dz
    if c >= 0:
        tz_min = (z0 - oz) * c
        tz_max = (z1 - oz) * c
    else:
        tz_min = (z1 - oz) * c
        tz_max = (z0 - oz) * c

    # find largest entering t value

    if tx_min > ty_min:
        t0 = tx_min
    else:
        t0 = ty_min

    if tz_min > t0:
        t0 = tz_min

    # find smallest exiting t value

    if tx_max < ty_max:
        t1 = tx_max
    else:
        t1 = ty_max

    if tz_max < t1:
        t1 = tz_max

    if t0 < t1 and t1 > 1e-6:
        if t0 > 1e-6:
            dr = t1-t0
        else:
            dr = t1
    else:
        dr = 0

    xe = ox + t1 * dx
    ye = oy + t1 * dy
    ze = oz + t1 * dz

    if dr == 0:
         raise Exception('Shouldnt be here')

    return dr, xe, ye, ze


def intersectEllipsoid(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):

    tempx = ox/sizex
    tempy = oy/sizey
    tempz = (oz - 0.5)/sizez

    dx = dx/sizex
    dy = dy/sizey
    dz = dz/sizez

    a = dx*dx + dy*dy + dz*dz

    b = 2.0 * (tempx*dx+tempy*dy+tempz*dz)

    c = (tempx*tempx+tempy*tempy+tempz*tempz) - 0.5*0.5
    disc = b * b - 4.0 * a * c

    if disc < 0.0:
        return 0
    else:
        e = sqrt(disc)
        denom = 2.0 * a

        t_small = (-b - e) / denom  # smaller root
        t_big = (-b + e) / denom  # larger root

        if t_small > 1e-6 or t_big > 1e-6:
            dr = abs(t_big-t_small)
            return dr

    return 0


def intersectCylinder(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez):

    tempx = ox / sizex
    tempy = oy / sizey
    tempz = oz / sizez

    dx = dx / sizex
    dy = dy / sizey
    dz = dz / sizez

    #19.5.3 of Suffern

    a = dx * dx + dy * dy

    b = 2.0 * (tempx * dx + tempy * dy)

    c = (tempx * tempx + tempy * tempy) - 0.5 * 0.5

    disc = b*b - 4 * a * c
    if disc >= 0 and a != 0:
        t0 = (-b + sqrt(disc)) / (2 * a)
        t1 = (-b - sqrt(disc)) / (2 * a)
    else:
        t0 = -1
        t1 = -1

    T = np.array([-1., -1., -1., -1.])

    # check if hit side surface of cylinder
    z0 = tempz + t0 * dz
    if t0 > 1e-6 and 1 >= z0 >= 0:
        T[0] = t0

    z1 = tempz + t1 * dz
    if t1 > 1e-6 and 1 >= z1 >= 0:
        T[1] = t1

    # check if it hits top of cylinder
    t_top = (1 - tempz) / dz
    hx = tempx + t_top * dx
    hy = tempy + t_top * dy
    if t_top > 0 and hx*hx + hy*hy <= 0.25:  # hits top
        T[2] = t_top

    # check if it hits bottom of cylinder
    t_bot = (0 - tempz) / dz
    hx = tempx + t_bot * dx
    hy = tempy + t_bot * dy
    if t_bot > 0 and hx*hx + hy*hy <= 0.25:  # hits bottom
        T[3] = t_bot

    t_int = T[T >= 0.0]
    if np.size(t_int) == 0:
        return 0
    Tin = np.min(t_int)
    Tout = np.max(T)

    if 0 < Tin < Tout and Tout > 0:
        dr = Tout-Tin
    else:
        dr = 0

    return dr


def importPLY(file):
    plydata = PlyData.read(file)
    return plydata


def intersectTriangle(ox, oy, oz, dx, dy, dz, vertices):

    kEpsilon = 1e-6

    a = vertices[0][0] - vertices[1][0]
    b = vertices[0][0] - vertices[2][0]
    c = dx
    d = vertices[0][0] - ox

    e = vertices[0][1] - vertices[1][1]
    f = vertices[0][1] - vertices[2][1]
    g = dy
    h = vertices[0][1] - oy

    i = vertices[0][2] - vertices[1][2]
    j = vertices[0][2] - vertices[2][2]
    k = dz
    l = vertices[0][2] - oz

    m = f * k - g * j
    n = h * k - g * l
    p = f * l - h * j

    q = g * i - e * k
    s = e * j - f * i

    if (a * m + b * q + c * s) == 0:
        inv_denom = 1e8
    else:
        inv_denom = 1.0 / (a * m + b * q + c * s)

    e1 = d * m - b * n - c * p
    beta = e1 * inv_denom

    if beta < 0.0:
        return 0

    r = e * l - h * i

    e2 = a * n + d * q + c * r

    gamma = e2 * inv_denom

    if gamma < 0.0:
        return 0

    if beta + gamma > 1.0:
        return 0

    e3 = a * p - b * r + d * s
    t = e3 * inv_denom

    if t < kEpsilon:
        return 0

    return t


def intersectPolymesh(ox, oy, oz, dx, dy, dz, sizex, sizey, sizez, plydata):

    ox = ox / sizex
    oy = oy / sizey
    oz = oz / sizez

    dx = dx / sizex
    dy = dy / sizey
    dz = dz / sizez

    vertices = plydata.elements[0].data
    faces = plydata.elements[1].data

    Nfaces = len(faces)

    face_verts = np.empty((3, 3))

    T = [0, 0]

    for face in range(0, Nfaces):
        f = faces[face][0]
        Nv = len(f)
        if Nv != 3:
            raise Exception('ERROR: only triangular elements are supported in PLY file geometry.')

        for v in range(0, 3):
            face_verts[0, v] = vertices[f[0]][v]
            face_verts[1, v] = vertices[f[1]][v]
            face_verts[2, v] = vertices[f[2]][v]

        t = intersectTriangle(ox, oy, oz, dx, dy, dz, face_verts)
        if t > 0:
            if T[0] > 0:
                T[1] = t
            else:
                T[0] = t

    if T[0] > 0 and T[1] > 0:
        return abs(T[0]-T[1])
    else:
        return 0


# ===========================================================================
# Mesh loading helper
# ===========================================================================

def _extract_faces(plydata):
    """Return a contiguous ``(Nfaces, 3, 3)`` float64 array of triangle verts.

    Raises if any face is not a triangle.
    """
    vertices = plydata.elements[0].data
    faces = plydata.elements[1].data

    verts = np.stack(
        [np.asarray(vertices['x'], dtype=np.float64),
         np.asarray(vertices['y'], dtype=np.float64),
         np.asarray(vertices['z'], dtype=np.float64)], axis=1)

    idx = np.empty((len(faces), 3), dtype=np.int64)
    for fi in range(len(faces)):
        f = faces[fi][0]
        if len(f) != 3:
            raise Exception(
                'ERROR: only triangular elements are supported in PLY file geometry.')
        idx[fi, 0] = f[0]
        idx[fi, 1] = f[1]
        idx[fi, 2] = f[2]

    return np.ascontiguousarray(verts[idx])


_SHAPE_CODE = {'prism': 0, 'ellipsoid': 1, 'cylinder': 2, 'polymesh': 3}


# ===========================================================================
# Public path-length API
# ===========================================================================

def pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith, ray_azimuth,
                nrays, plyfile='', outputfile='', degrees=False):
    """Compute the radiative path lengths through a shape for a beam direction.

    Launches a deterministic ``N x N`` grid of rays (``N = ceil(sqrt(nrays))``)
    from below the bounding box in the beam direction, using periodic side
    boundaries, and records the path length through the shape for each crossing.

    Parameters
    ----------
    shape : str
        'prism', 'ellipsoid' (alias 'sphere'), 'cylinder', or 'polymesh'
        (alias 'cone').  For a true sphere pass equal scales to 'ellipsoid'.
    scale_x, scale_y, scale_z : float
        Shape dimensions (m).  For ellipsoid these are the full axis lengths
        (semi-axes are scale/2); for cylinder scale_x/2 is the radius and
        scale_z the height; for prism they are the box side lengths.
    ray_zenith, ray_azimuth : float
        Beam zenith and azimuth.  Radians by default; degrees if ``degrees``.
    nrays : int
        Approximate number of rays; the actual grid is ``ceil(sqrt(nrays))**2``.
    plyfile : str, optional
        Path to a triangular-mesh PLY file (required for 'polymesh'/'cone').
    outputfile : str, optional
        If given, write all recorded path lengths to this file.
    degrees : bool, optional
        Interpret ``ray_zenith``/``ray_azimuth`` in degrees (default radians).

    Returns
    -------
    path_length : ndarray
        Path lengths of the rays that intersected the shape (m).
    projected_area : float
        Beam-normal silhouette area ``S(theta)`` of the crown (m^2).  At
        ``theta = 0`` this equals the true silhouette (e.g. ``pi R^2`` for a
        sphere).  At oblique angles the periodic side boundaries let a ray cross
        the tiled crown more than once, so this value is inflated relative to a
        single isolated crown; use :func:`silhouette_area` for the isolated
        ``S(theta)`` needed by the canopy binomial model.
    """

    shape = _normalize_shape(shape)
    ray_zenith = _to_radians(ray_zenith, degrees)
    ray_azimuth = _to_radians(ray_azimuth, degrees)

    kEpsilon = 1e-5

    N = int(ceil(sqrt(nrays)))

    # Ray direction Cartesian unit vector
    dx = sin(ray_zenith) * cos(ray_azimuth)
    dy = sin(ray_zenith) * sin(ray_azimuth)
    dz = cos(ray_zenith)

    faces = np.empty((0, 3, 3), dtype=np.float64)

    if shape == 'polymesh':
        if len(plyfile) == 0:
            raise Exception('Path to PLY file must be provided for polymesh intersection.')
        elif not os.path.exists(plyfile):
            raise Exception('PLY file does not exist.')
        plydata = PlyData.read(plyfile)
        faces = _extract_faces(plydata)

        vertices = plydata.elements[0].data
        Nvertices = len(vertices)
        bx_min = 1e6
        bx_max = -1e6
        by_min = 1e6
        by_max = -1e6
        z_min = 1e6
        z_max = -1e6
        for vert in range(0, Nvertices):
            vx = vertices[vert][0] * scale_x
            vy = vertices[vert][1] * scale_y
            vz = vertices[vert][2] * scale_z
            if vx < bx_min:
                bx_min = vx
            if vx > bx_max:
                bx_max = vx
            if vy < by_min:
                by_min = vy
            if vy > by_max:
                by_max = vy
            if vz < z_min:
                z_min = vz
            if vz > z_max:
                z_max = vz
        bbox_sizex = 2*max(abs(bx_max), abs(bx_min)) * (1.0 + kEpsilon)
        bbox_sizey = 2*max(abs(by_max), abs(by_min)) * (1.0 + kEpsilon)
        z_min = z_min
        z_max = z_max * (1.0 + kEpsilon)
    else:
        bbox_sizex = scale_x * (1.0 + kEpsilon)
        bbox_sizey = scale_y * (1.0 + kEpsilon)
        z_min = 0
        z_max = scale_z * (1.0 + kEpsilon)

    path_length = _march_kernel(
        _SHAPE_CODE[shape], faces, N, dx, dy, dz,
        float(bbox_sizex), float(bbox_sizey), float(z_min), float(z_max),
        float(scale_x), float(scale_y), float(scale_z), kEpsilon)

    path_length = np.asarray(path_length, dtype=float)

    projected_area = np.sum(path_length > kEpsilon) / (N * N) \
        * bbox_sizex * bbox_sizey * cos(ray_zenith)

    if outputfile != '':
        np.savetxt(outputfile, path_length, delimiter=',')

    return path_length[path_length > kEpsilon], projected_area


def pathlengthdistribution(shape, scale_x, scale_y, scale_z, ray_zenith,
                           ray_azimuth, nrays, plyfile='', bins=10,
                           normalize=True, degrees=False):
    """Probability density (or histogram) of path lengths through a shape."""

    path_lengths, _ = pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith,
                                  ray_azimuth, nrays, plyfile, degrees=degrees)

    hist, bin_edges = np.histogram(path_lengths, bins=bins, density=normalize)

    return hist, bin_edges


# ===========================================================================
# Crown-level interception and absorption
# ===========================================================================

def crown_interception(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                       ray_zenith, ray_azimuth, nrays,
                       path_multiplier=1.0, absorptivity=1.0,
                       plyfile='', degrees=False):
    r"""Per-crown probability of intercepting a leaf (Bailey et al. 2020).

    Monte-Carlo estimate of

    .. math::
        P_{\mathrm{leaf}}(\theta) = \int p(r|\theta)\,
            [1 - e^{-m\,\zeta\,G\,a\,r}]\,dr
        \approx \frac{1}{M}\sum_i [1 - e^{-m\,\zeta\,G\,a\,r_i}]

    over the ``M`` path lengths that intersect the crown, where ``m`` is
    ``path_multiplier`` and ``zeta`` is ``absorptivity``.

    Parameters
    ----------
    Gtheta : float
        Fraction of leaf area projected in the beam direction (Ross G).
    LAD : float
        Leaf area density ``a`` (m^2 m^-3).
    path_multiplier : float, optional
        Path-length multiplier ``m``.  Use 1.0 for the direct (single-pass)
        term and 2.0 for the mode-2 scattering approximation of Ponce de Leon
        et al. (2025).
    absorptivity : float, optional
        Leaf absorptivity ``zeta = 1 - rho_l - tau_l``.  Use 1.0 for total
        interception (no scattering).

    Returns
    -------
    float
        Interception probability in [0, 1].
    """
    path_length, _ = pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith,
                                 ray_azimuth, nrays, plyfile, degrees=degrees)
    if path_length.size == 0:
        return 0.0
    k = path_multiplier * absorptivity * Gtheta * LAD
    return float(np.mean(1.0 - np.exp(-k * path_length)))


def crownabsorptionfraction(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                            ray_zenith, ray_azimuth, nrays, plyfile='',
                            degrees=False):
    """Per-crown interception probability P_leaf (Bailey et al. 2020).

    This is the corrected form of the crown-scale interception probability and
    is equivalent to :func:`crown_interception` with unit multiplier and
    absorptivity.  (Earlier versions of this function erroneously divided by the
    projected area; that has been fixed so the return value is a probability in
    [0, 1].)
    """
    return crown_interception(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                              ray_zenith, ray_azimuth, nrays, plyfile=plyfile,
                              degrees=degrees)


def crown_sunlit_fraction(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                          ray_zenith, ray_azimuth, nrays, plyfile='',
                          degrees=False):
    r"""Fraction of crown leaf area that is directly sunlit (Bailey et al. 2020).

    A leaf element a distance ``l`` from the point where the beam enters the
    crown is sunlit (unshaded by leaves above it along the beam) with
    probability ``exp(-G a l)``.  Averaging that over the leaf area encountered
    along a beam of length ``r`` gives a per-ray sunlit fraction

    .. math::
        f_{\mathrm{sun,ray}} = \frac{1}{r}\int_0^r e^{-G a l}\,dl
        = \frac{1 - e^{-G a r}}{G a r},

    and the crown sunlit fraction is the path-length-weighted mean over rays,

    .. math::
        f_{\mathrm{sun}} = \frac{\sum_i (1 - e^{-G a r_i})/(G a)}{\sum_i r_i}.

    Limits: ``f_sun -> 1`` for a thin crown (``G a -> 0``) and ``f_sun -> 0``
    for a thick crown.
    """
    path_length, _ = pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith,
                                 ray_azimuth, nrays, plyfile, degrees=degrees)
    if path_length.size == 0:
        return 0.0
    Ga = Gtheta * LAD
    if Ga <= 0.0:
        return 1.0
    numer = np.sum((1.0 - np.exp(-Ga * path_length)) / Ga)
    denom = np.sum(path_length)
    if denom <= 0.0:
        return 0.0
    return float(numer / denom)


def crownsunlitfraction(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                        ray_zenith, ray_azimuth, nrays, plyfile='',
                        degrees=False):
    """Backward-compatible alias of :func:`crown_sunlit_fraction`.

    (Previously this was an accidental duplicate of the absorption function;
    it now computes the sunlit leaf-area fraction as intended.)
    """
    return crown_sunlit_fraction(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                                 ray_zenith, ray_azimuth, nrays, plyfile=plyfile,
                                 degrees=degrees)


# ===========================================================================
# Canopy-scale binomial interception, diffuse, and scattering
# ===========================================================================

def silhouette_area(shape, scale_x, scale_y, scale_z, ray_zenith,
                    ray_azimuth=0.0, nrays=10000, plyfile='', degrees=False):
    r"""Crown horizontal ground-shadow area ``S(theta, phi)`` (Bailey et al. 2020;
    Ponce de Leon et al. 2026, Table 2).

    This is the horizontal shadow the crown casts, matching the ``S(theta)`` used
    by the canopy binomial model's ``Nc = S(theta)/S(0)`` -- **not** the
    beam-normal silhouette (the two differ by a factor ``cos(theta)``).

    The shadow depends on the beam **azimuth** ``ray_azimuth`` whenever the crown
    is not horizontally circular (triaxial ellipsoid, elliptical cylinder,
    prism, or a general mesh).  ``ray_azimuth`` is measured from the ``+x`` axis
    in the same crown frame as the ray direction used by :func:`pathlengths` /
    :func:`crown_interception` (``dx = sin(theta) cos(azimuth)``,
    ``dy = sin(theta) sin(azimuth)``), so ``S(theta, phi)`` and the per-crown
    ``P_leaf`` are evaluated for a consistent beam direction.  At ``theta = 0``
    the shadow is azimuth-independent (``S(0) = `` the horizontal cross-section).

    Analytic where possible (semi-axes ``Rx = scale_x/2``, ``Ry = scale_y/2``,
    ``Rz = scale_z/2``; beam unit vector ``d = (sin th cos phi, sin th sin phi,
    cos th)``):

    * ellipsoid / sphere: the beam-normal silhouette of an ellipsoid is
      ``pi sqrt(Rx^2 Ry^2 dz^2 + Ry^2 Rz^2 dx^2 + Rx^2 Rz^2 dy^2)`` and the
      ground shadow is that divided by ``cos(theta) = dz``.  Reduces to
      ``pi R^2/cos(theta)`` for a sphere and, for ``Rx == Ry == R``, ``Rz == V``,
      to the Ponce de Leon et al. (2026) Table 2 form
      ``pi R^2 sqrt(1 + (V/R)^2 tan^2 theta)`` (azimuth-independent).
    * cylinder: ``pi Rx Ry + w(phi) H tan(theta)`` -- the horizontal top ellipse
      plus the side-wall band, whose width ``w(phi) = 2 sqrt(Rx^2 sin^2 phi +
      Ry^2 cos^2 phi)`` is the cross-section extent perpendicular to the beam's
      horizontal direction (``H = scale_z``).  Reduces to ``pi R^2 + 2 R H
      tan(theta)`` for a circular cylinder.
    * prism / polymesh: computed numerically from a single (non-periodic)
      shadow trace over an enlarged bounding box; the horizontal ground shadow
      of one isolated crown for the given ``(theta, phi)``.

    Notes
    -----
    The numeric route must **not** use :func:`pathlengths`' ``projected_area``:
    that sampler uses periodic side boundaries, so at oblique zenith angles a
    ray crosses the tiled crown more than once and the projected area is
    inflated.  A dedicated non-periodic first-hit trace is used instead.

    Returns
    -------
    float
        Silhouette area (m^2).
    """
    shape = _normalize_shape(shape)
    theta = _to_radians(ray_zenith, degrees)
    azimuth = _to_radians(ray_azimuth, degrees)

    if shape == 'ellipsoid':
        # Horizontal ground-shadow area of a triaxial ellipsoid with semi-axes
        # (Rx, Ry, Rz) for a beam at zenith ``theta`` and azimuth ``azimuth``.
        # The beam-normal projection of an ellipsoid along a unit direction
        # d = (dx, dy, dz) is
        #     pi sqrt(Rx^2 Ry^2 dz^2 + Ry^2 Rz^2 dx^2 + Rx^2 Rz^2 dy^2)
        # and the ground shadow is that divided by cos(theta) = dz.  For
        # Rx == Ry this is azimuth-independent and reduces to the Ponce de Leon
        # et al. (2026) Table 2 form; for a sphere it is pi R^2/cos(theta).
        Rx = 0.5 * scale_x
        Ry = 0.5 * scale_y
        Rz = 0.5 * scale_z
        dxh = sin(theta) * cos(azimuth)
        dyh = sin(theta) * sin(azimuth)
        dzh = cos(theta)
        perp = pi * sqrt(Rx ** 2 * Ry ** 2 * dzh ** 2 +
                         Ry ** 2 * Rz ** 2 * dxh ** 2 +
                         Rx ** 2 * Rz ** 2 * dyh ** 2)
        return perp / dzh
    if shape == 'cylinder':
        # (Elliptical) cylinder ground shadow: horizontal top ellipse pi Rx Ry
        # (constant area) plus a side-wall rectangle of length H tan(theta) and
        # width equal to the cross-section extent perpendicular to the beam's
        # horizontal direction, w(phi) = 2 sqrt(Rx^2 sin^2 phi + Ry^2 cos^2 phi).
        Rx = 0.5 * scale_x
        Ry = 0.5 * scale_y
        H = scale_z
        width = 2.0 * sqrt(Rx ** 2 * sin(azimuth) ** 2 +
                           Ry ** 2 * cos(azimuth) ** 2)
        return pi * Rx * Ry + width * H * np.tan(theta)

    # prism / polymesh: non-periodic shadow trace over an enlarged box.
    faces = np.empty((0, 3, 3), dtype=np.float64)
    if shape == 'polymesh':
        if len(plyfile) == 0:
            raise Exception('Path to PLY file must be provided for polymesh silhouette.')
        plydata = PlyData.read(plyfile)
        faces = _extract_faces(plydata)
        verts = plydata.elements[0].data
        vx = np.abs(np.asarray(verts['x'], dtype=float)) * scale_x
        vy = np.abs(np.asarray(verts['y'], dtype=float)) * scale_y
        vz = np.asarray(verts['z'], dtype=float) * scale_z
        halfx, halfy = vx.max(), vy.max()
        z_top = vz.max()
    else:  # prism
        halfx, halfy = 0.5 * scale_x, 0.5 * scale_y
        z_top = scale_z

    # Enlarge the launch box so the tilted shadow (shifted by up to
    # z_top*tan(theta) in the beam's horizontal direction) fits in both x and y
    # for any azimuth; the extra empty area cancels because we divide by the box
    # area to get the hit fraction.
    shift = z_top * np.tan(theta)
    bbox_sizex = 2.0 * (halfx + abs(shift) + 0.01) * 1.001
    bbox_sizey = 2.0 * (halfy + abs(shift) + 0.01) * 1.001

    dx = sin(theta) * cos(azimuth)
    dy = sin(theta) * sin(azimuth)
    dz = cos(theta)

    N = int(ceil(sqrt(nrays)))
    frac = _silhouette_shadow_kernel(
        _SHAPE_CODE[shape], faces, N, dx, dy, dz,
        float(bbox_sizex), float(bbox_sizey), -1e-6,
        float(scale_x), float(scale_y), float(scale_z))
    # Rays launched over a horizontal grid along the beam direction: the set of
    # origins whose ray hits the crown is the crown shadow cast onto the launch
    # plane.  Its area is the horizontal ground shadow S(theta) directly (the
    # paper's Table 2 convention), so return it without a cos(theta) factor.
    horizontal_shadow = frac * bbox_sizex * bbox_sizey
    return horizontal_shadow


def crown_volume(shape, scale_x, scale_y, scale_z, plyfile=''):
    r"""Volume of a single crown (m^3).

    Analytic for the primitive shapes; computed from the mesh (divergence
    theorem) for ``polymesh``/``cone``.  The ``scale_*`` arguments are the full
    crown extents in each axis, matching the convention used elsewhere in this
    module (e.g. an ellipsoid has horizontal radius ``scale_x/2``).

    * ellipsoid: ``(4/3) pi (scale_x/2)(scale_y/2)(scale_z/2)``
    * cylinder: elliptical cross-section ``pi (scale_x/2)(scale_y/2)`` times
      height ``scale_z``
    * prism / box: ``scale_x scale_y scale_z``
    * polymesh / cone: signed tetrahedron sum over the (scaled) triangle faces

    Returns
    -------
    float
        Crown volume (m^3).
    """
    shape = _normalize_shape(shape)

    if shape == 'ellipsoid':
        return float((4.0 / 3.0) * pi * (0.5 * scale_x) *
                     (0.5 * scale_y) * (0.5 * scale_z))
    if shape == 'cylinder':
        return float(pi * (0.5 * scale_x) * (0.5 * scale_y) * scale_z)
    if shape == 'prism':
        return float(scale_x * scale_y * scale_z)

    # polymesh / cone: sum signed volumes of tetrahedra formed by each triangle
    # face and the origin (divergence theorem).  Vertices are scaled per axis
    # to match the pathlengths/silhouette convention (vx = x * scale_x).
    if len(plyfile) == 0:
        raise Exception('Path to PLY file must be provided for polymesh volume.')
    plydata = PlyData.read(plyfile)
    faces = _extract_faces(plydata)  # (Nfaces, 3, 3)
    scale = np.array([scale_x, scale_y, scale_z], dtype=np.float64)
    v = faces * scale  # broadcast over the last axis
    v0, v1, v2 = v[:, 0, :], v[:, 1, :], v[:, 2, :]
    signed_six = np.einsum('ij,ij->i', v0, np.cross(v1, v2))
    return float(abs(signed_six.sum()) / 6.0)


def _crown_perp_width(shape, scale_x, scale_y, azimuth, plyfile=''):
    r"""Full crown width perpendicular to the beam's horizontal azimuth (m).

    This is the crown's cross-beam extent -- the ``w_perp`` used to set the
    canopy interception ceiling.  ``azimuth`` is the beam azimuth in the crown
    frame (radians, measured from +x, same convention as :func:`pathlengths`).
    """
    sa, ca = abs(sin(azimuth)), abs(cos(azimuth))
    if shape == 'prism':
        # Support width of the box footprint perpendicular to the beam.
        return scale_x * sa + scale_y * ca
    if shape in ('ellipsoid', 'cylinder'):
        # Support width of the (elliptical) horizontal cross-section.
        Rx, Ry = 0.5 * scale_x, 0.5 * scale_y
        return 2.0 * sqrt(Rx * Rx * sa * sa + Ry * Ry * ca * ca)
    # polymesh: perpendicular extent of the scaled vertices onto n=(-sin,cos).
    plydata = PlyData.read(plyfile)
    verts = plydata.elements[0].data
    vx = np.asarray(verts['x'], dtype=float) * scale_x
    vy = np.asarray(verts['y'], dtype=float) * scale_y
    proj = -vx * sin(azimuth) + vy * cos(azimuth)
    return float(proj.max() - proj.min())


def canopy_interception(Gtheta, LAD, shape, scale_x, scale_y, scale_z,
                        ray_zenith, ray_azimuth, nrays, sr, sp, phi=None,
                        path_multiplier=1.0, absorptivity=1.0,
                        plyfile='', degrees=False, P_leaf=None):
    r"""Canopy-level binomial interception probability (closure-corrected).

    .. math::
        N_c &= S(\theta,\phi) / S(0) \\
        C &= \min\!\left(1,\; w_\perp / s_\perp\right),\quad
            s_\perp = s_r \cos^2\phi + s_p \sin^2\phi \\
        P_c &= C\left[1 - \left(1 - \frac{S(0)}{C\, s_r s_p}
            P_{\mathrm{leaf}}\right)^{N_c}\right]

    This is the Bailey et al. (2020) binomial (Eq. 13) with a corrected
    saturation ceiling.  The published form pre-multiplies by ``s^2/(s_r s_p)``
    with ``s = s_r sin^2(phi) + s_p cos^2(phi)``, which pins ``P_c`` at an
    azimuth-dependent geometric ceiling ``s_p/s_r`` (at ``phi=0``); against an
    independent 3-D reference that ceiling is too low for a *closed* canopy
    (crowns that fill the cell should intercept ~all light at grazing) and lets
    ``P_c`` overshoot when ``s_r < s_p``.

    Here the ceiling is instead the fraction of the beam cross-section the row
    geometry can block -- the crown cross-beam width ``w_perp`` over the
    perpendicular spacing ``s_perp`` -- and the footprint normalization is
    carried inside the per-layer cover so the ceiling ``C`` does not disturb the
    (validated) sparse limit ``P_c -> S(theta) P_leaf/(s_r s_p)``.  The nadir
    limit (``N_c=1``) stays exact at ``S(0) P_leaf/(s_r s_p)``, and the model
    reduces to Eq. 13 in the open-canopy (sparse) regime.  Row-induced
    azimuthal anisotropy is retained -- it enters through ``C`` (via ``w_perp``,
    ``s_perp``) and through the azimuth-aware ``S(theta, phi)``.

    Parameters
    ----------
    sr, sp : float
        Row spacing and plant spacing (m).
    phi : float, optional
        Beam azimuth relative to the row direction.  Defaults to
        ``ray_azimuth`` (interpreted with the same ``degrees`` convention).
    P_leaf : float, optional
        Precomputed per-crown interception probability; if given, the crown
        interception is not recomputed (useful to share work across modes).

    Returns
    -------
    float
        Canopy interception probability in [0, 1].
    """
    shape = _normalize_shape(shape)
    theta = _to_radians(ray_zenith, degrees)
    azimuth = _to_radians(ray_azimuth, degrees)
    phi = azimuth if phi is None else _to_radians(phi, degrees)

    if P_leaf is None:
        P_leaf = crown_interception(
            Gtheta, LAD, shape, scale_x, scale_y, scale_z, theta, azimuth,
            nrays, path_multiplier=path_multiplier, absorptivity=absorptivity,
            plyfile=plyfile)

    # S(0) is azimuth-independent (theta = 0); S(theta) uses the beam azimuth
    # in the same crown frame as the per-crown ray tracing above, so Nc is the
    # ground-shadow ratio for the actual beam direction (matters for prisms and
    # triaxial ellipsoids; a no-op for horizontally circular crowns).
    S0 = silhouette_area(shape, scale_x, scale_y, scale_z, 0.0, nrays=nrays,
                         plyfile=plyfile)
    Stheta = silhouette_area(shape, scale_x, scale_y, scale_z, theta, azimuth,
                             nrays=nrays, plyfile=plyfile)
    Nc = Stheta / S0

    # Closure-corrected ceiling: the maximum interception the row geometry can
    # produce for this beam is the crown cross-beam width w_perp over the
    # perpendicular spacing s_perp (both azimuth-dependent).  w_perp is a crown
    # property (beam azimuth in the crown frame); s_perp is a lattice property
    # (beam azimuth relative to the rows).
    w_perp = _crown_perp_width(shape, scale_x, scale_y, azimuth, plyfile)
    s_perp = sr * cos(phi) ** 2 + sp * sin(phi) ** 2
    C = min(1.0, w_perp / s_perp) if s_perp > 0.0 else 1.0
    if C <= 0.0:
        return 0.0

    # Per-layer cover in the ceiling-scaled cell, pinned by the sparse limit so
    # C does not perturb it; compounded over Nc crown layers, then scaled to C.
    # At Nc=1 this is the exact nadir cover S(0) P_leaf/(sr sp); as Nc grows the
    # bracket saturates and P_c -> C.
    base = max(0.0, 1.0 - (S0 / (C * sr * sp)) * P_leaf)
    Pc = C * (1.0 - base ** Nc)
    return float(min(max(Pc, 0.0), 1.0))


def diffuse_interception(Gtheta, LAD, shape, scale_x, scale_y, scale_z, nrays,
                         sr=None, sp=None, fd=None, n_zenith=18, n_azimuth=1,
                         level='canopy', path_multiplier=1.0, absorptivity=1.0,
                         plyfile='', degrees=False):
    r"""Hemispherically-integrated (diffuse) interception (Bailey et al. 2020).

    Azimuthally-symmetric reduced form (``n_azimuth == 1``):

    .. math::
        P_{\mathrm{diff}} = 2 \int_0^{\pi/2}
            f_d(\theta)\, P(\theta)\, \cos\theta\, \sin\theta\, d\theta,
        \qquad \int_0^{\pi/2} f_d(\theta)\, d\theta = \pi/2.

    Full 2-D form (``n_azimuth > 1``), for row canopies that are not
    azimuthally symmetric:

    .. math::
        P_{\mathrm{diff}} = \frac{1}{\pi}
            \int_0^{2\pi}\!\!\int_0^{\pi/2}
            f_d(\theta,\phi)\, P(\theta,\phi)\, \cos\theta\, \sin\theta\,
            d\theta\, d\phi.

    Parameters
    ----------
    fd : callable, optional
        Sky-radiance weighting ``f_d(theta)`` or ``f_d(theta, phi)`` (radians).
        ``None`` gives an isotropic sky (``f_d == 1``).
    n_zenith, n_azimuth : int
        Number of Gauss-Legendre zenith nodes and (for the 2-D form) uniform
        azimuth nodes.
    level : {'canopy', 'crown'}
        Use :func:`canopy_interception` (requires ``sr``/``sp``) or
        :func:`crown_interception` as the directional interception ``P``.

    Returns
    -------
    float
        Diffuse interception fraction in [0, 1].
    """
    shape = _normalize_shape(shape)

    if level == 'canopy' and (sr is None or sp is None):
        raise ValueError("level='canopy' requires sr and sp.")

    # Gauss-Legendre nodes on [0, pi/2] for the zenith integral.
    x, w = np.polynomial.legendre.leggauss(n_zenith)
    a, b = 0.0, 0.5 * pi
    thetas = 0.5 * (b - a) * x + 0.5 * (b + a)
    wz = 0.5 * (b - a) * w

    def P_of(theta, phi):
        if level == 'crown':
            return crown_interception(
                Gtheta, LAD, shape, scale_x, scale_y, scale_z, theta, phi,
                nrays, path_multiplier=path_multiplier,
                absorptivity=absorptivity, plyfile=plyfile)
        return canopy_interception(
            Gtheta, LAD, shape, scale_x, scale_y, scale_z, theta, phi, nrays,
            sr, sp, phi=phi, path_multiplier=path_multiplier,
            absorptivity=absorptivity, plyfile=plyfile)

    if n_azimuth <= 1:
        total = 0.0
        for th, wth in zip(thetas, wz):
            fdv = 1.0 if fd is None else float(fd(th))
            total += wth * fdv * P_of(th, 0.0) * cos(th) * sin(th)
        return float(2.0 * total)

    phis = (np.arange(n_azimuth) + 0.5) * (2.0 * pi / n_azimuth)
    dphi = 2.0 * pi / n_azimuth
    total = 0.0
    for th, wth in zip(thetas, wz):
        for ph in phis:
            fdv = 1.0 if fd is None else float(fd(th, ph))
            total += wth * dphi * fdv * P_of(th, ph) * cos(th) * sin(th)
    return float(total / pi)


def absorbed_fraction(LAD, Gtheta, shape, scale_x, scale_y, scale_z,
                      ray_zenith, ray_azimuth, nrays, sr, sp,
                      rho_l, tau_l, rho_s, Q0=1.0, phi=None,
                      plyfile='', degrees=False):
    r"""Absorbed shortwave fraction via three-mode scattering (Ponce de Leon 2025).

    .. math::
        \zeta &= 1 - \rho_l - \tau_l \\
        Q/Q_0 &= P_{c2}
            + (1 - P_{c1})\, \rho_s\, \frac{S(0)}{s^2}\, P_{f1}

    where mode 1 is direct first-pass absorption (single path, absorptivity
    ``zeta``), mode 2 approximates scattered-then-absorbed radiation by doubling
    the path length, and mode 3 is ground-reflected radiation re-intercepted by
    the canopy.

    Parameters
    ----------
    rho_l, tau_l : float
        Leaf reflectivity and transmissivity for the band.
    rho_s : float
        Ground (soil) reflectivity.
    Q0 : float, optional
        Incident flux (default 1.0, so the return value is a fraction).

    Returns
    -------
    float
        Absorbed flux ``Q`` (equals the absorbed fraction when ``Q0 == 1``).
    """
    shape = _normalize_shape(shape)
    theta = _to_radians(ray_zenith, degrees)
    azimuth = _to_radians(ray_azimuth, degrees)
    phi = azimuth if phi is None else _to_radians(phi, degrees)

    zeta = 1.0 - rho_l - tau_l

    # Per-crown interception probabilities (share ray tracing across modes by
    # computing path lengths once).
    path_length, _ = pathlengths(shape, scale_x, scale_y, scale_z, theta,
                                 azimuth, nrays, plyfile)
    if path_length.size == 0:
        return 0.0

    def Pleaf(mult):
        k = mult * zeta * Gtheta * LAD
        return float(np.mean(1.0 - np.exp(-k * path_length)))

    P_f1 = Pleaf(1.0)
    P_f2 = Pleaf(2.0)

    P_c1 = canopy_interception(
        Gtheta, LAD, shape, scale_x, scale_y, scale_z, theta, azimuth, nrays,
        sr, sp, phi=phi, plyfile=plyfile, P_leaf=P_f1)
    P_c2 = canopy_interception(
        Gtheta, LAD, shape, scale_x, scale_y, scale_z, theta, azimuth, nrays,
        sr, sp, phi=phi, plyfile=plyfile, P_leaf=P_f2)

    S0 = silhouette_area(shape, scale_x, scale_y, scale_z, 0.0, nrays=nrays,
                         plyfile=plyfile)
    s = sr * sin(phi) ** 2 + sp * cos(phi) ** 2

    absorbed = P_c2 + (1.0 - P_c1) * rho_s * (S0 / s ** 2) * P_f1
    return float(Q0 * absorbed)
