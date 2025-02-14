import os
import numpy as np
from numpy import sqrt, sin, arcsin, cos, arccos, exp, pi, linspace, ceil
from plyfile import PlyData, PlyElement


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

    # xe = ox + t * dx
    # ye = oy + t * dy
    # ze = oz + t * dz

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

    Nvertices = len(vertices)
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

        # THIS CODE IS FOR POLYGONS WITH MORE THAN 3 VERTICES. IT DOESN'T SEEM TO WORK
        # for tri in range(0, Nv-2):
        #     for v in range(0, 3):
        #         face_verts[0, v] = vertices[f[0]][v]
        #         face_verts[1, v] = vertices[f[1+tri]][v]
        #         face_verts[2, v] = vertices[f[2+tri]][v]

        # print(face_verts[0, :])
        # ax.plot(face_verts[:, 0], face_verts[:, 1], face_verts[:, 2], '-')
        # if tri > 0:
        #     ax.plot_trisurf(face_verts[:, 0], face_verts[:, 1], face_verts[:, 2])
        # for i in range(0, 3):
        #     ax.plot(face_verts[i][0], face_verts[i][1], face_verts[i][2], '.')

        # plt.show()

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


def pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith, ray_azimuth, nrays, plyfile='', outputfile=''):

    kEpsilon = 1e-5

    N = int(ceil(sqrt(nrays)))

    # Ray direction Cartesian unit vector
    dx = sin(ray_zenith) * cos(ray_azimuth)
    dy = sin(ray_zenith) * sin(ray_azimuth)
    dz = cos(ray_zenith)

    path_length = np.zeros(N*N)

    plydata=[]
    if shape == 'polymesh':
        if len(plyfile) == 0:
            raise Exception('Path to PLY file must be provided for polymesh intersection.')
        # check that plyfile exists
        elif not os.path.exists(plyfile):
            raise Exception('PLY file does not exist.')
        plydata = PlyData.read(plyfile)

    if shape == 'polymesh':
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

    sx = bbox_sizex/N
    sy = bbox_sizey/N

    # loop over all rays, which originate at the bottom of the box
    for j in range(0, N):
        for i in range(0, N):

            # ray origin point
            ox = -0.5*bbox_sizex + (i+0.5)*sx
            oy = -0.5*bbox_sizey + (j+0.5)*sy
            oz = z_min-kEpsilon

            ze = 0
            dr = 0
            while ze <= z_max:

                # Intersect shape
                if shape == 'prism':
                    dr, _, _, _, = intersectBBox(ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                elif shape == 'ellipsoid':
                    dr = intersectEllipsoid(ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                elif shape == 'cylinder':
                    dr = intersectCylinder(ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z)
                elif shape == 'polymesh':
                    dr = intersectPolymesh(ox, oy, oz, dx, dy, dz, scale_x, scale_y, scale_z, plydata)
                else:
                    raise Exception('Invalid shape argument.')

                # Intersect bounding box walls
                _, xe, ye, ze = intersectBBox(ox, oy, oz, dx, dy, dz, bbox_sizex, bbox_sizey, 1e6)

                if ze <= z_max:  # intersection below object height -> record path length and periodically cycle ray

                    path_length = np.append(path_length, dr)

                    ox = xe
                    oy = ye
                    oz = ze

                    if abs(ox-0.5*bbox_sizex) < kEpsilon:  # hit +x wall
                        ox = ox - bbox_sizex + kEpsilon
                    elif abs(ox+0.5*bbox_sizex) < kEpsilon:  # hit -x wall
                        ox = ox + bbox_sizex - kEpsilon

                    if abs(oy-0.5*bbox_sizey) < kEpsilon:  # hit +y wall
                        oy = oy - bbox_sizey + kEpsilon
                    elif abs(oy + 0.5 * bbox_sizey) < kEpsilon:  # hit -y wall
                        oy = oy + bbox_sizey - kEpsilon

            path_length[i+j*N] = dr


    if( outputfile != '' ):
        np.savetxt(outputfile, path_length, delimiter=',')

    return path_length[path_length > kEpsilon]

def pathlengthdistribution(shape, scale_x, scale_y, scale_z, ray_zenith, ray_azimuth, nrays, plyfile='', bins=10, normalize=True):

    path_lengths = pathlengths(shape, scale_x, scale_y, scale_z, ray_zenith, ray_azimuth, nrays, plyfile)

    hist, bin_edges = np.histogram(path_lengths, bins=bins, density=normalize)

    return hist, bin_edges
