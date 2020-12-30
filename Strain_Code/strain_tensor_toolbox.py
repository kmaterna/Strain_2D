# A toolbox once we have strain tensor components. 
# We want to compute:
# Rotation vs. Volume Strain
# 2nd invariant
# Dilatation
# Max Shear Strain
# Eigenvectors and Eigenvalues


import numpy as np
import math as m


def second_invariant(exx, exy, eyy):
    e2nd = exx * eyy - exy * exy;
    return e2nd;


def eigenvector_eigenvalue(exx, exy, eyy):
    T = np.array([[exx, exy], [exy, eyy]]);  # the tensor
    w, v = np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
    return [w[0], w[1], v];


def max_shear_strain(exx, exy, eyy):
    T = np.array([[exx, exy], [exy, eyy]]);  # the tensor
    w, v = np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
    # w = eigenvalues
    # v = eigenvectors
    max_shear = (w[0] - w[1]) * 0.5;
    return max_shear;


def compute_displacement_gradients(up, vp, ur, vr, uq, vq, dx, dy):
    # up, vp describe velocity at a reference point P
    # R and Q are two other points: Q offset by dx in the x direction, and R offset by dy in the y direction.
    # In practical usage, these are in mm/yr and km.
    dudx = (uq - up) / dx;
    dvdx = (vq - vp) / dx;
    dudy = (ur - up) / dy;
    dvdy = (vr - vp) / dy;
    return [dudx, dvdx, dudy, dvdy];


def compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy):
    # Given a displacement tensor, compute the relevant parts of the strain and rotation tensors.
    # Also converts to nanostrain per year.
    # Rot is the off-diagonal element of the rotation tensor
    # http://www.engr.colostate.edu/~thompson/hPage/CourseMat/Tutorials/Solid_Mechanics/rotations.pdf
    exx = dudx * 1000;
    exy = (0.5 * (dvdx + dudy)) * 1000;
    eyy = dvdy * 1000;
    rot = (0.5 * (dvdx - dudy));
    rot = rot * 1000.0;
    return [exx, exy, eyy, rot];


def max_shortening_azimuth(e1, e2, v00, v01, v10, v11):
    az = np.zeros(e1.shape)
    for i in range(len(e1)):
        for j in range(len(e1[0])):
            az[i][j] = azimuth_math(e1[i][j], e2[i][j], v00[i][j], v01[i][j], v10[i][j], v11[i][j]);
            if az[i][j] > 179.999:
                print("OOPSSSSS %d %d" % (i, j))
    print("Minimum azimuth: %.3f degrees" % np.min(az))
    print("Maximum azimuth: %.3f degrees" % np.max(az))
    return az


def azimuth_math(e1, e2, v00, v01, v10, v11):
    if e1 < e2:
        maxv = np.array([v00, v10])
    else:
        maxv = np.array([v01, v11])
    strike = np.arctan2(maxv[1], maxv[0])
    theta = 90 - m.degrees(strike)
    if theta < 0:
        theta = 180 + theta
    elif theta > 180:
        theta = theta - 180
    return theta


# only used for tape
def max_shortening_azimuth_1d(e1, e2, v00, v01, v10, v11):
    az = np.zeros(len(e1))
    for i in range(len(e1)):
        az[i] = azimuth_math(e1[i], e2[i], v00[i], v01[i], v10[i], v11[i]);
        if az[i] > 180:
            print("OOPSSSSS %d " % i)
    print("Minimum azimuth: %.3f degrees" % np.min(az))
    print("Maximum azimuth: %.3f degrees" % np.max(az))
    return az


def compute_eigenvectors(exx, exy, eyy):
    # exx, eyy can be 1d arrays or 2D arrays
    e1 = np.zeros(np.shape(exx));
    e2 = np.zeros(np.shape(exx));  # eigenvalues
    v00 = np.zeros(np.shape(exx));
    v01 = np.zeros(np.shape(exx));
    v10 = np.zeros(np.shape(exx));
    v11 = np.zeros(np.shape(exx));  # eigenvectors
    dshape = np.shape(exx);
    if len(dshape) == 1:
        for i in range(len(exx)):
            [e11, e22, v] = eigenvector_eigenvalue(exx[i], exy[i], eyy[i]);
            e1[i] = -e11;  # the convention of this code returns negative eigenvalues compared to my other codes.
            e2[i] = -e22;
            v00[i] = v[0][0];
            v10[i] = v[1][0];
            v01[i] = v[0][1];
            v11[i] = v[1][1];
    elif len(dshape) == 2:
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                [e11, e22, v] = eigenvector_eigenvalue(exx[j][i], exy[j][i], eyy[j][i]);
                e1[j][i] = -e11;
                e2[j][i] = -e22;
                v00[j][i] = v[0][0];
                v01[j][i] = v[0][1];
                v10[j][i] = v[1][0];
                v11[j][i] = v[1][1];
    return [e1, e2, v00, v01, v10, v11];


def compute_derived_quantities(exx, exy, eyy):
    # Given the basic components of the strain tensor, compute the rest of the derived quantities
    # like 2nd invariant, azimuth of maximum strain, dilatation, etc.
    # exx, eyy can be 1d arrays or 2D arrays

    I2nd = np.zeros(np.shape(exx));
    max_shear = np.zeros(np.shape(exx));
    dilatation = np.zeros(np.shape(exx));
    azimuth = np.zeros(np.shape(exx));
    [e1, e2, v00, v01, v10, v11] = compute_eigenvectors(exx, exy, eyy);

    dshape = np.shape(exx);
    if len(dshape) == 1:
        datalength = dshape[0];
        print("Computing strain invariants for 1d dataset with length %d." % datalength);
        for i in range(datalength):
            dilatation[i] = e1[i] + e2[i];
            I2nd[i] = np.log10(np.abs(second_invariant(e1[i], 0, e2[i])));
            max_shear[i] = abs((-e1[i] + e2[i]) / 2);
            azimuth[i] = azimuth_math(e1[i], e2[i], v00[i], v01[i], v10[i], v11[i]);
    elif len(dshape) == 2:
        print("Computing strain invariants for 2d dataset.");
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                I2nd[j][i] = np.log10(np.abs(second_invariant(e1[j][i], 0, e2[j][i])));
                max_shear[j][i] = abs((e1[j][i] - e2[j][i]) / 2);
                dilatation[j][i] = e1[j][i] + e2[j][i];
                azimuth[j][i] = azimuth_math(e1[j][i], e2[j][i], v00[j][i], v01[j][i], v10[j][i], v11[j][i]);
    return [I2nd, max_shear, dilatation, azimuth];
