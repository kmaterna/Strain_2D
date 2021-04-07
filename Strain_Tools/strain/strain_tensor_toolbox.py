# A math toolbox for strain tensor components.
# Rotation, 2nd invariant, Dilatation, Shear Strain, Eigenvectors, Eigenvalues, etc.


import numpy as np
import math as m


def second_invariant(exx, exy, eyy):
    """
    :param exx: strain component
    :type exx: float
    :param exy: strain component
    :type exy: float
    :param eyy: strain component
    :type eyy: float
    :returns: second invariant, in units of (strain_component)^2
    :rtype: float
    """
    e2nd = exx * eyy - exy * exy;
    return e2nd;


def eigenvector_eigenvalue(exx, exy, eyy):
    """
    :param exx: strain component
    :type exx: float
    :param exy: strain component
    :type exy: float
    :param eyy: strain component
    :type eyy: float
    :returns: [eigenvalue1 eigenvalue2 eigenvectors]
    :rtype: list
    """
    if np.isnan(np.sum([exx, exy, eyy])):
        v = [[np.nan, np.nan], [np.nan, np.nan]];
        return [0, 0, v];
    T = np.array([[exx, exy], [exy, eyy]]);  # the tensor
    w, v = np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
    return [w[0], w[1], v];


def max_shear_strain(exx, exy, eyy):
    """
    :param exx: strain component
    :type exx: float
    :param exy: strain component
    :type exy: float
    :param eyy: strain component
    :type eyy: float
    :returns: max shear, in units of strain_component
    :rtype: float
    """
    if np.isnan(np.sum([exx, exy, eyy])):
        return 0;
    T = np.array([[exx, exy], [exy, eyy]]);  # the tensor
    w, v = np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
    # w = eigenvalues; v = eigenvectors
    max_shear = (w[0] - w[1]) * 0.5;
    return max_shear;


def compute_displacement_gradients(up, vp, ur, vr, uq, vq, dx, dy):
    """
    up, vp : velocity at a reference point P
    uq, vq, dx : velocity at point Q, which is offset from reference point P by dx in the x direction
    ur, vr, dy : velocity at point R, which is offset from reference point P by dy in the y direction
    In practical usage, these are in mm/yr and km.
    """
    dudx = (uq - up) / dx;
    dvdx = (vq - vp) / dx;
    dudy = (ur - up) / dy;
    dvdy = (vr - vp) / dy;
    return [dudx, dvdx, dudy, dvdy];


def compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy):
    """
    Given a displacement tensor, compute the relevant parts of the strain and rotation tensors.
    Also converts to nanostrain per year.
    Rot is the off-diagonal element of the rotation tensor
    Rot has native units radians/year. Here we return radians per 1000 yrs (easier to interpret numbers)
    http://www.engr.colostate.edu/~thompson/hPage/CourseMat/Tutorials/Solid_Mechanics/rotations.pdf

    :param dudx: displacement gradient
    :type dudx: float
    :param dvdx: displacement gradient
    :type dvdx: float
    :param dudy: displacement gradient
    :type dudy: float
    :param dvdy: displacement gradient
    :type dvdy: float
    :returns: strain and rotation components
    :rtype: list
    """
    exx = dudx * 1000;
    exy = (0.5 * (dvdx + dudy)) * 1000;
    eyy = dvdy * 1000;
    rot = (0.5 * (dvdx - dudy));
    rot = rot * 1000.0;
    return [exx, exy, eyy, rot];


def compute_max_shortening_azimuth(e1, e2, v00, v01, v10, v11):
    """
    :param e1: eigenvalue 1
    :type e1: float
    :param e2: eigenvalue 2
    :type e2: float
    :param v00: eigenvector 1
    :type v00: float
    :param v01: eigenvector 1
    :type v01: float
    :param v10: eigenvector 2
    :type v10: float
    :param v11: eigenvector 2
    :type v11: float
    :returns: azimuth of maximum shortening axis, in degrees CW from north
    :rtype: float
    """
    if e1 < e2:
        maxv = np.array([v00, v10])
    else:
        maxv = np.array([v01, v11])
    strike = np.arctan2(maxv[1], maxv[0])
    theta = 90 - m.degrees(strike)
    if np.isnan(theta):
        return np.nan;
    if theta < 0:
        theta = 180 + theta
    elif theta > 180:
        theta = theta - 180
    assert (theta < 180), ValueError("Error: computing an azimuth over 180 degrees.");
    return theta


def compute_eigenvectors(exx, exy, eyy):
    """
    exx, eyy can be 1d arrays or 2D arrays

    :param exx: strain component, float or 1d array
    :param exy: strain component, float or 1d array
    :param eyy: strain component, float or 1d array
    :rtype: list
    """
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
            e1[i] = e11;  # the convention of this code returns negative eigenvalues compared to my other codes.
            e2[i] = e22;
            v00[i] = v[0][0];
            v10[i] = v[1][0];
            v01[i] = v[0][1];
            v11[i] = v[1][1];
    elif len(dshape) == 2:
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                [e11, e22, v] = eigenvector_eigenvalue(exx[j][i], exy[j][i], eyy[j][i]);
                e1[j][i] = e11;
                e2[j][i] = e22;
                v00[j][i] = v[0][0];
                v01[j][i] = v[0][1];
                v10[j][i] = v[1][0];
                v11[j][i] = v[1][1];
    return [e1, e2, v00, v01, v10, v11];


def compute_derived_quantities(exx, exy, eyy):
    """
    Given the basic components of the strain tensor, compute the rest of the derived quantities
    like 2nd invariant, azimuth of maximum strain, dilatation, etc.
    exx, eyy can be 1d arrays or 2D arrays

    :param exx: strain component, float or 1d array
    :param exy: strain component, float or 1d array
    :param eyy: strain component, float or 1d array
    :rtype: list
    """

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
            azimuth[i] = compute_max_shortening_azimuth(e1[i], e2[i], v00[i], v01[i], v10[i], v11[i]);
    elif len(dshape) == 2:
        print("Computing strain invariants for 2d dataset.");
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                I2nd[j][i] = np.log10(np.abs(second_invariant(e1[j][i], 0, e2[j][i])));
                max_shear[j][i] = abs((e1[j][i] - e2[j][i]) / 2);
                dilatation[j][i] = e1[j][i] + e2[j][i];
                azimuth[j][i] = compute_max_shortening_azimuth(e1[j][i], e2[j][i], v00[j][i], v01[j][i],
                                                               v10[j][i], v11[j][i]);
    return [I2nd, max_shear, dilatation, azimuth];


def angle_mean_math(azimuth_values):
    """
    :param azimuth_values: azimuths in degrees
    :type azimuth_values: list
    :returns: an average azimuth and standard deviation of azimuths, in degrees
    :rtype: float
    """
    sin_list, cos_list = [], [];
    for phi in azimuth_values:
        sin_list.append(np.sin(2 * np.radians(90 - phi)));
        cos_list.append(np.cos(2 * np.radians(90 - phi)));
    s = np.nanmean(sin_list);
    c = np.nanmean(cos_list);
    R = ((s ** 2 + c ** 2) ** .5)
    sd = np.degrees((-2 * np.log(R)) ** .5) / 2
    # V = 1 - R
    # sd = np.degrees((2*V)**.5)
    # t = np.arctan2(s, c)
    # strike = R*math.e**(math.i*t)
    strike = np.arctan2(s, c) / 2
    theta = 90 - np.degrees(strike)
    if theta < 0:
        theta = 180 + theta
    elif theta > 180:
        theta = theta - 180
    return theta, sd;
