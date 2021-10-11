# A math toolbox for strain tensor components.
# Rotation, 2nd invariant, Dilatation, Shear Strain, Eigenvectors, Eigenvalues, etc.


import numpy as np
import math as m


def strain_on_regular_grid(dx, dy, V1, V2):
    """Compute strain rate on a regular grid.
    :param dx: 1d array
    :param dy: 1d array
    :param V1: 2d array of velocities in x-direction
    :param V2: 2d array of velocities in y-direction
    :returns: 2d arrays of derived quantities
    """
    e11, dV1dx2 = np.gradient(V1, dx, dy)
    dV2dx1, e22 = np.gradient(V2, dx, dy)

    # shear strain rate is the symmetric part of the gradient
    e12 = 0.5*(dV1dx2 + dV2dx1)

    # Rotation is the anti-symmetric component of the displacement gradient tensor
    rot = 0.5*(dV1dx2 - dV2dx1)

    # all the strain rates in each grid location
    return e11, e22, e12, rot


def compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy):
    """
    Given a displacement tensor, compute the relevant parts of the strain and rotation tensors.
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
    exx = dudx;
    exy = (0.5 * (dvdx + dudy));
    eyy = dvdy;
    rot = (0.5 * (dvdx - dudy));
    return [exx, exy, eyy, rot];


def compute_derived_quantities(exx, exy, eyy):
    """
    Given the basic components of the strain tensor, compute the rest of the derived quantities
    like 2nd invariant, azimuth of maximum strain, dilatation, etc.
    exx, eyy can be 1d arrays or 2D arrays

    :param exx: strain component, float or 1d array or 2d array
    :param exy: strain component, float or 1d array or 2d array
    :param eyy: strain component, float or 1d array or 2d array
    :rtype: list
    """
    # Since exx etc. are numpy arrays, we can use numpy's vectorized math
    I2nd = np.log10(np.abs(exx*eyy - np.square(exy)))
    max_shear = np.sqrt(np.square(exx - eyy) + 4*np.square(exy))
    dilatation = exx + eyy

    # Azimuth is tricky so leaving it as a for-loop for now
    azimuth = np.zeros(np.shape(exx));
    [e1, e2, v00, v01, v10, v11] = compute_eigenvectors(exx, exy, eyy);

    dshape = np.shape(exx);
    if len(dshape) == 1:
        datalength = dshape[0];
        print("Computing strain invariants for 1d dataset with length %d." % datalength);
        for i in range(datalength):
            azimuth[i] = compute_max_shortening_azimuth(e1[i], e2[i], v00[i], v01[i], v10[i], v11[i]);
    elif len(dshape) == 2:
        print("Computing strain invariants for 2d dataset.");
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                azimuth[j][i] = compute_max_shortening_azimuth(e1[j][i], e2[j][i], v00[j][i], v01[j][i],
                                                               v10[j][i], v11[j][i]);
    return [I2nd, max_shear, dilatation, azimuth];


def compute_eigenvectors(exx, exy, eyy):
    """
    exx, eyy can be 1d arrays or 2D arrays

    :param exx: strain component, float or 1d array
    :param exy: strain component, float or 1d array
    :param eyy: strain component, float or 1d array
    :rtype: list
    """
    e1, e2 = np.zeros(np.shape(exx)), np.zeros(np.shape(exx));  # eigenvalues
    v00, v01 = np.zeros(np.shape(exx)), np.zeros(np.shape(exx));
    v10, v11 = np.zeros(np.shape(exx)), np.zeros(np.shape(exx));  # eigenvectors
    dshape = np.shape(exx);
    if len(dshape) == 1:
        for i in range(len(exx)):
            [e11, e22, v] = eigenvector_eigenvalue(exx[i], exy[i], eyy[i]);
            e1[i], e2 = e11, e22;  # convention of this code returns negative eigenvalues compared to my other codes
            v00[i], v10[i] = v[0][0], v[1][0];
            v01[i], v11[i] = v[0][1], v[1][1];
    elif len(dshape) == 2:
        for j in range(dshape[0]):
            for i in range(dshape[1]):
                [e11, e22, v] = eigenvector_eigenvalue(exx[j][i], exy[j][i], eyy[j][i]);
                e1[j][i], e2[j][i] = e11, e22;
                v00[j][i], v01[j][i] = v[0][0], v[0][1];
                v10[j][i], v11[j][i] = v[1][0], v[1][1];
    return [e1, e2, v00, v01, v10, v11];


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
