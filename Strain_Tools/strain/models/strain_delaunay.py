# # delaunay_strain on a sphere (c/o of Bill Hammond)
# # Based on Savage (2001)
# # 
# % [e_phiphi,e_thetaphi,e_thetatheta,omega_r,U_theta,U_phi,...
# %          s_omega_r,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,s_U_theta,...
# %          s_U_phi,chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p,r_PHITHETA,u_phi_p,u_theta_p]=...
# %          strain_sphere(phi,theta,u_phi,u_theta,s_phi,s_theta,weight,paramsel)
# %
# % INPUTS:
# %   phi is longitude, entered in DEGREES
# %   theta is COlatitude, entered in DEGREES (note: same as latitude minus 90 degrees).
# %   u_phi is velocity in phi direction, entered in meters per yr.
# %   u_theta is velocity in theta direction, entered in meters per yr.
# %   s_phi is uncertainty in velocity in phi direction, entered in meters per yr.
# %   s_theta is uncertainty in velocity in theta direction, entered in meters per yr.
# %   weight = 1 if you want to use uncertainties as weights in the inversion.
# %   paramsel = 1 if you want to invert for solid body rotations only (not strains).
# %   paramsel = 2 if you want to invert for strains only.
# %   paramsel = anything else, if you want to invert for both
# %
# % OUTPUTS:
# %   e_phiphi = extensional strain rate along phi direction
# %   e_thetaphi = shear strain rate 
# %   e_thetatheta = extensional strain rate along theta direction
# %   omega_r = average rotation rate around vertical axis
# %   U_phi = average translation rate phi-ward (meters!)            
# %   U_theta = average translation rate theta-ward (meters!)          
# %   s_e_phiphi = uncertainty in extensional strain rate along phi direction
# %   s_e_thetaphi = uncertainty in shear strain rate 
# %   s_e_thetatheta = uncertainty in extensional strain rate along theta direction
# %   s_omega_r = uncertainty in average rotation rate along vertical axis
# %   chi2 is the misfit between the strain rate estimate and the data
# %   OMEGA  = total rotation rate associated with Euler pole
# %   THETA_p = colatitude of Euler Pole
# %   PHI_p = longitude of Euler Pole
# %   s_OMEGA  = unc. in total rotation rate associated with Euler pole
# %   s_THETA_p = unc. in colatitude of Euler Pole
# %   s_PHI_p = unc. in longitude of Euler Pole
# %   r_PHITHETA = is correllation coefficient between PHI & THETA for
# %               plotting uncertainty ellipses in gmt.
# %   u_phi_p, u_theta_p are the predicted velocities at each location
# %
# % note that theta is COlatitude AND that thetaward is Southing
# %                    ^^
# %
# % see Savage et al., JGR October 2001, p.22,005 for derivation
# %
# % Version Oct. 12, 2004


import numpy as np
from scipy.spatial import Delaunay
from .. import output_manager, produce_gridded
from strain.models.strain_2d import Strain_2d


class delaunay(Strain_2d):
    """ Delaunay class for 2d strain rate """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir)
        self._Name = 'delaunay'

    def compute(self, myVelfield):
        print("------------------------------\nComputing strain via Delaunay on a sphere, and converting to a grid.");

        [xcentroid, ycentroid, triangle_verts, rot, exx, exy, eyy] = compute_with_delaunay_polygons(myVelfield);

        lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd = produce_gridded.tri2grid(self._grid_inc, self._strain_range,
                                                                                  triangle_verts, rot, exx, exy, eyy);

        # Here we output convenient things on polygons, since it's intuitive for the user.
        # output_manager.outputs_1d(xcentroid, ycentroid, triangle_verts, rot, exx, exy, eyy, self._strain_range,
        #                           myVelfield, self._outdir);

        # Velocities aren't used in Delaunay
        Ve, Vn = np.nan*np.empty(exx_grd.shape), np.nan*np.empty(exx_grd.shape)

        print("Success computing strain via Delaunay method.\n");
        return [lons, lats, Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd];


def compute_with_delaunay_polygons(myVelfield):
    elon = [x.elon for x in myVelfield];
    nlat = [x.nlat for x in myVelfield];
    e = [x.e for x in myVelfield];
    n = [x.n for x in myVelfield];
    se = [x.se for x in myVelfield];
    sn = [x.sn for x in myVelfield];
    z = np.array([elon, nlat]);
    z = z.T;
    tri = Delaunay(z);

    triangle_vertices = z[tri.simplices];
    trishape = np.shape(triangle_vertices);  # 516 x 3 x 2, for example
    print("Number of triangle elements: %d" % (trishape[0]));

    # We are going to solve for the velocity gradient tensor at the centroid of each triangle.
    centroids = [];
    for i in range(trishape[0]):
        xcor_mean = np.mean(
                [
                    triangle_vertices[i, 0, 0], 
                    triangle_vertices[i, 1, 0], 
                    triangle_vertices[i, 2, 0]
                ]
            );
        ycor_mean = np.mean(
                [
                    triangle_vertices[i, 0, 1], 
                    triangle_vertices[i, 1, 1], 
                    triangle_vertices[i, 2, 1]
                ]
            );
        centroids.append([xcor_mean, ycor_mean]);
    xcentroid = [x[0] for x in centroids];
    ycentroid = [x[1] for x in centroids];

    # Initialize arrays.
    rot = [];
    exx, exy, eyy = [], [], [];

    # for each triangle:
    for i in range(trishape[0]):
        # Get the velocities of each vertex (VE1, VN1, VE2, VN2, VE3, VN3)
        # Get velocities for Vertex 1 (triangle_vertices[i,0,0] and triangle_vertices[i,0,1])
        xindex1 = np.where(elon == triangle_vertices[i, 0, 0])
        yindex1 = np.where(nlat == triangle_vertices[i, 0, 1])
        index1 = int(np.intersect1d(xindex1, yindex1)[0]);
        xindex2 = np.where(elon == triangle_vertices[i, 1, 0])
        yindex2 = np.where(nlat == triangle_vertices[i, 1, 1])
        index2 = int(np.intersect1d(xindex2, yindex2)[0]);
        xindex3 = np.where(elon == triangle_vertices[i, 2, 0])
        yindex3 = np.where(nlat == triangle_vertices[i, 2, 1])
        index3 = int(np.intersect1d(xindex3, yindex3)[0]);

        phi = np.array(
                [
                    triangle_vertices[i, 0, 0], 
                    triangle_vertices[i, 1, 0], 
                    triangle_vertices[i, 2, 0]
                ]
            );
        theta = np.array(
                [
                    triangle_vertices[i, 0, 1], 
                    triangle_vertices[i, 1, 1], 
                    triangle_vertices[i, 2, 1]
                ]
            );
        theta = [i - 90 for i in theta];
        u_phi = np.array([e[index1], e[index2], e[index3]]);
        u_theta = np.array([n[index1], n[index2], n[index3]]);
        u_theta = np.array([-i for i in u_theta]);  # colatitude needs negative theta values.
        s_phi = np.array([se[index1], se[index2], se[index3]]);
        s_theta = np.array([sn[index1], sn[index2], sn[index3]]);


        # HERE WE PLUG IN BILL'S CODE!
        weight = 1;
        paramsel = 0;
        [e_phiphi, e_thetaphi, e_thetatheta, omega_r, U_theta, U_phi, s_omega_r, \
            s_e_phiphi, s_e_thetaphi, s_e_thetatheta, s_U_theta, s_U_phi, chi2, \
            OMEGA, THETA_p, PHI_p, s_OMEGA, s_THETA_p, s_PHI_p, r_PHITHETA, u_phi_p, \
            u_theta_p] = strain_sphere(
                phi, 
                theta, 
                u_phi, 
                u_theta, 
                s_phi, 
                s_theta, 
                weight, 
                paramsel
            );

        # The components that are easily computed
        # Units: nanostrain per year.
        # There might be a sign issue here compared to other codes.
        exx.append(-e_phiphi * 1e6);
        exy.append(e_thetaphi * 1e6);
        eyy.append(-e_thetatheta * 1e6);

        # # Compute a number of values based on tensor properties.
        rot.append(OMEGA * 1000 * 1000);

    return [xcentroid, ycentroid, triangle_vertices, rot, exx, exy, eyy];


def strain_sphere(phi, theta, u_phi, u_theta, s_phi, s_theta, weight, paramsel):
    theta = [i * np.pi / 180 for i in theta];  # convert to radians
    phi = [i * np.pi / 180 for i in phi];
    r0 = 6.378e6;  # mean equitorial Earth radius

    theta_0 = np.mean(theta);
    phi_0 = np.mean(phi);

    m = len(np.shape(u_phi));  # 1 for a column vector
    m2 = len(np.shape(u_theta));  # n was the first element here.
    n = np.shape(u_theta)[0];
    if m != 1 or m2 != 1:
        print('u_phi and u_theta must be column vectors');

        # check for colinearity
    colin = 0;
    dc = [90 - i for i in theta];
    Gc = np.column_stack((np.ones([len(phi), 1]), phi));  # make a 2d array
    mc = np.linalg.lstsq(Gc, dc, rcond=None)[0];
    resc = np.dot(Gc, mc) - dc;  # residuals
    if max(abs(resc)) < 1e-5:
        colin = 1;
    if all(abs(phi - np.mean(phi)) < 1e-5) or all(abs(theta - np.mean(theta)) < 1e-5):
        colin = 2;
    if colin == 1 or colin == 2:
        print(['Warning:  Points are collinear. Type=' + str(colin)]);

    # #%% Make the matrix needed for least squares inversion
    temparray = np.array([u_phi.T, u_theta.T])
    newd = temparray.reshape((6, 1), order='F')
    d = newd.T;
    d = d.flatten();  # this performs the matlab reshape() function.

    if paramsel != 2 and paramsel != 1:  # if we're solving for both strain and rotation.
        G = np.zeros([2 * n, 6]);
    else:
        G = np.zeros([2 * n, 3]);  # if paramsel == 2
    if paramsel != 2:
        for i in range(n):
            del_phi = phi[i] - phi_0;
            del_theta = theta[i] - theta_0;

            G[2 * i, 0] = -r0;
            G[2 * i, 1] = -r0 * np.cos(theta_0) * del_phi;
            G[2 * i, 2] = r0 * del_theta;
            if paramsel != 1:
                G[2 * i, 3] = r0 * np.sin(theta_0) * del_phi;
                G[2 * i, 4] = r0 * del_theta;
                G[2 * i, 5] = 0;

            G[2 * i + 1, 0] = -r0 * np.cos(theta_0) * del_phi;
            G[2 * i + 1, 1] = r0;
            G[2 * i + 1, 2] = -r0 * np.sin(theta_0) * del_phi;
            if paramsel != 1:
                G[2 * i + 1, 3] = 0;
                G[2 * i + 1, 4] = r0 * np.sin(theta_0) * del_phi;
                G[2 * i + 1, 5] = r0 * del_theta;
    else:
        for i in range(n):
            del_phi = phi[i] - phi_0;
            del_theta = theta[i] - theta_0;

            G[2 * i, 0] = r0 * np.sin(theta_0) * del_phi;
            G[2 * i, 1] = r0 * del_theta;
            G[2 * i, 2] = 0;

            G[2 * i + 1, 0] = 0;
            G[2 * i + 1, 1] = r0 * np.sin(theta_0) * del_phi;
            G[2 * i + 1, 2] = r0 * del_theta;

            # #% perform the inversion as in an overdetermined system

    covd = np.diag(np.hstack((np.square(s_phi), np.square(s_theta))));
    # Same as Matlab's

    if weight == 1:
        W = np.diag(1. / np.diag(covd));  # this checks out.
        M = np.dot(np.linalg.inv(np.dot(np.dot(G.T, W), G)), np.dot(G.T, W));
    else:
        M = np.dot(np.linalg.inv(np.dot(G.T, G)), G.T);

    m = np.dot(M, d);
    # m is the same between python and matlab

    dpred = np.dot(G, m);

    # #the predicted u_phi_p, u_theta_p;
    u_phi_p = dpred[::2];  # elements 0, 2, 4...
    u_theta_p = dpred[1::2];  # elements 1, 3, 5...

    if paramsel == 2:
        omega_theta = np.nan;
        omega_phi = np.nan;
        omega_r = np.nan;
        e_phiphi = m[0];
        e_thetaphi = m[1];
        e_thetatheta = m[2];
    elif paramsel == 1:
        omega_theta = m[0];
        omega_phi = m[1];
        omega_r = m[2];
        e_phiphi = np.nan;
        e_thetaphi = np.nan;
        e_thetatheta = np.nan;
    else:
        omega_theta = m[0];
        omega_phi = m[1];
        omega_r = m[2];
        e_phiphi = m[3];
        e_thetaphi = m[4];
        e_thetatheta = m[5];

    U_theta = r0 * omega_phi;
    U_phi = -r0 * omega_theta;

    covm = np.dot(np.dot(M, covd), M.T);

    # # % obtain the uncertainties
    if paramsel == 2:
        s_e_phiphi = np.sqrt(covm[0, 0]);
        s_e_thetaphi = np.sqrt(covm[1, 1]);
        s_e_thetatheta = np.sqrt(covm[2, 2]);
        s_U_phi = 0;
        s_U_theta = 0;
        s_omega_r = np.nan;
        s_omega_phi = np.nan;
        s_omega_theta = np.nan;
    elif paramsel == 1:
        s_e_phiphi = 0;
        s_e_thetaphi = 0;
        s_e_thetatheta = 0;
        s_U_phi = r0 * np.sqrt(covm[0, 0]);
        s_U_theta = r0 * np.sqrt(covm[1, 1]);
        s_omega_r = np.sqrt(covm[2, 2]);
        s_omega_phi = np.sqrt(covm[1, 1]);
        s_omega_theta = np.sqrt(covm[0, 0]);
    else:
        s_U_phi = r0 * np.sqrt(covm[0, 0]);
        s_U_theta = r0 * np.sqrt(covm[1, 1]);
        s_omega_r = np.sqrt(covm[2, 2]);
        s_omega_phi = np.sqrt(covm[1, 1]);
        s_omega_theta = np.sqrt(covm[0, 0]);
        s_e_phiphi = np.sqrt(covm[3, 3]);
        s_e_thetaphi = np.sqrt(covm[4, 4]);
        s_e_thetatheta = np.sqrt(covm[5, 5]);

    N = len(d);
    first_term = np.linalg.lstsq(covd.T, (d - np.dot(G, m)).T, rcond=None)[0].T;
    if paramsel == 1 or paramsel == 2:
        chi2 = (np.dot(first_term, (d - np.dot(G, m)))) / (N);  # N-3 but a divide by zero happened
    else:
        chi2 = (np.dot(first_term, (d - np.dot(G, m)))) / (N);  # N-6 but a divide by zero happened

    # # % Compute the Euler Vectors for the solid body rotation
    # # % from (A7) of Savage et al., October 2001, JGR appendix
    if paramsel != 2:
        OMEGA = np.sqrt(np.square(omega_r) + np.square(omega_phi) + np.square(omega_theta));
        y = omega_r * np.cos(theta_0) - omega_theta * np.sin(theta_0);
        THETA_p = np.arccos(y / OMEGA);
        A = omega_r * np.sin(theta_0) * np.sin(phi_0) + omega_theta * np.cos(theta_0) * np.sin(
            phi_0) + omega_phi * np.cos(phi_0);
        B = omega_r * np.sin(theta_0) * np.cos(phi_0) + omega_theta * np.cos(theta_0) * np.cos(
            phi_0) - omega_phi * np.sin(phi_0);
        PHI_p = np.arctan2(A, B);
        z = 1 / np.sqrt(1 - (np.square(y / OMEGA)));

        # %cast the problem of finding the uncertainties as a linear
        # %inverse problem (linearized) of
        J = np.zeros([3, 3]);
        J[0, 0] = omega_r / OMEGA;
        J[0, 1] = omega_phi / OMEGA;
        J[0, 2] = omega_theta / OMEGA;

        J[1, 0] = (np.sin(theta_0) / (A * A + B * B)) * (B * np.sin(phi_0) - A * np.cos(phi_0));
        J[1, 1] = (1 / (A * A + B * B)) * (B * np.cos(phi_0) + A * np.sin(phi_0));
        J[1, 2] = (np.cos(theta_0) / (A * A + B * B)) * (B * np.sin(phi_0) - A * np.cos(phi_0));

        J[2, 0] = -(z / OMEGA) * (np.cos(theta_0) - y * omega_r / (OMEGA * OMEGA));
        J[2, 1] = z * y * omega_phi / (OMEGA * OMEGA * OMEGA);
        J[2, 2] = (z / OMEGA) * (np.sin(theta_0) + y * omega_theta / (OMEGA * OMEGA));

        covp = np.array([[covm[2, 2], covm[2, 1], covm[2, 0]], [covm[2, 1], covm[1, 1], covm[0, 1]],
                         [covm[2, 0], covm[0, 1], covm[0, 0]]]);
        covE = np.dot(np.dot(J, covp), J.T);

        s_OMEGA = np.sqrt(covE[0, 0]);
        s_PHI_p = np.sqrt(covE[1, 1]);
        s_THETA_p = np.sqrt(covE[2, 2]);

        r_PHITHETA = covE[1, 2] / np.sqrt(covE[1, 1] * covE[2, 2]);

    else:
        OMEGA = np.nan;
        THETA_p = np.nan
        PHI_p = np.nan;
        r_PHITHETA = np.nan;
        s_OMEGA = np.nan;
        s_PHI_p = np.nan;
        s_THETA_p = np.nan;

    return [
            e_phiphi, 
            e_thetaphi, 
            e_thetatheta, 
            omega_r, 
            U_theta, 
            U_phi, 
            s_omega_r, 
            s_e_phiphi, 
            s_e_thetaphi,
            s_e_thetatheta, 
            s_U_theta, 
            s_U_phi, 
            chi2, 
            OMEGA, 
            THETA_p, 
            PHI_p, 
            s_OMEGA, 
            s_THETA_p, 
            s_PHI_p, 
            r_PHITHETA,
            u_phi_p, 
            u_theta_p
        ];


def print_all_values(
        e_phiphi, 
        e_thetaphi, 
        e_thetatheta,
        omega_r, 
        U_theta, 
        U_phi, 
        s_omega_r, 
        s_e_phiphi, 
        s_e_thetaphi,
        s_e_thetatheta, 
        s_U_theta, 
        s_U_phi, 
        chi2, 
        OMEGA, 
        THETA_p, 
        PHI_p, 
        s_OMEGA, 
        s_THETA_p, 
        s_PHI_p,
        r_PHITHETA, 
        u_phi_p, 
        u_theta_p
    ):
    print("e_phiphi is:     " + str(e_phiphi));
    print("e_thetaphi is:   " + str(e_thetaphi));
    print("e_thetatheta is: " + str(e_thetatheta));
    print("omega_r is:      " + str(omega_r));
    print("U_theta is:      " + str(U_theta));
    print("U_phi is:        " + str(U_phi));
    print("s_omega_r is:    " + str(s_omega_r));
    print("s_e_phiphi is:   " + str(s_e_phiphi));
    print("s_e_thetaphi is: " + str(s_e_thetaphi));
    print("s_e_thetatheta:  " + str(s_e_thetatheta));
    print("s_U_theta is:    " + str(s_U_theta));
    print("s_U_phi is:      " + str(s_U_phi));
    print("chi2 is:         " + str(chi2));
    print("OMEGA is:        " + str(OMEGA));
    print("THETA_p is:      " + str(THETA_p));
    print("PHI_p is:        " + str(PHI_p));
    print("s_OMEGA is:      " + str(s_OMEGA));
    print("s_THETA_p is:    " + str(s_THETA_p));
    print("s_PHI_p is:      " + str(s_PHI_p));
    print("r_PHITHETA is:   " + str(r_PHITHETA));
    print("u_phi_p is:      " + str(u_phi_p));
    print("u_theta_p is:    " + str(u_theta_p));


if __name__ == "__main__":
    phi = np.array([-123, -123, -123.5]);
    theta_lat = np.array([40, 40.25, 40.5]);
    theta = [90 - i for i in theta_lat];
    u_phi = np.array([0.023, 0.025, 0.027]);
    u_theta = np.array([0.013, 0.011, 0.011]);
    s_phi = np.array([0.002, 0.002, 0.002]);
    s_theta = np.array([0.002, 0.002, 0.002]);
    weight = 1;
    paramsel = 0;
    # [e_phiphi,e_thetaphi,e_thetatheta,omega_r,U_theta,U_phi,s_omega_r,s_e_phiphi,s_e_thetaphi,s_e_thetatheta,s_U_theta,s_U_phi,chi2,OMEGA,THETA_p,PHI_p,s_OMEGA,s_THETA_p,s_PHI_p,r_PHITHETA,u_phi_p,u_theta_p] = strain_sphere(phi,theta,u_phi,u_theta,s_phi,s_theta,weight,paramsel);
