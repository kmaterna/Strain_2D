# Local Average Gradient method, courtesy of Mong-han Huang
# Strain calculation tool based on a certain number of nearby stations

import numpy as np
from strain.models.strain_2d import Strain_2d
from .. import utilities

class loc_avg_grad(Strain_2d):
    """ loc_avg_grad class for 2d strain rate, with general strain_2d behavior """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.xdata, params.ydata,
                           params.outdir);
        self._Name = 'loc_avg_grad'
        self._radiuskm, self._nstations = verify_inputs_loc_avg_grad(params.method_specific);

    def compute(self, myVelfield):
        [Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd] = compute_loc_avg_grad(myVelfield, self._xdata, self._ydata,
                                                                            self._radiuskm, self._nstations);
        # Report observed and residual velocities within bounding box
        filtered_velfield = utilities.filter_by_bounding_box(myVelfield, self._strain_range);
        model_velfield = utilities.create_model_velfield(self._xdata, self._ydata, Ve, Vn, filtered_velfield);
        residual_velfield = utilities.subtract_two_velfields(filtered_velfield, model_velfield);
        return [Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd, filtered_velfield, residual_velfield];


def verify_inputs_loc_avg_grad(method_specific_dict):
    # Takes a dictionary and verifies that it contains the right parameters for loc_avg_grad method
    if 'estimateradiuskm' not in method_specific_dict.keys():
        raise ValueError("\nloc_avg_grad requires estimateradiuskm. Please add to method_specific config. Exiting.\n");
    if 'nstations' not in method_specific_dict.keys():
        raise ValueError("\nloc_avg_grad requires nstations. Please add to method_specific config. Exiting.\n");
    radiuskm = float(method_specific_dict["estimateradiuskm"]);
    nstations = int(method_specific_dict["nstations"]);
    return radiuskm, nstations;


def compute_loc_avg_grad(myVelfield, xlons, ylats, radiuskm, nstations):
    print("------------------------------\nComputing strain via loc_avg_grad method.");

    # Set up grids for the computation
    gx = len(xlons);  # number of x - grid
    gy = len(ylats);  # number of y - grid

    [elon, nlat, e, n, _, _] = velfield_to_LAG_non_utm(myVelfield);
    reflon = np.min([item.elon for item in myVelfield]);
    reflat = np.min([item.nlat for item in myVelfield]);

    # set up a local coordinate reference
    refx = np.min(elon);
    refy = np.min(nlat);
    elon = elon - refx;
    nlat = nlat - refy;

    # Setting calculation parameters
    EstimateRadius = radiuskm * 1000;  # convert to meters
    ns = nstations;  # number of selected stations

    # 2. The main loop, getting displacement gradients around ns nearest stations
    Uxx = np.zeros((gy, gx));
    Uyy = np.zeros((gy, gx));
    Uxy = np.zeros((gy, gx));
    Uyx = np.zeros((gy, gx));
    exx = np.zeros((gy, gx));
    exy = np.zeros((gy, gx));
    eyy = np.zeros((gy, gx));
    rot = np.zeros((gy, gx));
    Ve = np.zeros((gy, gx));
    Vn = np.zeros((gy, gx));
    for i in range(gx):
        for j in range(gy):
            [gridX_loc, gridY_loc] = convert_to_local_planar(xlons[i], ylats[j], reflon, reflat);
            X = elon;   # in local coordinates, m
            Y = nlat;   # in local coordiantes, m
            l1 = len(elon);

            # the distance from stations to grid
            r = np.zeros((l1,));
            for ii in range(l1):
                r[ii] = np.sqrt((gridX_loc - X[ii]) * (gridX_loc - X[ii]) + (gridY_loc - Y[ii]) * (gridY_loc - Y[ii]));
                # for a particular grid point, we are computing the distance to every station

            stations = np.zeros((l1, 5));
            stations[:, 0] = r;
            stations[:, 1] = elon.reshape((l1,));
            stations[:, 2] = nlat.reshape((l1,));
            stations[:, 3] = e.reshape((l1,));
            stations[:, 4] = n.reshape((l1,));

            iX = stations[stations[:, 0].argsort()];  # sort data, iX represented the order of the sorting of 1st column
            # we choose the first ns data
            SelectStations = np.zeros((ns, 5));
            for iiii in range(ns):
                SelectStations[iiii, :] = iX[iiii, :];
            # print(SelectStations);
            Px, Py = 0, 0;
            Px2, Py2, Pxy = 0, 0, 0;
            dU, dV = 0, 0;
            dxU, dyU = 0, 0;
            dxV, dyV = 0, 0;  # should be inside the if statement or not?
            # print(SelectStations[ns-1, 0]);  # the distance of the cutoff station (in m)
            if SelectStations[ns-1, 0] <= EstimateRadius:
                for iii in range(ns):
                    Px = Px + SelectStations[iii, 1];
                    Py = Py + SelectStations[iii, 2];
                    Px2 = Px2 + SelectStations[iii, 1] * SelectStations[iii, 1];
                    Py2 = Py2 + SelectStations[iii, 2] * SelectStations[iii, 2];
                    Pxy = Pxy + SelectStations[iii, 1] * SelectStations[iii, 2];
                    dU = dU + SelectStations[iii, 3];
                    dxU = dxU + SelectStations[iii, 1] * SelectStations[iii, 3];
                    dyU = dyU + SelectStations[iii, 2] * SelectStations[iii, 3];
                    dV = dV + SelectStations[iii, 4];
                    dxV = dxV + SelectStations[iii, 1] * SelectStations[iii, 4];
                    dyV = dyV + SelectStations[iii, 2] * SelectStations[iii, 4];
            G = [[ns, Px, Py], [Px, Px2, Pxy], [Py, Pxy, Py2]];
            if np.sum(G) == ns:
                Uxx[j, i], Uyy[j, i], Uxy[j, i], Uyx[j, i] = 0, 0, 0, 0;
            else:
                dx = np.array([[dU], [dxU], [dyU]]);
                dy = np.array([[dV], [dxV], [dyV]]);
                modelx = np.dot(np.linalg.inv(G), dx);
                modely = np.dot(np.linalg.inv(G), dy);
                Uxx[j, i] = modelx[1];
                Uyy[j, i] = modely[2];
                Uxy[j, i] = modelx[2];
                Uyx[j, i] = modely[1];

                Ve[j, i] = 1000 * (modelx[0] + modelx[1]*SelectStations[0, 1] + modelx[2]*SelectStations[0, 2]);  # mm
                Vn[j, i] = 1000 * (modely[0] + modely[1]*SelectStations[0, 1] + modely[2]*SelectStations[0, 2]);  # mm

            # skipping misfit right now
            # misfit estimation   d = m1 + m2 x + m3 y

            # 3. Moving on to strain calculation
            sxx = Uxx[j, i];
            syy = Uyy[j, i];
            sxy = .5 * (Uxy[j, i] + Uyx[j, i]);
            omega = .5 * (Uxy[j, i] - Uyx[j, i]);
            exx[j, i] = sxx * 1e9;
            exy[j, i] = sxy * 1e9;
            eyy[j, i] = syy * 1e9;
            rot[j, i] = omega * 1e9;

    print("Success computing strain via loc_avg_grad method.\n");

    return [Ve, Vn, rot, exx, exy, eyy];


def velfield_to_LAG_non_utm(myVelfield):
    elon, nlat = [], [];
    e, n, esig, nsig = [], [], [], [];
    elon_all = [item.elon for item in myVelfield];
    nlat_all = [item.nlat for item in myVelfield];
    reflon = np.min(elon_all);
    reflat = np.min(nlat_all);
    for item in myVelfield:
        [x_meters, y_meters] = convert_to_local_planar(item.elon, item.nlat, reflon, reflat)
        elon.append(x_meters);
        nlat.append(y_meters);
        e.append(item.e*0.001);
        n.append(item.n*0.001);
        esig.append(item.se*0.001);
        nsig.append(item.sn*0.001);
    return [np.array(elon), np.array(nlat), np.array(e), np.array(n), np.array(esig), np.array(nsig)];


def convert_to_local_planar(lon, lat, reflon, reflat):
    earth_r = 6371000;  # earth's radius
    x_deg = np.subtract(lon, reflon);
    y_deg = np.subtract(lat, reflat);
    x_meters = x_deg * earth_r * np.cos(np.deg2rad(lat)) * 2 * np.pi / 360;
    y_meters = y_deg * earth_r * 2 * np.pi / 360;
    return [x_meters, y_meters];
