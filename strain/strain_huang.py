# Courtesy of Mong-han Huang
# Strain calculation tool based on a certain number of nearby stations

import numpy as np
from Tectonic_Utils.geodesy import utm_conversion


# ----------------- COMPUTE -------------------------
def compute(myVelfield, MyParams):
    print("------------------------------\nComputing strain via Huang method.");

    # Set up grids for the computation
    xlons = np.arange(MyParams.range_strain[0], MyParams.range_strain[1], MyParams.inc[0]);
    ylats = np.arange(MyParams.range_strain[2], MyParams.range_strain[3], MyParams.inc[1]);
    gx = len(xlons);  # number of x - grid
    gy = len(ylats);  # number of y - grid

    [elon, nlat, e, n, esig, nsig] = velfield_to_huang_format(myVelfield);

    # set up a local coordinate reference
    refx = np.min(elon);
    refy = np.min(nlat);
    elon = elon - refx;
    nlat = nlat - refy;

    # Setting calculation parameters
    EstimateRadius = float(MyParams.method_specific["estimateradiuskm"]) * 1000;  # convert to meters
    ns = int(MyParams.method_specific["nstations"]);  # number of selected stations

    # 2. The main loop, getting displacement gradients around stations
    # find at least 5 smallest distance stations
    Uxx = np.zeros((gy, gx));
    Uyy = np.zeros((gy, gx));
    Uxy = np.zeros((gy, gx));
    Uyx = np.zeros((gy, gx));
    exx = np.zeros((gy, gx));
    exy = np.zeros((gy, gx));
    eyy = np.zeros((gy, gx));
    rot = np.zeros((gy, gx));
    for i in range(gx):
        for j in range(gy):
            [gridX_loc, gridY_loc] = coord_to_local_utm(xlons[i], ylats[j], refx, refy);
            X = elon;
            Y = nlat;
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
                Uxx[j, i] = 0;
                Uyy[j, i] = 0;
                Uxy[j, i] = 0;
                Uyx[j, i] = 0;
            else:
                dx = np.array([[dU], [dxU], [dyU]]);
                dy = np.array([[dV], [dxV], [dyV]]);
                modelx = np.dot(np.linalg.inv(G), dx);
                modely = np.dot(np.linalg.inv(G), dy);
                Uxx[j, i] = modelx[1];
                Uyy[j, i] = modely[2];
                Uxy[j, i] = modelx[2];
                Uyx[j, i] = modely[1];

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

    return [xlons, ylats, rot, exx, exy, eyy];


def velfield_to_huang_format(myVelfield):
    elon, nlat = [], [];
    e, n, esig, nsig = [], [], [], [];
    for item in myVelfield:
        [x, y, _] = utm_conversion.deg2utm([item.nlat], [item.elon]);
        elon.append(x);
        nlat.append(y);
        e.append(item.e*0.001);
        n.append(item.n*0.001);
        esig.append(item.se*0.001);
        nsig.append(item.sn*0.001);

    return [np.array(elon), np.array(nlat), np.array(e), np.array(n), np.array(esig), np.array(nsig)];

def coord_to_local_utm(lon, lat, utm_xref, utm_yref):
    [x, y, _] = utm_conversion.deg2utm([lat], [lon]);
    local_utmx = x - utm_xref;
    local_utmy = y - utm_yref;
    return [local_utmx, local_utmy];
