# Strain calculate - Velmap
# From:  Wang, H., & Wright, T. J. (2012). 
# AUthor: Yi-Chieh Lee
# Last modified: 02.07.2023
# Note: Only G_gps_plane and only consider Ve & Vn for now

import numpy as np
from scipy.sparse import block_diag
from .. import velocity_io, strain_tensor_toolbox, utilities
from strain.utilities import getVels
from strain.models.strain_2d import Strain_2d
import sys

class velmap(Strain_2d):
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir);
        self._Name = 'Velmap'
        self._tempdir = params.outdir;
        self._nrows, self._ncols, self._smoothing_constant, self._grid_size_lon, self._grid_size_lat = verify_inputs_velmap(params.method_specific);

    def compute(self, myVelfield):
        [Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd] = compute_velmap(myVelfield, self._range_strain,
                                                                          self._nrows, self._ncols, self._smoothing_constant,
                                                                          self._grid_size_lon, self._grid_size_lat);
        # Report observed and residual velocities within bounding box
        velfield_within_box = utilities.filter_by_bounding_box(myVelfield, self._strain_range);
        model_velfield = utilities.create_model_velfield(self._xdata, self._ydata, Ve, Vn, velfield_within_box);
        residual_velfield = utilities.subtract_two_velfields(velfield_within_box, model_velfield);
        return [Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd, velfield_within_box, residual_velfield];
   

def compute_velmap(myVelfield, range_strain, nrows, ncols, smoothing_constant, grid_size_lon, grid_size_lat):
    '''Compute the interpolated velocity field'''

    # Read File
    dlon, dlat, ve, vn, ve, se, sn, su = getVels(myVelfield)
    mask = (dlon > utilities.get_float_range(range_strain)[0] & (dlon < utilities.get_float_range(range_strain)[1]) & (dlat > utilities.get_float_range(range_strain)[2]) & (dlat < utilities.get_float_range(range_strain)[3]))
    ve = ve[mask], vn = vn[mask], vu = vu[mask]
    se = se[mask], se = se[mask], su = su[mask]

    # GPS_plane
    lon = np.repeat(np.linspace(utilities.get_float_range(range_strain)[1], utilities.get_float_range(range_strain)[0], num = nrows), ncols).reshape(-1, 1)
    lat = np.tile(np.linspace(utilities.get_float_range(range_strain)[2], utilities.get_float_range(range_strain)[3], num = ncols), nrows).reshape(-1, 1)
    xy_gps = np.hstack((lon, lat))
    G_plane_gps = np.hstack((xy_gps, np.ones((xy_gps.shape[0], 1))))
    G_plane = np.block([ [G_plane_gps, np.zeros((G_plane_gps.shape[0], G_plane_gps.shape[1]))], [np.zeros((G_plane_gps.shape[0], G_plane_gps.shape[1])), G_plane_gps] ])

    # G-gps matrix
    Ngps = len(ve)
    G_gps_x = np.zeros((Ngps, nrows*ncols))
    G_gps_y = np.zeros((Ngps, nrows*ncols))
    #G_gps_z = np.zeros((Ngps, nrows*ncols))
    

    for i in range(Ngps):
        for j in range(ncols*nrows):
            if (np.round(dlon[i], 2) == G_plane_gps[j, 0]) & (np.round(dlat[i], 2) == G_plane_gps[j, 1]):
                index = i

        G_gps_x[i, j] = 1
        G_gps_y[i, j] = 1
        #G_gps_z[i, j] = 1
    
    G_gps = block_diag((G_gps_x, G_gps_y)).toarray()
    #G_gps = block_diag((G_gps_x, G_gps_y, G_gps_z)).toarray()

    # Laplacian smoothing matrix 
    Lap = Laplacian_backslip(nrows, ncols, grid_size_lon, grid_size_lat, 0)
    full_Lap = np.block([[Lap, np.zeros((ncols*nrows, ncols*nrows))], [np.zeros((ncols*nrows, ncols*nrows)), Lap]])

    # Uncertainities
    SIG_gps = np.append(se, sn)
    #SIG_gps = np.append(se, sn, su)
    COV_gps = np.diag(SIG_gps**2)


    # Solving inverse problem
    d = np.concatenate((ve, vn, np.zeros((full_Lap.shape[0]))))

    G = np.vstack(( np.hstack((G_gps, np.zeros((G_gps.shape[0], G_plane.shape[1])))), np.hstack((full_Lap, np.zeros((full_Lap.shape[0], G_plane.shape[1])))) ))
    # should add Gmatrix for ENU

    L2SIG = block_diag((COV_gps, smoothing_constant * np.eye(full_Lap.shape[0]))).toarray()

    mask3 = np.isnan(d) | np.isnan(np.sum(G, axis=1))
    d = d[~mask3]
    G = G[~mask3, :]

    L2SIG = L2SIG[~mask3, :]
    L2SIG = L2SIG[:, ~mask3]

    try:
        Lsig = np.linalg.cholesky(L2SIG)
    except np.linalg.LinAlgError as e:
        print("Cholesky decomposition failed:", e)
        sys.exit()

    dd = np.linalg.solve(Lsig, d)
    GG = np.linalg.solve(Lsig, G)

    vhat, resid, rank, s = np.linalg.lstsq(GG, dd, rcond=None)

    Nc = (len(vhat) - len(G_plane[0])) // 2
    Ve_pred = vhat[:Nc]
    Vn_pred = vhat[Nc:2*Nc]  

    Ve = Ve_pred.reshape(ncols, nrows)
    Vn = Vn_pred.reshape(ncols, nrows)

    # Calculate strain rate
    exx, eyy, exy, rot = strain_tensor_toolbox.strain_on_regular_grid(ncols, nrows, Ve, Vn)

    # Return the strain rates etc. in the same units as other methods
    return Ve, Vn, rot*1000, exx*1000, exy*1000, eyy*1000



def verify_inputs_velmap(method_specific_dict):
    if 'nrows' not in method_specific_dict.keys():
        raise ValueError("\nvelmap requires the number of row. Please add to method_specific config. Exiting.\n");
    if 'ncols' not in method_specific_dict.keys():
        raise ValueError("\nvelmap requires the number of column. Please add to method_specific config. Exiting.\n");
    if 'smoothing_constrant' not in method_specific_dict.keys():
        raise ValueError("\nvelmap requires the value of smoothing constrain. Please add to method_specific config. Exiting.\n");
    if 'grid_size_lon' not in method_specific_dict.keys():
        raise ValueError("\nvelmap requires the grid size of longitude. Please add to method_specific config. Exiting.\n");
    if 'grid_size_lat' not in method_specific_dict.keys():
        raise ValueError("\nvelmap requires the grid size of latitude. Please add to method_specific config. Exiting.\n");
    
    nrows = method_specific_dict["nrows"];
    ncols = method_specific_dict["ncols"];
    smoothing_constrant = method_specific_dict["smoothing_constrant"];
    grid_size_lon = method_specific_dict["grid_size_lon"];
    grid_size_lat = method_specific_dict["grid_size_lat"];

    return nrows, ncols, smoothing_constrant, grid_size_lon, grid_size_lat;


def Laplacian_backslip(nve, nhe, delx, dely, surf):
    ngrid = nhe * nve
    Lap = np.zeros([ngrid, ngrid])
    xpartd = np.zeros([ngrid, ngrid])
    ypartd = np.zeros([ngrid, ngrid])
    temp = np.zeros([nve, nve])

    # x-derivative for central part of grid (exclude left & right edges)
    for i in range(nve, (nhe*nve-nve)): 
        xpartd[i, i-nve] = 1
        xpartd[i, i] = -2
        xpartd[i, i+nve] = 1

    # y-derivative for central part of grid (exclude top & bottom edges)
    for i in range(1, nve-1):
        temp[i, i-1:i+2] = [1, -2, 1]

    for j in range(nhe):
        k = (j)*nve
        ypartd[k:k+nve, k:k+nve] = temp

    # need to do top edge y-derivative for slip breaking the surface
    for i in range(nhe):
        k = (i-1)*nve+1

        if surf == 1:
            ypartd[k, k:k+2] = [-1, 1]
        else:
            ypartd[k, k:k+2] = [-2, 1]

    # bottom edge y-derivative - smooth to zero
    for i in range(nhe):
        k = i*nve
        if k == 0 :
            ypartd[0, 0:2] = [-2, 1]
        else:
            ypartd[k, k-1:k+1] = [1, -2]
        
    Lap = xpartd/delx**2 + ypartd/dely**2
    # Lap_inv = np.linalg.inv(Lap)
    return Lap

