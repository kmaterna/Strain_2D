# Use GPS Gridder from GMT to interpolate between GPS stations
# The algorithm is based on the greens functions for elastic sheets with a given Poisson's ratio. 
# From: Sandwell, D. T., and P. Wessel (2016),
# Interpolation of 2-D vector data using constraints from elasticity, Geophys. Res.Lett. 


import numpy as np
import subprocess
import xarray as xr
import os
import shutil
from strain.strain_tensor_toolbox import strain_on_regular_grid
from strain.velocity_io import write_gmt_format
from strain.utilities import (
    make_grid, filter_by_bounding_box, create_model_velfield,
    subtract_two_velfields, get_string_range, get_string_inc,

)
from strain.models.strain_2d import Strain_2d


class gpsgridder(Strain_2d):
    """ gps_gridder class for 2d strain rate """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.xdata, params.ydata,
                           params.outdir)
        self._Name = 'gpsgridder'
        self._tempdir = params.outdir
        self._poisson, self._fd, self._eigenvalue = verify_inputs_gpsgridder(params.method_specific)

    def compute(self, myVelfield):
        [Ve, Vn, rot_grd, exx_grd, exy_grd, eyy_grd] = compute_gpsgridder(
            myVelfield, self._strain_range,
            self._grid_inc, self._poisson, self._fd,
            self._eigenvalue, self._tempdir
        )

        ########################################################
        # check that the grids are not off-by-one
        reduce_x, reduce_y = False, False
        lenx, leny = len(self._xdata), len(self._ydata)
        if Ve.shape[0] != leny:
            if Ve.shape[0] == leny - 1:
                reduce_y = True
            else:
                raise RuntimeError('Velocity is not the right shape')
        if Ve.shape[1] != lenx:
            if Ve.shape[1] == lenx + 1:
                reduce_x = True
            else:
                raise RuntimeError('Velocity is not the right shape')

        if reduce_x:
            Ve = Ve[:,:-1]
            Vn = Vn[:,:-1]
            rot_grd = rot_grd[:,:-1]
            exx_grd = exx_grd[:,:-1]
            exy_grd = exy_grd[:,:-1]
            eyy_grd = eyy_grd[:,:-1]
        if reduce_y:
            Ve = Ve[:-1,:]
            Vn = Vn[:-1,:]
            rot_grd = rot_grd[:-1,:]
            exx_grd = exx_grd[:-1,:]
            exy_grd = exy_grd[:-1,:]
            eyy_grd = eyy_grd[:-1,:]
        ########################################################

        # Report observed and residual velocities within bounding box
        velfield_within_box = filter_by_bounding_box(myVelfield, self._strain_range)
        model_velfield = create_model_velfield(self._xdata, self._ydata, Ve, Vn, velfield_within_box)
        residual_velfield = subtract_two_velfields(velfield_within_box, model_velfield)
        return [Ve, Vn, np.empty(Ve.shape), np.empty(Vn.shape), rot_grd, exx_grd, exy_grd, eyy_grd, velfield_within_box, residual_velfield]


def verify_inputs_gpsgridder(method_specific_dict):
    if 'poisson' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires poisson's ratio. Please add to method_specific config. Exiting.\n")
    if 'fd' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires fudge factor fd. Please add to method_specific config. Exiting.\n")
    if 'eigenvalue' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires eigenvalue. Please add to method_specific config. Exiting.\n")
    poisson = method_specific_dict["poisson"]
    fd = method_specific_dict["fd"]
    eigenvalue = method_specific_dict["eigenvalue"]
    return poisson, fd, eigenvalue


# ----------------- COMPUTE -------------------------

def compute_gpsgridder(myVelfield, range_strain, inc, poisson, fd, eigenvalue, tempoutdir):
    print("------------------------------\nComputing strain via gpsgridder method.")
    write_gmt_format(myVelfield, "tempgps.txt")
    command = "gmt gpsgridder tempgps.txt" + \
              " -R" + get_string_range(range_strain, x_buffer=inc[0]/2, y_buffer=inc[1]/2) + \
              " -I" + get_string_inc(inc) + \
              " -S" + poisson + \
              " -Fd" + fd + \
              " -C" + eigenvalue + \
              " -Emisfitfile.txt -fg -r -Gnc_%s.nc"
    print(command)
    subprocess.call(command, shell=True)  # makes a netcdf grid file
    # -R = range. -I = interval. -E prints the model and data fits at the input stations (very useful).
    # -S = poisson's ratio. -Fd = fudge factor. -C = eigenvalues below this value will be ignored.
    # -fg = flat earth approximation. -G = output netcdf files (x and y displacements).
    # -r is pixel node registration
    # You should experiment with Fd and C values to find something that you like (good fit without overfitting).
    # For Northern California, I like -Fd0.01 -C0.005. -R-125/-121/38/42.2

    if os.path.isfile('tempgps.txt'):
        os.remove('tempgps.txt')
    if os.path.isfile('gmt.history'):
        os.remove('gmt.history')
    shutil.move(src='misfitfile.txt', dst=os.path.join(tempoutdir, 'misfitfile.txt'))
    shutil.move(src='nc_u.nc', dst=os.path.join(tempoutdir, 'nc_u.nc'))
    shutil.move(src='nc_v.nc', dst=os.path.join(tempoutdir, 'nc_v.nc'))

    # Get ready to do strain calculation.
    file1 = tempoutdir+"nc_u.nc"
    file2 = tempoutdir+"nc_v.nc"
    
    # gpsgridder won't give back the original grid in some cases
    lons, lats, _ = make_grid(range_strain, inc)
    
    ds = xr.open_dataset(file1)
    udata = ds["z"].interp(lat=lats,lon=lons, method='nearest', kwargs={"fill_value": "extrapolate"}).to_numpy()
    ds = xr.open_dataset(file2)
    vdata = ds["z"].interp(lat=lats,lon=lons, method='nearest', kwargs={"fill_value": "extrapolate"}).to_numpy()

    xinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $8}\'', shell=True))  # x-inc
    yinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $9}\'', shell=True))  # y-inc
    xinc = xinc * 111.000 * np.cos(np.deg2rad(range_strain[2]))  # in km (not degrees)
    yinc = yinc * 111.000   # in km (not degrees)

    [exx, eyy, exy, rot] = strain_on_regular_grid(xinc, yinc, udata * 1000, vdata * 1000)

    print("Success computing strain via gpsgridder method.\n")
    return [udata, vdata, abs(rot), exx, exy, eyy]
