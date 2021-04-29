# Use GPS Gridder from GMT to interpolate between GPS stations
# The algorithm is based on the greens functions for elastic sheets with a given Poisson's ratio. 
# From: Sandwell, D. T., and P. Wessel (2016),
# Interpolation of 2-D vector data using constraints from elasticity, Geophys. Res.Lett. 


import numpy as np
import subprocess
from Tectonic_Utils.read_write import netcdf_read_write
from .. import velocity_io, strain_tensor_toolbox, utilities
from strain.models.strain_2d import Strain_2d


class gpsgridder(Strain_2d):
    """ gps_gridder class for 2d strain rate """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir);
        self._Name = 'gpsgridder'
        self._tempdir = params.outdir;
        self._poisson, self._fd, self._eigenvalue = verify_inputs_gpsgridder(params.method_specific);

    def compute(self, myVelfield):
        [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd] = compute_gpsgridder(myVelfield, self._strain_range,
                                                                              self._grid_inc, self._poisson, self._fd,
                                                                              self._eigenvalue, self._tempdir);
        return [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd];


def verify_inputs_gpsgridder(method_specific_dict):
    if 'poisson' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires poisson's ratio. Please add to method_specific config. Exiting.\n");
    if 'fd' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires fudge factor fd. Please add to method_specific config. Exiting.\n");
    if 'eigenvalue' not in method_specific_dict.keys():
        raise ValueError("\ngps_gridder requires eigenvalue. Please add to method_specific config. Exiting.\n");
    poisson = method_specific_dict["poisson"];
    fd = method_specific_dict["fd"];
    eigenvalue = method_specific_dict["eigenvalue"];
    return poisson, fd, eigenvalue;

# ----------------- COMPUTE -------------------------
def compute_gpsgridder(myVelfield, range_strain, inc, poisson, fd, eigenvalue, tempoutdir):
    print("------------------------------\nComputing strain via gpsgridder method.");
    velocity_io.write_gmt_format(myVelfield, "tempgps.txt");
    command = "gmt gpsgridder tempgps.txt" + \
              " -R" + utilities.get_string_range(range_strain, x_buffer=0.02, y_buffer=0.02) + \
              " -I" + utilities.get_string_inc(inc) + \
              " -S" + poisson + \
              " -Fd" + fd + \
              " -C" + eigenvalue + \
              " -Emisfitfile.txt -fg -r -Gnc_%s.nc";
    print(command);
    subprocess.call(command, shell=True);  # makes a netcdf grid file
    # -R = range. -I = interval. -E prints the model and data fits at the input stations (very useful).
    # -S = poisson's ratio. -Fd = fudge factor. -C = eigenvalues below this value will be ignored.
    # -fg = flat earth approximation. -G = output netcdf files (x and y displacements).
    # You should experiment with Fd and C values to find something that you like (good fit without overfitting).
    # For Northern California, I like -Fd0.01 -C0.005. -R-125/-121/38/42.2

    subprocess.call(['rm', 'tempgps.txt'], shell=False);
    subprocess.call(['rm', 'gmt.history'], shell=False);
    subprocess.call(['mv', 'misfitfile.txt', tempoutdir], shell=False);
    subprocess.call(['mv', 'nc_u.nc', tempoutdir], shell=False);
    subprocess.call(['mv', 'nc_v.nc', tempoutdir], shell=False);

    # Get ready to do strain calculation.
    file1 = tempoutdir+"nc_u.nc";
    file2 = tempoutdir+"nc_v.nc";
    [xdata, ydata, udata] = netcdf_read_write.read_any_grd(file1);
    [_, _, vdata] = netcdf_read_write.read_any_grd(file2);
    udata = udata.T;
    vdata = vdata.T;

    xinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $8}\'', shell=True));  # x-inc
    yinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $9}\'', shell=True));  # y-inc
    xinc = xinc * 111.000 * np.cos(np.deg2rad(range_strain[2]));  # in km (not degrees)
    yinc = yinc * 111.000;   # in km (not degrees)

    [exx, eyy, exy, rot] = strain_tensor_toolbox.strain_on_regular_grid(xinc, yinc, udata, vdata)
    # Lastly we multiply by 1000 for units
    exx = exx.T;
    eyy = eyy.T;
    exy = exy.T;
    rot = rot.T;
    exx = np.multiply(exx, 1000);
    exy = np.multiply(exy, 1000);
    eyy = np.multiply(eyy, 1000);
    rot = abs(np.multiply(rot, 1000));

    print("Success computing strain via gpsgridder method.\n");
    return [xdata, ydata, rot, exx, exy, eyy];
