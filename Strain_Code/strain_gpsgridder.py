# Use GPS Gridder from GMT to interpolate between GPS stations
# The algorithm is based on the greens functions for elastic sheets with a given Poisson's ratio. 
# From: Sandwell, D. T., and P. Wessel (2016),
# Interpolation of 2-D vector data using constraints from elasticity, Geophys. Res.Lett. 


import numpy as np
import subprocess
from Tectonic_Utils.read_write import netcdf_read_write
import strain_tensor_toolbox
import output_manager
import configure_functions


# ----------------- COMPUTE -------------------------
def compute(myVelfield, MyParams):
    print("------------------------------\nComputing strain via gpsgridder method.");
    output_manager.write_simple_gmt_format(myVelfield, "tempgps.txt");
    command = "gmt gpsgridder tempgps.txt" + \
              " -R"+configure_functions.get_string_range(MyParams.range_strain) + \
              " -I"+configure_functions.get_string_inc(MyParams.inc) + \
              " -S"+MyParams.method_specific['poisson'] + \
              " -Fd"+MyParams.method_specific['fd'] + \
              " -C"+MyParams.method_specific['eigenvalue'] + \
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
    subprocess.call(['mv', 'misfitfile.txt', MyParams.outdir], shell=False);
    subprocess.call(['mv', 'nc_u.nc', MyParams.outdir], shell=False);
    subprocess.call(['mv', 'nc_v.nc', MyParams.outdir], shell=False);

    # Get ready to do strain calculation.
    file1 = MyParams.outdir+"nc_u.nc";
    file2 = MyParams.outdir+"nc_v.nc";
    [xdata, ydata, udata] = netcdf_read_write.read_any_grd(file1);
    [_, _, vdata] = netcdf_read_write.read_any_grd(file2);
    xinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $8}\'', shell=True));  # x-inc
    yinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $9}\'', shell=True));  # y-inc
    xinc = xinc * 111.000 * np.cos(np.deg2rad(MyParams.range_strain[2]));  # in km (not degrees)
    yinc = yinc * 111.000;   # in km (not degrees)
    [ydim, xdim] = np.shape(udata)
    exx = np.zeros(np.shape(vdata));
    exy = np.zeros(np.shape(vdata));
    eyy = np.zeros(np.shape(vdata));
    rot = np.zeros(np.shape(vdata));  # 2nd invariant of rotation rate tensor

    # the strain calculation
    for j in range(ydim-1):
        for i in range(xdim-1):
            up = udata[j][i];
            vp = vdata[j][i];
            uq = udata[j][i+1];
            vq = vdata[j][i+1];
            ur = udata[j+1][i];
            vr = vdata[j+1][i];

            [dudx, dvdx, dudy, dvdy] = strain_tensor_toolbox.compute_displacement_gradients(up, vp, ur, vr, uq, vq,
                                                                                            xinc, yinc);

            # The basic strain tensor components (units: nanostrain per year)
            [exx1, exy1, eyy1, rot1] = strain_tensor_toolbox.compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy);
            rot[j][i] = abs(rot1);
            exx[j][i] = exx1;
            exy[j][i] = exy1;
            eyy[j][i] = eyy1;

    print("Success computing strain via gpsgridder method.\n");

    return [xdata, ydata, rot, exx, exy, eyy];
