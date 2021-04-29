# Use the interpolation scheme of Zheng-Kang Shen et al.
# Shen, Z.-K., M. Wang, Y. Zeng, and F. Wang, Strain determination using spatially discrete geodetic data, 
# Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015.
# http://scec.ess.ucla.edu/~zshen/visr/visr.html
# The fortran files must be compiled and linked like this: 
# gfortran -c voronoi_area_version.f90
# gfortran visr.f voronoi_area_version.o -o visr.exe


import numpy as np
from .. import produce_gridded
from strain.models.strain_2d import Strain_2d
import subprocess, sys, os


class visr(Strain_2d):
    """ Visr class for 2d strain rate, with general strain_2d behavior """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir);
        self._Name = 'visr';
        self._tempdir = params.outdir;
        self._distwgt, self._spatwgt, self._smoothincs, self._exec = verify_inputs_visr(params.method_specific);

    def compute(self, myVelfield):
        [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd] = compute_visr(myVelfield, self._strain_range, self._grid_inc,
                                                                        self._distwgt, self._spatwgt,
                                                                        self._smoothincs, self._exec, self._tempdir);
        return [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd];


def verify_inputs_visr(method_specific_dict):
    if 'distance_weighting' not in method_specific_dict.keys():
        raise ValueError("\nvisr requires distance weighting. Please add to method_specific config. Exiting.\n");
    if 'spatial_weighting' not in method_specific_dict.keys():
        raise ValueError("\nvisr requires spatial weighting. Please add to method_specific config. Exiting.\n");
    if 'min_max_inc_smooth' not in method_specific_dict.keys():
        raise ValueError("\nvisr requires smoothing information. Please add to method_specific config. Exiting.\n");
    if 'executable' not in method_specific_dict.keys():
        raise ValueError("\nvisr requires path to executable. Please add to method_specific config. Exiting.\n");
    distance_weighting = method_specific_dict["distance_weighting"];
    spatial_weighting = method_specific_dict["spatial_weighting"];
    min_max_inc_smooth = method_specific_dict["min_max_inc_smooth"];
    executable = method_specific_dict["executable"];
    return distance_weighting, spatial_weighting, min_max_inc_smooth, executable;


def compute_visr(myVelfield, strain_range, inc, distwgt, spatwgt, smoothincs, executable, tempdir):
    print("------------------------------\nComputing strain via Visr method.");
    strain_config_file = 'visr_strain.drv';
    strain_data_file = 'strain_input.txt';  # can only be 20 characters long bc fortran!
    strain_output_file = 'strain_output.txt';  # can only be 20 characters long bc fortran!
    write_fortran_config_file(strain_config_file, strain_data_file, strain_output_file, strain_range, inc, distwgt, spatwgt, smoothincs);
    write_fortran_data_file(strain_data_file, myVelfield);
    check_fortran_executable(executable);
    call_fortran_compute(strain_config_file, executable);

    # We convert that text file into grids, which we will write as GMT grd files.
    [xdata, ydata, rot, exx, exy, eyy] = make_output_grids_from_strain_out(strain_output_file, strain_range, inc);
    subprocess.call(['mv', strain_config_file, tempdir], shell=False);
    subprocess.call(['mv', strain_data_file, tempdir], shell=False);
    subprocess.call(['mv', strain_output_file, tempdir], shell=False);
    print("Success computing strain via Visr method.\n");
    return [xdata, ydata, rot, exx, exy, eyy];


def write_fortran_config_file(strain_config_file, strain_data_file, strain_output_file, range_strain, inc, distwgt, spatwgt, smoothincs):
    # The config file will have the following components.
    """
    visr/visr_drive_strain.drv contains:
    visr/velh.cmm4                             ! Station coordinate and velocity solution file
    visr/strain.out                            ! Strain rate output file
    1                                          ! distance weighting scheme: 1=gaussian, 2=quadratic
    2                                          ! spatial weighting scheme: 1=azimuth, 2=voronoi area
    1 100 1                                    ! minimum, maximum, and incremental spatial smoothing constants (km)
    24                                         ! weighting threshold Wt
    0.5                                        ! uncertainty threshold for reset
    3                                          ! function: 1=velocity compatibility checking; 2=velocity interpolation; 3=strain rate interpolation
    -122.5 -114.0 32.0 37.5 0.04 0.04          ! Lon_min, Lon_max, Lat_min, Lat_max, dLon, dLat
    0                                          ! number of creep faults
    crp.dat                                    ! creep fault data file
    """
    # Begin by parsing params
    if distwgt == 'quadratic':
        dist_code = '2';
    else:
        dist_code = '1';
    if spatwgt == 'azimuth':
        spac_code = '1';
    else:
        spac_code = '2';
    min_smooth = smoothincs.split('/')[0];
    max_smooth = smoothincs.split('/')[1];
    inc_smooth = smoothincs.split('/')[2];

    # Write the fortran-readable config file.
    ofile = open(strain_config_file, 'w');
    ofile.write(strain_data_file+'                              ! Station coordinate and velocity solution file\n');
    ofile.write(strain_output_file+'                             ! Strain rate output file\n');
    ofile.write(dist_code + '                                          ! distance weighting scheme: 1=gaussian, 2=quadratic\n');
    ofile.write(spac_code + '                                          ! spatial weighting scheme: 1=azimuth, 2=voronoi area\n');
    ofile.write(min_smooth+' '+max_smooth+' '+inc_smooth+'                                    ! minimum, maximum, and incremental spatial smoothing constants (km)\n');
    ofile.write('2                                          ! weighting threshold Wt\n');
    ofile.write('0.05                                       ! uncertainty threshold for reset\n');
    ofile.write('3                                          ! function: 1=velocity compatibility checking; 2=velocity interpolation; 3=strain rate interpolation\n');
    ofile.write(str(range_strain[0])+' '+str(range_strain[1])+' '+str(range_strain[2])+' '+str(range_strain[3])+' '+str(inc[0])+' '+str(inc[1])+'                ! Lon_min, Lon_max, Lat_min, Lat_max, dLon, dLat\n');
    ofile.write('0                                          ! number of creep faults\n');
    ofile.write('crp.dat                                    ! creep fault data file\n');
    ofile.close();
    return;


def write_fortran_data_file(data_file, Velfield):
    ofile = open(data_file, 'w');
    # 35    format(a8,2f10.4,2(f7.2,f5.2),f7.3)
    # 0102_GPS -119.2642   34.5655 -29.02 0.79  22.96 0.73  0.082     4   7.2  1994.4
    for item in Velfield:
        if len(item.name) == 4:
            ofile.write(item.name + "_GPS ");
        elif len(item.name) == 8:
            ofile.write(item.name + " ");
        else:
            ofile.write("         ");  # 8 characters plus space
        ofile.write("%9.4f %9.4f " % (item.elon, item.nlat) );
        se = np.min((item.se, 9.99));  # uncertainty can only have one digit before decimal place :(
        sn = np.min((item.sn, 9.99));  # uncertainty can only have one digit before decimal place :(
        ofile.write("%6.2f %4.2f %6.2f %4.2f %6.3f " % (item.e, se, item.n, sn, 0.001) );
        ofile.write("    5   2.1  2005.0\n");
    ofile.close();
    return;


def call_fortran_compute(config_file, executable):
    # Here we will call the strain compute function, using visr's fortran code.
    # It will output a large text file.
    print("Calling visr.exe fortran code to compute strain: ");
    print(executable + ' < ' + config_file);
    subprocess.call(executable + ' < '+config_file, shell=True);
    return;


def make_output_grids_from_strain_out(infile, range_strain, inc):
    x, y = [], [];
    rotation, exx, exy, eyy = [], [], [], [];
    ifile = open(infile, 'r');
    for line in ifile:
        temp = line.split();
        if 'index' in temp or 'longitude' in temp or 'deg' in temp:
            continue;
        else:
            x.append(float(temp[0]));
            y.append(float(temp[1]));
            rotation.append(float(line[53:60]));
            exx.append(float(temp[9]));
            exy.append(float(temp[11]));
            eyy.append(float(temp[13]));
    ifile.close();

    if len(set(x)) == 0 and len(set(y)) == 0:
        print("ERROR! No valid strains have been computed. Try again.")
        sys.exit(0);

    lons, lats, zgrid = produce_gridded.make_grid(range_strain, inc);

    # Loop through x and y lists, find the index of coordinates in xaxis and yaxis sets, place them into 2d arrays.
    rot_grd = np.zeros(np.shape(zgrid));
    exx_grd = np.zeros(np.shape(zgrid));
    exy_grd = np.zeros(np.shape(zgrid));
    eyy_grd = np.zeros(np.shape(zgrid));
    lons = np.round(lons, 6);
    lats = np.round(lats, 6);
    x = np.round(x, 6);
    y = np.round(y, 6);

    for i in range(len(x)):
        xindex = np.where(lons == x[i])[0];
        yindex = np.where(lats == y[i])[0];
        xindex = xindex[0];
        yindex = yindex[0];
        rot_grd[yindex][xindex] = rotation[i];
        exx_grd[yindex][xindex] = exx[i];
        exy_grd[yindex][xindex] = exy[i];
        eyy_grd[yindex][xindex] = eyy[i];

    return [lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd];


def check_fortran_executable(path_to_executable):
    if os.path.isfile(path_to_executable):
        print("VISR executable found at %s " % path_to_executable);
    else:
        raise FileNotFoundError("VISR executable not found on your system. Check config file for path.");
    return;
