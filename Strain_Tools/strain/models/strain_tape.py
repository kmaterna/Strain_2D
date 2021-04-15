import numpy as np
import scipy.interpolate as interp
from . import strain_2d


class tape(strain_2d.Strain_2d):
    """ Tape class for 2d strain rate, with general strain_2d behavior """
    def __init__(self, params):
        strain_2d.Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir);
        self._Name = 'tape'
        self._code_dir = verify_inputs_tape(params.method_specific);

    def compute(self, myVelfield):
        write_to_tape_vel_format(myVelfield, params.outdir+"/vel_tape.txt");
        sys.exit(0);
        # HERE YOU CALL THE CODE.
        # x, y, tt, tp, pp = input_tape(self._code_dir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat",
        #                               "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
        # [exx, exy, eyy, rot] = compute_tape(tt, tp, pp);
        # lons, lats, exx = nn_interp(x, y, exx, self._strain_range, self._grid_inc);
        # _, _, exy = nn_interp(x, y, exy, self._strain_range, self._grid_inc);
        # _, _, eyy = nn_interp(x, y, eyy, self._strain_range, self._grid_inc);
        return [lons, lats, rot, exx, exy, eyy];


def verify_inputs_tape(method_specific_dict):
    # Takes a dictionary and verifies that it contains the right parameters for Tape method
    if 'code_dir' not in method_specific_dict.keys():
        raise ValueError("\nTape requires code_dir. Please add to method_specific config. Exiting.\n");
    code_dir = method_specific_dict["code_dir"];
    return code_dir;


# This code was created to work with matlab scripts published on Github by Carl Tape under the name surfacevel2strain.
# The Tape-format columns are: lon, lat, ve, vn, vu, se, sn, su, ren, reu, rnu, start, finish, name
# ren, reu, rnu, start, finish, and name are not used
def write_to_tape_vel_format(velfield, outfile):
    print("Writing file %s" % outfile);
    ofile = open(outfile, 'w');
    for item in velfield:
        ofile.write("%f %f %f %f %f %f %f %f 0 0 0 0 0 name\n" % (item.elon, item.nlat, item.e, item.n, item.u,
                                                                  item.se, item.sn, item.su) );
    return;


# This code works with .dat files outputted by Carl Tape's matlab code, surfacevel2strain
# First, run tape code, selecting to output gmt files.
# this code inputs the "strain" .dat file and the "D tensor 6 entries" .dat
# and outputs the coordinates and strain tensor components.
def input_tape(indir, coordsfile, datafile):
    incoords = np.loadtxt(indir+coordsfile, usecols=(0, 1), unpack=True)
    lon = incoords[0]
    lat = incoords[1]
    infile = np.loadtxt(indir+datafile, skiprows=1, usecols=(3, 4, 5), unpack=True)
    thth = infile[0]
    thph = infile[1]
    phph = infile[2]
    return lon, lat, thth, thph, phph


# This function computes symmetric strain tensor components
# (in spherical coords, from tape)
def compute_tape(thth, thph, phph):
    exx, eyy, exy, rot = [], [], [], [];
    for i in range(len(thth)):
        eyy.append(1e9*phph[i])
        exy.append(-1e9*thph[i])
        exx.append(1e9*thth[i])
    return [exx, exy, eyy, rot]


# This function performs scipy nearest-neighbor interpolation on the data
# it assumes Tape scripts were run on a finer grid (try npts = 250)
# the mins, maxes, and incriment should match that of other methods for easy comparison.
def nn_interp(x, y, vals, coord_box, inc):
    xmin, xmax = coord_box[0], coord_box[1];
    ymin, ymax = coord_box[2], coord_box[3];
    newx = np.arange(xmin, xmax, inc[0])
    newy = np.arange(ymin, ymax, inc[1])
    tempvals = []

    nn_interpolator = interp.NearestNDInterpolator((x, y), vals)
    for i in range(len(newx)):
        for j in range(len(newy)):
            val = nn_interpolator(newx[i], newy[j])
            tempvals.append(val)

    newvals = []
    for i in np.arange(0, len(tempvals), len(newy)):
        newvals.append(tempvals[i:i+len(newy)])

    newvals = np.transpose(newvals)

    return newx, newy, newvals


# # outputs a netcdf to the desired results directory
# def output_tape(lon, lat, vals, outdir, file):
#     netcdf_read_write.produce_output_netcdf(lon, lat, vals, 'per yr', outdir+file);
#     print("Success fitting wavelet-generated data to the required grid!");
#     return;
# # Takes reformatted data and outputs it as a .txt for use in matlab scripts.
# # Outdir should refer to location accessed by matlab scripts, and outdir and outfile should be strings.
# def output_to_tape(data, outdir, outfile):
#     print("Creating %s in %s" % (outfile, outdir))
#     np.savetxt(outdir+outfile, data, delimiter=" ", fmt="%s")
#     return
# for PBO/NAM08:
# infile = input_to_tape("../Example_data/NAM08_pbovelfile_feb2018.vel")
# output_to_tape(infile, "../../compearth/surfacevel2strain/data/", "NAM08.txt")


"""
Steps for Tape:
Create compearth/ somewhere in your Software directory
Inside compearth, git clone Tape's surfacevel2strain
first time: open matlab. 
>> setenv('REPOS','/Users/kzm/Documents/Software')   # where compearth lives  // maybe not necessary in the end? not sure. 

# Set up matlab python integration with python3.7 in a pygmt environment (see computer setups)
go to "$matlabroot/extern/engines/python"
python setup.py install

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'/Users/kzm/Documents/Software/compearth/surfacevel2strain/matlab', nargout=0)  // matlab folder within source directory
eng.surfacevel2strain(nargout=0)


go to surfacevel2strain in matlab prompt
I needed to create matlab_outputs manually for this to work.
from code place in matlab: call surfacevel2strain
I needed to remove one "error" in line 809 of surfacevel2strain, maybe? 
"""






