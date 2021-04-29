import numpy as np
import scipy.interpolate as interp
import subprocess
from strain.models.strain_2d import Strain_2d


class tape(Strain_2d):
    """ Tape class for 2d strain rate, with general strain_2d behavior """
    def __init__(self, params):
        Strain_2d.__init__(self, params.inc, params.range_strain, params.range_data, params.outdir);
        self._Name = 'tape'
        self._code_dir, self._qmin, self._qmax, self._qsec = verify_inputs_tape(params.method_specific);

    def compute(self, myVelfield):
        # Setup for Matlab calculation
        configure_file = self._outdir + "surfacevel2strain_config_params.txt"
        velocity_file = self._outdir + "vel_tape.txt";
        subprocess.call(['mkdir', '-p', self._code_dir+'/matlab_output'], shell=False);
        write_to_tape_vel_format(myVelfield, velocity_file);
        write_tape_parameter_file(self._data_range, self._code_dir, self._qmin, self._qmax, self._qsec, velocity_file, configure_file);

        # Now go away and do your matlab calculation. Come back with the output stem. 
        output_tag = input("\n\nNow go away, run Matlab using the instructions and parameters in "+configure_file+"\n\n. When done, enter the output stem of your favorite model: ");
        # Format is like: /Users/kzm/Documents/Software//compearth/surfacevel2strain/matlab_output/_d-01_q04_q07_b1_2D_s1_u1
        
        # Parse the results
        x, y, tt, tp, pp, rot = input_tape(output_tag+"_strain.dat", output_tag+"_Dtensor_6entries.dat", output_tag+"_Wtensor_3entries.dat");
        exx, exy, eyy, rot = compute_tape(tt, tp, pp, rot);
        lons, lats, exx = nn_interp(x, y, exx, self._strain_range, self._grid_inc);
        _, _, exy = nn_interp(x, y, exy, self._strain_range, self._grid_inc);
        _, _, eyy = nn_interp(x, y, eyy, self._strain_range, self._grid_inc);
        _, _, rot = nn_interp(x, y, rot, self._strain_range, self._grid_inc); 
        return [lons, lats, rot, exx, exy, eyy];


def verify_inputs_tape(method_specific_dict):
    # Takes a dictionary and verifies that it contains the right parameters for Tape method
    if 'code_dir' not in method_specific_dict.keys():
        raise ValueError("\nTape requires code_dir. Please add to method_specific config. Exiting.\n");
    if 'qmin' not in method_specific_dict.keys():
        raise ValueError("\nTape requires qmin. Please add to method_specific config. Exiting.\n"); 
    if 'qmax' not in method_specific_dict.keys():
        raise ValueError("\nTape requires qmax. Please add to method_specific config. Exiting.\n");
    if 'qsec' not in method_specific_dict.keys():
        raise ValueError("\nTape requires qsec. Please add to method_specific config. Exiting.\n"); 
    code_dir = method_specific_dict["code_dir"];
    qmin = method_specific_dict["qmin"];
    qmax = method_specific_dict["qmax"];
    qsec = method_specific_dict["qsec"];
    return code_dir, qmin, qmax, qsec;


def write_to_tape_vel_format(velfield, outfile):
    """
    Interface with matlab scripts published on Github by Carl Tape under the name surfacevel2strain.
    The Tape-format columns are: lon, lat, ve, vn, vu, se, sn, su, ren, reu, rnu, start, finish, name
    ren, reu, rnu, start, finish, and name are not used
    """
    print("Writing file %s" % outfile);
    ofile = open(outfile, 'w');
    for item in velfield:
        ofile.write("%f %f %f %f %f %f %f %f 0 0 0 0 0 name\n" % (item.elon, item.nlat, item.e, item.n, item.u,
                                                                  item.se, item.sn, item.su) );
    return;


def write_tape_parameter_file(range_data, code_dir, qmin, qmax, qsec, velocity_file, outfile):
    """Write the file that tells you how to operate Tape's compearth code"""
    print("Writing parameter file %s " % outfile);
    ofile = open(outfile, 'w');
    ofile.write("Total steps for computing strain with Tape method:\n");
    ofile.write("Manual:   Create directory called compearth/ somewhere in your Software directory.\n");
    ofile.write("Manual:   Inside compearth, git clone Tape's surfacevel2strain repository.\n");
    ofile.write("Manual:   Open matlab.  \n");
    ofile.write("Manual:   In matlab, navigate to your strain experiment directory. \n");
    ofile.write("Manual:   In the matlab prompt, \n");
    ofile.write("Manual:   Run: >> addpath('" + code_dir + "/matlab')\n");
    ofile.write("Manual:   Run: >> setenv('REPOS','" + code_dir.split('compearth')[0] + "')   # Or wherever is your software directory, where compearth lives \n");
    ofile.write("Manual:   Run: >> surfacevel2strain\n");
    ofile.write("Manual:   Try models until you're happy. Start with the parameters below.\n");
    ofile.write("Manual:   Record your favorite parameters in the master config file for good measure.\n");
    ofile.write("Manual:   Collect name of output_stem of your favorite model (printed near the end of outputs from surfacevel2strain) and return to Python.\n");
    ofile.write("\n\n");

    ofile.write('MANUAL OPTIONS FOR INPUT() IN TAPE MATLAB PROGRAM\n');
    ofile.write('1  # new calculation\n');
    ofile.write('1  # \n');
    ofile.write('1  # \n');
    ofile.write('2  # 2d velocities\n');
    ofile.write('0  # no mask\n');
    ofile.write('1  # yes gmt outputs\n');
    ofile.write('-1 # ropt = overridden for user-specified velocity file\n');
    ofile.write('-1 # dopt = overridden for user-specified velocity file\n');
    # ofile.write('\'/Users/kzm/Documents/Projects/Strain/output/tape/vel_tape.txt\'\n');
    ofile.write('\''+velocity_file+'\'\n');
    ofile.write('\''+str(range_data[0])+' '+str(range_data[1])+' '+str(range_data[2])+' '+str(range_data[3])+'\'\n');
    ofile.write('1  # remove rotation\n');
    ofile.write(str(qmin) + '  # PLEASE_CHOOSE qmin\n');
    ofile.write(str(qmax) + '  # PLEASE_CHOOSE, qmax\n');
    ofile.write(str(qsec) + '  # PLEASE_CHOOSE, qsec\n');
    ofile.write('-1 # default lambda\n');
    ofile.write('-1 # default lambda\n');
    ofile.write('0  # rotation\n');
    ofile.close();
    return;


def input_tape(coordsfile, datafile, wfile):
    """
    Reads .dat files outputted by Carl Tape's matlab code, surfacevel2strain
    First, run tape code, selecting to output gmt files.
    this code inputs the "strain" .dat file and the "D tensor 6 entries" .dat and the "W 3 entries" .dat
    and outputs the coordinates and strain tensor components.    
    """
    incoords = np.loadtxt(coordsfile, usecols=(0, 1), unpack=True)
    lon = incoords[0]
    lat = incoords[1]
    infile = np.loadtxt(datafile, skiprows=1, usecols=(3, 4, 5), unpack=True)
    thth = infile[0]
    thph = infile[1]
    phph = infile[2]
    infile = np.loadtxt(wfile, skiprows=1, usecols=(0, 1, 2), unpack=True) 
    rot = infile[0]+infile[1]+infile[2];
    return lon, lat, thth, thph, phph, rot


def compute_tape(thth, thph, phph, rot_sph):
    """
    Computes symmetric strain tensor components
    (in spherical coords, from tape)
    """
    exx, eyy, exy, rot = [], [], [], [];
    for i in range(len(thth)):
        eyy.append(1e9*thth[i])
        exy.append(-1e9*thph[i])
        exx.append(1e9*phph[i])
        rot.append(1e9*rot_sph[i])
    return exx, exy, eyy, rot


def nn_interp(x, y, vals, coord_box, inc):
    """
    Performs scipy nearest-neighbor interpolation on the data
    Assumes Tape scripts were run on a finer grid (try npts = 250)
    the mins, maxes, and incriment should match that of other methods for easy comparison.
    """
    xmin, xmax = coord_box[0], coord_box[1];
    ymin, ymax = coord_box[2], coord_box[3];
    newx = np.arange(xmin, xmax + inc[0]/2, inc[0]);  # needed to match other grids somehow
    newy = np.arange(ymin, ymax + inc[1], inc[1])
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


"""
Steps for Tape:
Create compearth/ somewhere in your Software directory
Inside compearth, git clone Tape's surfacevel2strain
Configure and run a strain calculation. Some default options will helpfully be printed for you. 
Open matlab. 
>> setenv('REPOS','/Users/kzm/Documents/Software')   # where compearth lives 
>> addpath(self._code_dir+'/matlab');

If using python (No, BAD, this while thing didn't work because of matlab's input() function): 
# Set up matlab python integration with python3.7 in a pygmt environment (see computer setups)
go to "$matlabroot/extern/engines/python"
python setup.py install
python: 
>> import matlab.engine
>> eng = matlab.engine.start_matlab()
>> eng.setenv('REPOS',self._code_dir.split('compearth')[0], nargout=0);   # the Software directory
>> eng.addpath(self._code_dir+'/matlab', nargout=0)    # matlab folder within source directory
>> eng.surfacevel2strain(nargout=0)    # start the entire calculation loop. 
"""
