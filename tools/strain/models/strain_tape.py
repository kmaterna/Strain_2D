import numpy as np
import gps_io_functions
from Tectonic_Utils.read_write import netcdf_read_write
import scipy.interpolate as interp
from .. import strain_tensor_toolbox


# 2021: THIS IS NOT UPDATED WITH NEW ARCHITECTURE OF STRAIN LIBRARY.

# This code was created to work with matlab scripts published on Github by Carl Tape under the name surfacevel2strain.
# The functions reformulate .txt files to be read by Tape's code.

# This function inputs a NAM velo file and selects the columns needed for the tape strain method, in correct order.
# These columns are: lon, lat, ve, vn, vu, se, sn, su, ren, reu, rnu, start, finish, name
def input_to_tape(filename):
    if "unr" in filename or "UNR" in filename or "MIDAS" in filename or "midas" in filename:
        [velfield] = gps_io_functions.read_unr_vel_file(filename)
        infile = np.array((velfield[2], velfield[1], velfield[4], velfield[3], velfield[5], velfield[7], velfield[6], velfield[8]))
        name_field = velfield[0]

        lon = np.array(infile[0])
        for i in range(len(lon)):
            if lon[i] >= 180:
                lon[i] = lon[i] - 360

        starts = []
        ends = []
        for i in range(len(velfield[9])):
            start = float(velfield[9][i].date().strftime("%Y%m%d"))
            end = float(velfield[10][i].date().strftime("%Y%m%d"))
            # start = datetime.strptime(velfield[9], %Y%m%d)
            # end = datetime.strptime(velfield[10], %Y%m%d)
            starts.append(start)
            ends.append(end)

        infile = np.vstack((lon, np.array(infile[1]), np.array(infile[2]), np.array(infile[3]), np.array(infile[4]), np.array(infile[5]), np.array(infile[6]), np.array(infile[7]), np.zeros(len(lon)), np.zeros(len(lon)), np.zeros(len(lon)), starts, ends))

    elif "nam" in filename or "NAM" in filename:
        infile = np.loadtxt(filename, skiprows = 37, usecols = (8, 7, 20, 19, 21, 23, 22, 24, 25, 27, 26, 28, 29), unpack = True)

        temp_file = open(filename, 'r')
        data = temp_file.readlines()[37:]
        temp_file.close()

        name_field = []
        for line in data:
            name_field.append(line.split()[0])

        lon = np.array(infile[0])
        for i in range(len(lon)):
            if lon[i] >= 180:
                lon[i] = lon[i] - 360
        infile = np.vstack(lon, infile[1:])

        infile[2:10] = infile[2:10]*1e3

    else:
        print("Cannot read file format")

    outfile = np.vstack((infile, name_field))
    outfile = np.transpose(outfile)

    return outfile


def drive_tape(myParams):
    # generic steps
    indir = "../compearth/surfacevel2strain/matlab_output/"
    print("Producing gridded dataset of: ")
    x, y, tt, tp, pp = input_tape(indir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat",
                                  "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
    [I2nd, max_shear, _, e1, e2, v00, v01, v10, v11, dilatation] = compute_tape(tt, tp, pp);
    # second invariant
    newx, newy, newI2nd = nn_interp(x, y, I2nd, myParams.coord_box, myParams.grid_inc);
    output_tape(newx, newy, newI2nd, myParams.outdir, "I2nd.nc");
    # dilatation
    newx, newy, newdila = nn_interp(x, y, dilatation, myParams.coord_box, myParams.grid_inc);
    output_tape(newx, newy, newdila, myParams.outdir, "dila.nc");
    # max shear
    newx, newy, newmax = nn_interp(x, y, max_shear, myParams.coord_box, myParams.grid_inc);
    output_tape(newx, newy, newmax, myParams.outdir, "max_shear.nc");
    # azimuth
    azimuth = strain_tensor_toolbox.max_shortening_azimuth(e1, e2, v00, v01, v10, v11)
    newx, newy, newaz = nn_interp(x, y, azimuth, myParams.coord_box, myParams.grid_inc);
    netcdf_read_write.produce_output_netcdf(newx, newy, newaz, 'degrees', myParams.outdir+'azimuth.nc');
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


# This function computes second invariant, max shear strain, etc from symmetric strain tensor components
# (in spherical coords, from tape)
def compute_tape(thth, thph, phph):
    I2nd, max_shear, rot = [], [], [];
    e1, e2 = [], [];  # eigenvalues
    v00, v01, v10, v11 = [], [], [], [];  # eigenvectors
    dilatation = [];  # dilatation	= e1+e2

    for i in range(len(thth)):
        eyy = 1e9*phph[i]
        exy = -1e9*thph[i]
        exx = 1e9*thth[i]
        [e11, e22, v] = strain_tensor_toolbox.eigenvector_eigenvalue(exx, exy, eyy);
        e1.append(e11)
        e2.append(e22)
        v00.append(v[0][0])
        v01.append(v[0][1])
        v10.append(v[1][0])
        v11.append(v[1][1])
        max_shear.append(abs((-e11 + e22)/2))
        dilatation.append(e11+e22)
        I2nd.append(np.log10(np.abs(strain_tensor_toolbox.second_invariant(exx, exy, eyy))))
    return [I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation]


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


# outputs a netcdf to the desired results directory
def output_tape(lon, lat, vals, outdir, file):
    netcdf_read_write.produce_output_netcdf(lon, lat, vals, 'per yr', outdir+file);
    print("Success fitting wavelet-generated data to the required grid!");
    return;




# Takes reformatted data and outputs it as a .txt for use in matlab scripts.
# Outdir should refer to location accessed by matlab scripts, and outdir and outfile should be strings.
def output_to_tape(data, outdir, outfile):
    print("Creating %s in %s" % (outfile, outdir))
    np.savetxt(outdir+outfile, data, delimiter=" ", fmt="%s")
    return


# for PBO/NAM08:
# infile = input_to_tape("../Example_data/NAM08_pbovelfile_feb2018.vel")
# output_to_tape(infile, "../../compearth/surfacevel2strain/data/", "NAM08.txt")

# for UNR:
infile = input_to_tape("../Example_data/midas.NA12.txt")
output_to_tape(infile, "../../compearth/surfacevel2strain/data/", "UNR.txt")
