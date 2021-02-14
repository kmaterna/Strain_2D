# Convert triangulation polygon values into gridded netcdf

import numpy as np
import matplotlib.path
import scipy.interpolate as interp
from . import strain_tensor_toolbox
from Tectonic_Utils.read_write import netcdf_read_write


def tri2grid(grid_inc, range_strain,  triangle_vertices, rot, exx, exy, eyy):
    # steps to bring delaunay 1-D quantities into the same 2-D form as the other methods
    lons, lats, grid = make_grid(range_strain, grid_inc);
    print("Producing gridded dataset of Exx")
    exx_grd = find_in_triangles(triangle_vertices, exx, lons, lats, grid);
    print("Producing gridded dataset of: Exy")
    exy_grd = find_in_triangles(triangle_vertices, exy, lons, lats, grid);
    print("Producing gridded dataset of: Eyy")
    eyy_grd = find_in_triangles(triangle_vertices, eyy, lons, lats, grid);
    print("Producing gridded dataset of: Rot")
    rot_grd = find_in_triangles(triangle_vertices, rot, lons, lats, grid);
    return lons, lats, rot_grd, exx_grd, exy_grd, eyy_grd;


# makes grid for delaunay
def make_grid(coordbox, inc):
    # coordbox is [float, float, float, float] [W, E, S, N]
    # inc is a float
    # return value is a 2d array of zeros
    lonmin = coordbox[0]
    lonmax = coordbox[1]
    latmin = coordbox[2]
    latmax = coordbox[3]
    lons = np.arange(lonmin, lonmax+0.00001, inc[0])
    lats = np.arange(latmin, latmax+0.00001, inc[1])
    grid = np.zeros((len(lats), len(lons)));
    return lons, lats, grid


# searches path created by triangle vertices for each gridpoint, then assigns that triangle's value to the gridpoint
def find_in_triangles(triangles, values, lons, lats, grid):
    val_arr = np.nan*np.ones(np.shape(grid));
    for j in range(np.shape(grid)[0]):
        for k in range(np.shape(grid)[1]):
            for i in range(len(triangles)):
                tripath = matplotlib.path.Path(triangles[i])
                if tripath.contains_point((lons[k], lats[j])):
                    val_arr[j][k] = values[i];
                    break;
    return val_arr


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
    write_tape_eigenvectors(x, y, e1, e2, v00, v01, v10, v11)
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


# not working yet: decimates incorrectly!
def write_tape_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11):
    positive_file = open("Results/Results_Tape/"+"positive_eigs.txt", 'w');
    negative_file = open("Results/Results_Tape/"+"negative_eigs.txt", 'w');
    coord_box = [min(xdata), max(xdata), min(ydata), max(ydata)];
    inc = [0.02, 0.02]

    w1 = nn_interp(xdata, ydata, w1, coord_box, inc)[2]
    w2 = nn_interp(xdata, ydata, w2, coord_box, inc)[2]
    v00 = nn_interp(xdata, ydata, v00, coord_box, inc)[2]
    v01 = nn_interp(xdata, ydata, v01, coord_box, inc)[2]
    v10 = nn_interp(xdata, ydata, v10, coord_box, inc)[2]
    xdata, ydata, v11 = nn_interp(xdata, ydata, v11, coord_box, inc);

    eigs_dec = 12;

    do_not_print_value = 200;
    overmax_scale = 200;

    for j in range(len(ydata)):
        for k in range(len(xdata)):
            if np.mod(j, eigs_dec) == 0 and np.mod(k, eigs_dec) == 0:
                if w1[j][k] > 0:
                    scale = w1[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale,
                                                                 v10[j][k]*scale) );
                    positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale,
                                                                 -v10[j][k]*scale) );
                if w1[j][k] < 0:
                    scale = w1[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale,
                                                                 v10[j][k]*scale) );
                    negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale,
                                                                 -v10[j][k]*scale) );
                if w2[j][k] > 0:
                    scale = w2[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale,
                                                                 v11[j][k]*scale) );
                    positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale,
                                                                 -v11[j][k]*scale) );
                if w2[j][k] < 0:
                    scale = w2[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale,
                                                                 v11[j][k]*scale) );
                    negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale,
                                                                 -v11[j][k]*scale) );
    positive_file.close();
    negative_file.close();

    return;


# outputs a netcdf to the desired results directory
def output_tape(lon, lat, vals, outdir, file):
    netcdf_read_write.produce_output_netcdf(lon, lat, vals, 'per yr', outdir+file);
    print("Success fitting wavelet-generated data to the required grid!");
    return;
