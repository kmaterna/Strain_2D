import numpy as np 
import matplotlib.path
import matplotlib.pyplot as plt
import netcdf_functions
import strain_tensor_toolbox
import netCDF4
import scipy.interpolate as interp

# this function works with .txt files outputtd by the delaunay or hammond triangulation methods
def input_delaunay(outdir, filename):
	triangle_file = []
	index_file = []
	with open(outdir+filename, 'r') as file:
		for line in file:
			if line[3] == 'Z':
				index_file.append(float(line[4:-2]))
			else:
				string = line.split()
				triangle_file.append(float(string[0]))
				triangle_file.append(float(string[1]))
	return triangle_file, index_file

# This code works with .dat files outputted by Carl Tape's matlab code, surfacevel2strain
# First, run tape code, selecting to output gmt files.
# this code inputs the "strain" .dat file and the "D tensor 6 entries" .dat
# and outputs the coordinates and strain tensor components.
def input_tape(indir, coordsfile, datafile):
	incoords = np.loadtxt(indir+coordsfile, usecols = (0, 1), unpack = True)
	lon = incoords[0]
	lat = incoords[1]
	infile = np.loadtxt(indir+datafile, skiprows = 1, usecols = (3, 4, 5), unpack = True)
	thth = infile[0]
	thph = infile[1]
	phph = infile[2]
	return lon, lat, thth, thph, phph

# makes grid for delaunay
def make_grid(lonmin, lonmax, latmin, latmax, stepsize):
	lons = np.arange(lonmin, lonmax, stepsize)
	lats = np.arange(latmin, latmax, stepsize)
	grid = []
	for i in range(len(lats)):
		for j in range(len(lons)):
			grid.append([lons[j], lats[i]])
	return grid, lons, lats

# gathers delaunay-triangulated points into coordinate triples representing the triangle vertices
def make_triangles(triangles):
	new_triangles = []
	for i in range(0, len(triangles), 6):
		triangle = np.empty([3, 2])
		triangle[0][0] = triangles[i]
		triangle[0][1] = triangles[i+1] 
		triangle[1][0] = triangles[i+2]
		triangle[1][1] = triangles[i+3]
		triangle[2][0] = triangles[i+4]
		triangle[2][1] = triangles[i+5]
		new_triangles.append(triangle)
	return new_triangles

# searches path created by triangle vertices for each gridpoint, then assigns that triangle's value to the gridpoint
def find_in_triangles(triangles, index_file, grid):
	val_arr = np.nan*np.ones(len(grid))
	for j in range(len(grid)):
		for i in range(len(triangles)):
			tripath = matplotlib.path.Path(triangles[i])
			if tripath.contains_point(grid[j]):
				val_arr[j] = index_file[i]
				break
	return val_arr

# arranges values for outputting a netcdf file for gmt plotting
def configure_vals(val_arr, lons, lats):
	new_vals = []
	for i in range(0, len(val_arr), len(lons)):
		new_vals.append(np.array(val_arr[i:i+len(lons)]))
	return new_vals

# This function computes second invariant, max shear strain, etc from symmetric strain tensor components (in spherical coords, from tape)
def compute_tape(thth, thph, phph):
	I2nd=[];
	rot=[];
	max_shear=[];
	e1=[]; # eigenvalues
	e2=[];
	v00=[];  # eigenvectors
	v01=[];
	v10=[];
	v11=[];
	dilatation=[]; # dilatation	= e1+e2

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
def nn_interp(x, y, vals, xmin, xmax, ymin, ymax, inc):

	newx = np.arange(xmin, xmax, inc)
	newy = np.arange(ymin, ymax, inc)

	tempvals = []
	newgrid = []

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
	positive_file=open("Results/Results_Tape/"+"positive_eigs.txt",'w');
	negative_file=open("Results/Results_Tape/"+"negative_eigs.txt",'w');
	

	w1 = nn_interp(xdata, ydata, w1, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)[2]
	w2 = nn_interp(xdata, ydata, w2, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)[2]
	v00 = nn_interp(xdata, ydata, v00, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)[2]
	v01 = nn_interp(xdata, ydata, v01, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)[2]
	v10 = nn_interp(xdata, ydata, v10, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)[2]
	xdata, ydata, v11 = nn_interp(xdata, ydata, v11, min(xdata), max(xdata), min(ydata), max(ydata), 0.02)

		
	eigs_dec=12;

	do_not_print_value=200;
	overmax_scale=200;

	for j in range(len(ydata)):
		for k in range(len(xdata)):
			if np.mod(j,eigs_dec)==0 and np.mod(k,eigs_dec)==0:
				if w1[j][k]>0:
					scale=w1[j][k];
					if abs(scale)>do_not_print_value:
						scale=overmax_scale;
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale, v10[j][k]*scale) );
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale, -v10[j][k]*scale) );
				if w1[j][k]<0:
					scale=w1[j][k];
					if abs(scale)>do_not_print_value:
						scale=overmax_scale;					
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k]*scale, v10[j][k]*scale) );
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k]*scale, -v10[j][k]*scale) );
				if w2[j][k]>0:
					scale=w2[j][k];
					if abs(scale)>do_not_print_value:
						scale=overmax_scale;
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale, v11[j][k]*scale) );
					positive_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale, -v11[j][k]*scale) );
				if w2[j][k]<0:
					scale=w2[j][k];
					if abs(scale)>do_not_print_value:
						scale=overmax_scale;
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k]*scale, v11[j][k]*scale) );
					negative_file.write("%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k]*scale, -v11[j][k]*scale) );
	positive_file.close();
	negative_file.close();


	return;



# outputs a netcdf to the desired results directory to be used in GMT
def output_tape(lon, lat, vals, outdir, file):
	netcdf_functions.produce_output_netcdf(lon, lat, vals, 'per yr', outdir+file);
	print("Success fitting wavelet-generated data to the required grid!");
	return;

# outputs a netcdf to the desired results directory to be used in GMT
def output_delaunay(xdata, ydata, vals, outdir, indexfile, outfile):
	valfile = open(outdir+indexfile, 'w');
	for val in vals:
		valfile.write(str(val)+'\n');
	valfile.close();
	netcdf_functions.produce_output_netcdf(xdata, ydata, vals, 'per year', outdir+outfile);
	print("Success fitting triangulated data to the required grid!");
	return;