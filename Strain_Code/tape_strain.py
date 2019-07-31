import numpy as np
import strain_tensor_toolbox
import matplotlib.pyplot as plt
import netCDF4
import netcdf_functions

# This code was created to work with matlab scripts published on Github by Carl Tape under the name surfacevel2strain.
# The first 2 functions reformulate txt files to be read by Tape's code.
# The other functions take the outputted files from the Tape code and compute strain properties and prepare them for mapping.

# This functioninputs a NAM velo file and selects the columns needed for the tape strain method, in correct order.
# These columns are: lon, lat, ve, vn, vu, se. sn, su, ren, reu, rnu, start, finish, name
def input_to_tape(filename):
	infile = np.loadtxt(filename, skiprows = 37, usecols = (8, 7, 20, 19, 21, 23, 22, 24, 25, 27, 26, 28, 29), unpack = True)
	
	lon = infile[0]
	for i in range(len(lon)):
		if lon[i] >=180:
			lon[i] = lon[i] - 360
	infile = np.vstack((lon, infile[1:]))

	infile[2:10] = infile[2:10]*1e3

	temp_file = open(filename, 'r')
	data = temp_file.readlines()[37:]
	temp_file.close()

	name_field = []
	for line in data:
	    name_field.append(line.split()[0])

	outfile = np.vstack((infile, name_field))
	outfile = np.transpose(outfile)

	return outfile


# Takes reformatted data and outputs it as a .txt for use in matlab scripts.
# Outdir should refer to location accessed by matlab scripts, and outdir and outfile should be strings.
def output_to_tape(data, outdir, outfile):
	print("Creating %s in %s" % (outfile, outdir))
	np.savetxt(outdir+outfile, data, delimiter=" ", fmt="%s")
	return


# This code works with .dat files outputted by Carl Tape's matlab code, surfacevel2strain
# First, run tape code, selecting to output gmt files.
# this code inputs the "strain" .dat file and the "D tensor 6 entries" .dat
# and outputs the coordinates and strain tensor components.
def input_from_tape(coordsfile, datafile):
	incoords = np.loadtxt(coordsfile, usecols = (0, 1), unpack = True)
	lon = incoords[0]
	lat = incoords[1]
	infile = np.loadtxt(datafile, skiprows = 1, usecols = (3, 4, 5), unpack = True)
	thth = infile[0]
	thph = infile[1]
	phph = infile[2]
	return lon, lat, thth, thph, phph

# This function computes second invariant, max shear strain, etc from symmetric strain tensor components
def compute(thth, thph, phph):
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
		max_shear.append(strain_tensor_toolbox.max_shear_strain(exx, exy, eyy))
		dilatation.append(e11+e22)
		I2nd.append(np.log10(np.abs(strain_tensor_toolbox.second_invariant(exx, exy, eyy))))
		# I2nd.append(10*np.abs(strain_tensor_toolbox.second_invariant(exx, exy, eyy)))
	return [I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation]

# This function configures longitude, latitude, and whichever strain values you'd like to plot 
# and returns them in a format that can be written to netcdf.
def configure(x, y, vals):
	ymult = np.count_nonzero(y == np.min(y)) # number of loops of y values
	print(ymult)
	ylen = int(len(y)/ymult) # length of each loop of y values
	newx = []
	newvals = []
	newy = []

	for i in range(0, len(x), ylen):
		elem = x[i]
		newx.append(elem)

	for i in range(0, ylen):
		newy.append(y[i])

	for i in np.arange(0, len(vals), len(newy)):
		newvals.append(vals[i:i+len(newy)])
	newvals = np.transpose(newvals)

	return newx, newy, newvals

# outputs a netcdf to the desired results directory to be used in GMT
def output_nc(lon, lat, vals, outdir, file):
	netcdf_functions.produce_output_netcdf(lon, lat, vals, 'per yr', outdir+file+".nc");
	return

# infile = input_to_tape("../Example_data/NAM08_pbovelfile_feb2018.vel")
# output_to_tape(infile, "../../compearth/surfacevel2strain/data/", "NAM08.txt")

myx, myy, mytt, mytp, mypp = input_from_tape("../../compearth/surfacevel2strain/matlab_output/cascadia_d02_q04_q05_b1_2D_s1_u1_strain.dat", "../../compearth/surfacevel2strain/matlab_output/cascadia_d02_q04_q05_b1_2D_s1_u1_Dtensor_6entries.dat")
[I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = compute(mytt, mytp, mypp)
newx, newy, newI2nd = configure(myx, myy, I2nd)
output_nc(newx, newy, newI2nd,"../Results/Results_Tape/", "I2nd")

# myx, myy, mytt, mytp, mypp = input_from_tape("../../compearth/surfacevel2strain/matlab_output/socal_d01_q03_q05_b1_2D_s1_u1_strain.dat", "../../compearth/surfacevel2strain/matlab_output/socal_d01_q03_q05_b1_2D_s1_u1_Dtensor_6entries.dat")
# [I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = compute(mytt, mytp, mypp)
# newx, newy, newI2nd = configure(myx, myy, I2nd)
# output_nc(newx, newy, newI2nd, "../Results/Results_Tape/", "I2nd_socal")

# print(newx)
# print(newy)
# print(max(max_shear))
# print(min(max_shear))
print("second invariant max:")
print(max(I2nd))
print("second invariant min:")
print(min(I2nd))
# print(np.mean(I2nd))
# print(max(dilatation))
# print(min(dilatation))

# print(myx[:200])
# print(myy[:200])
# print(newx)
# print(newy)
# print(len(newx))
# print(len(newy))
# print(newI2nd.shape)
# print(newI2nd)

