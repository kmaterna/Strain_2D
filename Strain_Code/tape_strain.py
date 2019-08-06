import numpy as np
import strain_tensor_toolbox
import matplotlib.pyplot as plt
import netCDF4
import netcdf_functions
import scipy.interpolate as interp
import matplotlib.path

# This code was created to work with matlab scripts published on Github by Carl Tape under the name surfacevel2strain.
# The functions reformulate .txt files to be read by Tape's code.

# This function inputs a NAM velo file and selects the columns needed for the tape strain method, in correct order.
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



# infile = input_to_tape("../Example_data/NAM08_pbovelfile_feb2018.vel")
# output_to_tape(infile, "../../compearth/surfacevel2strain/data/", "NAM08.txt")
