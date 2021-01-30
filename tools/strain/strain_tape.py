import numpy as np
import gps_io_functions
import datetime as dt

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

# if reading in UNR: use gps_input_functions

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
