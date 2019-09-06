# The input manager for GPS Strain analysis. 

import numpy as np 
import subprocess, sys, os
import collections
import gps_input_functions

Params=collections.namedtuple("Params",['strain_method','input_file','map_range','coord_box','coord_box_data','num_years','max_sigma','grid_inc','outdir','gmtfile']);

# ----------------- INPUTS -------------------------
def inputs(MyParams):
	# Purpose: generate input velocity field. 
	print("Reading %s" % MyParams.input_file);
	if 'PBO' in MyParams.input_file or 'pbo' in MyParams.input_file:
		[myVelfield]=gps_input_functions.read_pbo_vel_file(MyParams.input_file);  # read the raw velfield from file. 
	elif 'MAGNET' in MyParams.input_file or 'unr' in MyParams.input_file or 'midas' in MyParams.input_file:
		[myVelfield]=gps_input_functions.read_unr_vel_file(MyParams.input_file);  # read the raw velfield from file. 
	elif 'ETS' in MyParams.input_file:
		[myVelfield]=gps_input_functions.read_ETS_vel_file(MyParams.input_file);
	else:
		print("Error! Cannot read %s " % MyParams.input_file);
		sys.exit(1);
	print("%d stations before applying coord_box." % (len(myVelfield.name)) );
	[myVelfield]=gps_input_functions.remove_blacklist(myVelfield, MyParams.blacklist_file);
	[myVelfield]=gps_input_functions.clean_velfield(myVelfield, MyParams.num_years, MyParams.max_sigma, MyParams.coord_box_data);
	[myVelfield]=gps_input_functions.remove_duplicates(myVelfield);
	print("%d stations after selection criteria.\n" % (len(myVelfield.name)) );
	return [myVelfield];

