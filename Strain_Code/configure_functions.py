
import numpy as np 
import subprocess, sys, os
import collections
import gps_input_functions
import netcdf_functions


Params=collections.namedtuple("Params",['strain_method','input_file','map_range','coord_box','coord_box_data','num_years','max_sigma','grid_inc','outdir','gmtfile']);


# ----------------- CONFIGURE -------------------------
def configure(strain_method):
	input_file="Example_data/NAM08_pbovelfile_feb2018.vel";
	map_range=[-125, -121, 37.0, 42.2]; # Northern California
	map_range_string = str(map_range[0])+'/'+str(map_range[1])+'/'+str(map_range[2])+'/'+str(map_range[3]);
	num_years=3.0;
	max_sigma=2.0;
	[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = get_tunable_options(strain_method, map_range);

	print("\n------------------------------");
	print("Hello! We are...");
	print("   Computing strain using : %s " % strain_method);
	print("   Input data from        : %s" % input_file );
	print("   Coordinates            : [%.2f, %.2f, %.2f, %.2f]" % (coord_box[0], coord_box[1], coord_box[2], coord_box[3]) );
	print("   Putting the outputs    : %s \n" % outdir);
	subprocess.call(['mkdir','-p',outdir],shell=False);
	MyParams=Params(strain_method=strain_method, input_file=input_file, map_range=map_range_string, coord_box=coord_box, coord_box_data=coord_box_data, num_years=num_years, max_sigma=max_sigma, grid_inc=grid_inc, outdir=outdir, gmtfile=gmtfile);
	return [MyParams];


# A couple of options that change based on strain method. 
def get_tunable_options(strain_method, map_range):
	if strain_method=="gpsgridder":
		grid_inc   =0.02;
		coord_box  =[map_range[0]-1, map_range[1]+3, map_range[2]-2, map_range[3]+2];
		coord_box_data = coord_box;
		outdir     ="Results/Results_GPSgridder/";
		gmtfile    ="Strain_Code/GMT_mapping_codes/gpsgridder_gmt.gmt";

	elif strain_method=="visr":
		grid_inc   =0.04;
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		outdir     ="Results/Results_Visr/";
		gmtfile    ="Strain_Code/GMT_mapping_codes/visr_gmt.gmt";

	elif strain_method=="delaunay":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		outdir     ="Results/Results_Delaunay/"
		gmtfile    ="Strain_Code/GMT_mapping_codes/delaunay_gmt.gmt"

	elif strain_method=="hammond":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		outdir     ="Results/Results_Hammond/"
		gmtfile    ="Strain_Code/GMT_mapping_codes/hammond_gmt.gmt"

	elif strain_method=="spline":
		grid_inc   =0.02; 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];  # big box, small box
		coord_box_data = [-126, -114, 32, 48,];
		outdir     ="Results/Results_Numpy_Spline/"
		gmtfile    ="Strain_Code/GMT_mapping_codes/spline_gmt.gmt"

	else:
		print("ERROR: "+strain_method+" is not a known strain method. ");
		sys.exit(1);

	return [grid_inc, coord_box, coord_box_data, outdir, gmtfile];
