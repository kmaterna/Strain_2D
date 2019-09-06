
import numpy as np 
import subprocess, sys, os
import collections
import gps_input_functions
import netcdf_functions


Params=collections.namedtuple("Params",['strain_method','input_file','map_range','coord_box','coord_box_data','num_years','max_sigma','grid_inc','outdir','gmtfile','blacklist_file']);

# Inputs that we've tried and tested: 
# Bartlow_ETSvels.txt
# Bartlow_interETSvels.txt
# NAM08_pbovelfile_feb2018.txt
# midas.NA12.txt


# ----------------- CONFIGURE -------------------------
def configure(strain_method):
	
	input_file="Vel_Data/midas.NA12.txt";
	# input_file="Vel_Data/Other_vels/Bartlow_ETSvels.txt";
	# input_file="Vel_Data/NAM08_pbovelfile_feb2018.vel";

	map_range=[-125, -121, 37.0, 42.2]; # Northern California
	map_range_string = str(map_range[0])+'/'+str(map_range[1])+'/'+str(map_range[2])+'/'+str(map_range[3]);
	num_years=3.0;
	max_sigma=2.0;
	[grid_inc, coord_box, coord_box_data, gmtfile] = get_tunable_options(strain_method, map_range);
	blacklist_file = "Vel_Data/blacklist_stations.txt";
	outdir = get_output_dir(input_file, strain_method);

	print("\n------------------------------");
	print("Hello! We are...");
	print("   Computing strain using : %s " % strain_method);
	print("   Input data from        : %s" % input_file );
	print("   Coordinates            : [%.2f, %.2f, %.2f, %.2f]" % (coord_box[0], coord_box[1], coord_box[2], coord_box[3]) );
	print("   Putting the outputs    : %s \n" % outdir);
	MyParams=Params(strain_method=strain_method, input_file=input_file, map_range=map_range_string, coord_box=coord_box, 
		coord_box_data=coord_box_data, num_years=num_years, max_sigma=max_sigma, grid_inc=grid_inc, outdir=outdir, 
		gmtfile=gmtfile, blacklist_file=blacklist_file);
	print_configure_options(MyParams); # record-keeping
	return [MyParams];


def get_output_dir(input_file, strain_method):
	# Making a nested output directory structure for each different input file. 
	outer_name=input_file.split('/')[-1].split('.')[0];  # example: for Vel_Data/unr_velo_clean.txt, we have "unr_velo_clean"
	inner_name="Results_"+strain_method;
	outdir_outer="Results/"+outer_name+"/";
	outdir_inner=outdir_outer+inner_name+"/";
	subprocess.call(['mkdir','-p',outdir_outer],shell=False); 
	subprocess.call(['mkdir','-p',outdir_inner],shell=False); 
	return outdir_inner; 


# A couple of options that change based on strain method. 
def get_tunable_options(strain_method, map_range):
	if strain_method=="gpsgridder":
		grid_inc   =0.04;
		coord_box  =[map_range[0]-1, map_range[1]+3, map_range[2]-2, map_range[3]+2];
		coord_box_data = coord_box;
		gmtfile    ="Strain_Code/GMT_mapping_codes/gpsgridder_gmt.gmt";

	elif strain_method=="visr":
		grid_inc   =0.04;
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		gmtfile    ="Strain_Code/GMT_mapping_codes/visr_gmt.gmt";

	elif strain_method=="delaunay":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		gmtfile    ="Strain_Code/GMT_mapping_codes/delaunay_gmt.gmt"

	elif strain_method=="hammond":
		grid_inc   =0.04; # larger interval for convenience, because it's slow. 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];
		coord_box_data = coord_box;
		gmtfile    ="Strain_Code/GMT_mapping_codes/hammond_gmt.gmt"

	elif strain_method=="spline":
		grid_inc   =0.04; 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];  # big box, small box
		coord_box_data = [-126, -114, 32, 48];
		gmtfile    ="Strain_Code/GMT_mapping_codes/spline_gmt.gmt"

	elif strain_method=="ND_interp":
		grid_inc   =0.04; 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];  # big box, small box
		coord_box_data = [-126, -114, 32, 48,];
		gmtfile    ="Strain_Code/GMT_mapping_codes/ND_interp_gmt.gmt"

	elif strain_method=="tape":
		grid_inc   =0.04; 
		coord_box  =[map_range[0]-0.5, map_range[1]+0.5, map_range[2]-0.5, map_range[3]+0.5];  # big box, small box
		coord_box_data = coord_box;
		gmtfile    ="Strain_Code/GMT_mapping_codes/tape_gmt.gmt"

	else:
		print("ERROR: "+strain_method+" is not a known strain method. ");
		sys.exit(1);

	return [grid_inc, coord_box, coord_box_data, gmtfile];

def print_configure_options(MyParams):
	ofile=open(MyParams.outdir+"config_options.txt",'w');
	ofile.write("Strain Method      : %s\n" % MyParams.strain_method);
	ofile.write("Input file         : %s\n" % MyParams.input_file);
	ofile.write("Calc. Coordinates  : [%.2f, %.2f, %.2f, %.2f]\n" % (MyParams.coord_box[0], MyParams.coord_box[1], MyParams.coord_box[2], MyParams.coord_box[3]) );	
	ofile.write("Data Range         : [%.2f, %.2f, %.2f, %.2f]\n" % (MyParams.coord_box_data[0], MyParams.coord_box_data[1], MyParams.coord_box_data[2], MyParams.coord_box_data[3]) );	
	ofile.write("grid_inc           : %.2f\n" % MyParams.grid_inc);
	ofile.write("Num_years          : %.2f\n" % MyParams.num_years);
	ofile.write("Max Sigma          : %.2f\n" % MyParams.max_sigma);
	ofile.write("Blacklist file     : %s\n" % MyParams.blacklist_file);
	ofile.write("Outdir             : %s\n" % MyParams.outdir);
	ofile.close();
	return;
