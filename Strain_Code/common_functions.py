# The configure, input, and output steps for GPS Strain analysis. 
# map_range=[-125, -121, 37.0, 42.2]; # Northern California
# map_range=[-121, -115, 32, 36]; # Southern California
# map_range=[-123, -118, 34, 38]; # Central California
# map_range=[-125, -114, 32.0, 42.0]; # ALL California	
# map_range=[-125, -110, 33.0, 48.2]; # Western US


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



# ----------------- INPUTS -------------------------
def inputs(MyParams):
	# Purpose: generate input velocity field. 
	print("Reading %s" % MyParams.input_file);
	if 'PBO' in MyParams.input_file or 'pbo' in MyParams.input_file:
		[myVelfield]=gps_input_functions.read_pbo_vel_file(MyParams.input_file);  # read the raw velfield from file. 
	elif 'MAGNET' in MyParams.input_file or 'unr' in MyParams.input_file:
		[myVelfield]=gps_input_functions.read_unr_vel_file(MyParams.input_file);  # read the raw velfield from file. 
	else:
		print("Error! Cannot read %s " % MyParams.input_file);
		sys.exit(1);
	print("%d stations before applying coord_box." % (len(myVelfield.name)) );
	[myVelfield]=gps_input_functions.clean_velfield(myVelfield, MyParams.num_years, MyParams.max_sigma, MyParams.coord_box_data);
	[myVelfield]=gps_input_functions.remove_duplicates(myVelfield);
	print("%d stations after selection criteria.\n" % (len(myVelfield.name)) );
	return [myVelfield];




# ----------------- OUTPUTS -------------------------

def outputs_2d(xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation, myVelfield, MyParams):
	print("Writing 2d outputs:");
	outfile=open(MyParams.outdir+"tempgps.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();
	netcdf_functions.produce_output_netcdf(xdata, ydata, I2nd, 'per yr', MyParams.outdir+'I2nd.nc');
	netcdf_functions.flip_if_necessary(MyParams.outdir+'I2nd.nc');
	netcdf_functions.produce_output_netcdf(xdata, ydata, rot, 'per yr', MyParams.outdir+'rot.nc');
	netcdf_functions.flip_if_necessary(MyParams.outdir+'rot.nc');
	netcdf_functions.produce_output_netcdf(xdata, ydata, dilatation, 'per yr', MyParams.outdir+'dila.nc');
	netcdf_functions.flip_if_necessary(MyParams.outdir+'dila.nc');
	print("Max I2: %f " % (np.amax(I2nd)));
	print("Max rot: %f " % (np.amax(rot)));
	print("Min rot: %f " % (np.amin(rot)));
	write_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11, MyParams);
	print("../../"+MyParams.gmtfile+" "+MyParams.map_range);
	subprocess.call("../../"+MyParams.gmtfile+" "+MyParams.map_range,shell=True,cwd=MyParams.outdir);
	return;

def write_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11, MyParams):
	# Need eigs_interval and outdir from MyParams. 
	positive_file=open(MyParams.outdir+"positive_eigs.txt",'w');
	negative_file=open(MyParams.outdir+"negative_eigs.txt",'w');
	if MyParams.strain_method=='visr':
		eigs_dec=8;
	elif MyParams.strain_method=='gpsgridder':
		eigs_dec=12;
	elif MyParams.strain_method=='spline':
		eigs_dec=12;
	else:
		print("Error! strain method not recognized for eigenvector plotting.");

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






def outputs_1d(xcentroid, ycentroid, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation, myVelfield, MyParams):
	print("Writing 1d outputs:");

	rotfile=open(MyParams.outdir+"rotation.txt",'w');
	I2ndfile=open(MyParams.outdir+"I2nd.txt",'w');
	Dfile=open(MyParams.outdir+"Dilatation.txt",'w');
	positive_file=open(MyParams.outdir+"positive_eigs.txt",'w');
	negative_file=open(MyParams.outdir+"negative_eigs.txt",'w');

	outfile=open(MyParams.outdir+"tempgps.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();

	for i in range(len(I2nd)):
		# Write the triangle's rotation
		rotfile.write("> -Z"+str(rot[i])+"\n");
		rotfile.write(str(polygon_vertices[i,0,0])+" "+str(polygon_vertices[i,0,1])+"\n");
		rotfile.write(str(polygon_vertices[i,1,0])+" "+str(polygon_vertices[i,1,1])+"\n");
		rotfile.write(str(polygon_vertices[i,2,0])+" "+str(polygon_vertices[i,2,1])+"\n");

		# Write the triangle's I2
		I2ndfile.write("> -Z"+str(I2nd[i])+"\n");
		I2ndfile.write(str(polygon_vertices[i,0,0])+" "+str(polygon_vertices[i,0,1])+"\n");
		I2ndfile.write(str(polygon_vertices[i,1,0])+" "+str(polygon_vertices[i,1,1])+"\n");
		I2ndfile.write(str(polygon_vertices[i,2,0])+" "+str(polygon_vertices[i,2,1])+"\n");
		
		# Write the dilatation
		Dfile.write("> -Z"+str(dilatation[i])+"\n"); 
		Dfile.write(str(polygon_vertices[i,0,0])+" "+str(polygon_vertices[i,0,1])+"\n");
		Dfile.write(str(polygon_vertices[i,1,0])+" "+str(polygon_vertices[i,1,1])+"\n");
		Dfile.write(str(polygon_vertices[i,2,0])+" "+str(polygon_vertices[i,2,1])+"\n");

		# Write the eigenvectors and eigenvalues
		write_single_eigenvector(positive_file, negative_file, e1[i], v00[i], v10[i], xcentroid[i], ycentroid[i]);
		write_single_eigenvector(positive_file, negative_file, e2[i], v01[i], v11[i], xcentroid[i], ycentroid[i]);
	
	print("Max I2: %f " % (max(I2nd)));
	print("Max rot: %f " % (max(rot)));
	print("Min rot: %f " % (min(rot)));

	rotfile.close();
	I2ndfile.close();
	Dfile.close();
	positive_file.close();
	negative_file.close();
	print("../../"+MyParams.gmtfile+" "+MyParams.map_range);
	subprocess.call("../../"+MyParams.gmtfile+" "+MyParams.map_range,shell=True,cwd=MyParams.outdir);
	return;


def write_single_eigenvector(positive_file, negative_file, e, v0, v1, x, y):
	# e = eigenvalue, [v0, v1] = eigenvector. 
	# Writes a single eigenvector eigenvalue pair. 
	# Also has functionality to saturate eigenvectors so they don't blow up. 
	overall_max=40.0;
	scale=0.4*e;

	vx=v0*scale;
	vy=v1*scale;
	if np.sqrt(vx*vx+vy*vy)>overall_max:
		scale=scale*(overall_max/np.sqrt(vx*vx+vy*vy))
		vx=v0*scale;
		vy=v1*scale;

	if e>0:
		positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy ));
		positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy ));
	else:	
		negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy ));
		negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy ));
	return;



