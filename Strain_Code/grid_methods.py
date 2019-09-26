# Methods for turning a 1D computation into a 2D computation

import produce_gridded as pg
import configure_functions as config
import strain_tensor_toolbox
import netcdf_functions
import numpy as np



def drive(method):

	if method == "hammond":

		# generic steps
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		grid, lons, lats = pg.make_grid(myParams.coord_box[0], myParams.coord_box[1], myParams.coord_box[2], myParams.coord_box[3], myParams.grid_inc);
		print("Producing gridded dataset of: ")

		# second invariant
		print("... 2nd Invariant...");
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "I2nd.txt");
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, myParams.outdir, "I2nd_index_file.txt","I2nd.nc");

		# dilatation
		print("... Dilatation...");
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "dilatation.txt");
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, myParams.outdir, "dilatation_index_file.txt","dila.nc");

		# max shear
		print("... Max Shear...");
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "max_shear.txt");
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, myParams.outdir, "max_shear_index_file.txt","max_shear.nc");

		# azimuth
		print("... Azimuth...");
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "azimuth.txt");
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, myParams.outdir, "azimuth_index_file.txt","azimuth.nc");


	elif method == "tape":
		# generic steps
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		coord_box=myParams.coord_box;
		indir = "../compearth/surfacevel2strain/matlab_output/"
		print("Producing gridded dataset of: ")
		x, y, tt, tp, pp = pg.input_tape(indir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat", "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
		[I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = pg.compute_tape(tt, tp, pp);
		# second invariant
		newx, newy, newI2nd = pg.nn_interp(x, y, I2nd, coord_box[0], coord_box[1], coord_box[2], coord_box[3], myParams.grid_inc);
		pg.output_tape(newx, newy, newI2nd, myParams.outdir, "I2nd.nc");
		# dilatation
		newx, newy, newdila = pg.nn_interp(x, y, dilatation, coord_box[0], coord_box[1], coord_box[2], coord_box[3], myParams.grid_inc);
		pg.output_tape(newx, newy, newdila, myParams.outdir, "dila.nc");
		# max shear
		newx, newy, newmax = pg.nn_interp(x, y, max_shear, coord_box[0], coord_box[1], coord_box[2], coord_box[3], myParams.grid_inc);
		pg.output_tape(newx, newy, newmax, myParams.outdir, "max_shear.nc");		
		# azimuth
		azimuth = strain_tensor_toolbox.max_shortening_azimuth_1d(e1, e2, v00, v01, v10, v11)
		newx, newy, newaz = pg.nn_interp(x, y, azimuth, coord_box[0], coord_box[1], coord_box[2], coord_box[3], myParams.grid_inc);
		netcdf_functions.produce_output_netcdf(newx, newy, newaz, 'degrees', myParams.outdir+'azimuth.nc');
		pg.write_tape_eigenvectors(x, y, e1, e2, v00, v01, v10, v11)

	else:
		print("Oops! Method not recognized, try again.")
	return;