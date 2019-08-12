import produce_gridded as pg
import configure_functions as config
import strain_tensor_toolbox
import netcdf_functions



def drive(method):

	if method == "hammond":

		# second invariant
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "I2nd.txt");
		grid, lons, lats = pg.make_grid(coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, outdir, "I2nd_index_file.txt","I2nd.nc");

		# azimuth
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "azimuth.txt");
		grid, lons, lats = pg.make_grid(coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, outdir, "azimuth_index_file.txt","azimuth.nc");

		# dilatation
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "Dilatation.txt");
		grid, lons, lats = pg.make_grid(coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, outdir, "dilatation_index_file.txt","dilatation.nc");

		# max shear
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "max_shear.txt");
		grid, lons, lats = pg.make_grid(coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, outdir, "max_shear_index_file.txt","max_shear.nc");


	elif method == "tape":
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		indir = "../compearth/surfacevel2strain/matlab_output/"
		x, y, tt, tp, pp = pg.input_tape(indir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat", "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
		[I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = pg.compute_tape(tt, tp, pp);
		# second invariant
		newx, newy, newI2nd = pg.nn_interp(x, y, I2nd, coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		pg.output_tape(newx, newy, newI2nd, outdir, "I2nd.nc");
		# dilatation
		newx, newy, newdila = pg.nn_interp(x, y, dilatation, coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		pg.output_tape(newx, newy, newdila, outdir, "dilatation.nc");
		# max shear
		newx, newy, newmax = pg.nn_interp(x, y, max_shear, coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		pg.output_tape(newx, newy, newmax, outdir, "max_shear.nc");		
		# azimuth
		azimuth = strain_tensor_toolbox.max_shortening_azimuth_1d(e1, e2, v00, v01, v10, v11)
		newx, newy, newaz = pg.nn_interp(x, y, azimuth, coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		netcdf_functions.produce_output_netcdf(newx, newy, newaz, 'degrees', outdir+'azimuth.nc');


	else:
		print("Oops! Method not recognized, try again.")
	return;