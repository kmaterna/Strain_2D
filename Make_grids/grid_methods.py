import produce_gridded as pg
import configure_functions as config




def drive(method):

	if method == "delaunay":
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		trifile, indexfile = pg.input_delaunay(myParams.outdir, "I2nd.txt");
		grid, lons, lats = pg.make_grid(coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		triangles = pg.make_triangles(trifile);
		gridvals = pg.find_in_triangles(triangles, indexfile, grid);
		newvals = pg.configure_vals(gridvals, lons, lats);
		pg.output_delaunay(lons, lats, newvals, outdir, "I2nd_index_file.txt","I2nd.nc");

	elif method == "tape":
		[myParams] = config.configure(method);
		map_range = [float(i) for i in myParams.map_range.split("/")]
		[grid_inc, coord_box, coord_box_data, outdir, gmtfile] = config.get_tunable_options(method, map_range);
		indir = "../compearth/surfacevel2strain/matlab_output/"
		x, y, tt, tp, pp = pg.input_tape(indir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat", "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
		[I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = pg.compute_tape(tt, tp, pp);
		newx, newy, newI2nd = pg.nn_interp(x, y, I2nd, coord_box[0], coord_box[1], coord_box[2], coord_box[3], grid_inc);
		pg.output_tape(newx, newy, newI2nd, outdir, "I2nd.nc");

	else:
		print("Oops! Method not recognized, try again.")
	return;