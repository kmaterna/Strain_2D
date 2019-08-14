import compare_methods as comp 



def drive(component):

	if component == "azimuth":
		print("Comparing %s across all methods" % component)
		# directories to netcdfs for each method
		file1 = "Results/Results_Hammond/azimuth.nc"
		file2 = "Results/Results_Numpy_Spline/azimuth.nc"
		file3 = "Results/Results_Visr/azimuth.nc"
		file4 = "Results/Results_ND_interp/azimuth.nc"
		file5 = "Results/Results_Tape/azimuth.nc"

		lon1, lat1, val1 = comp.input_netcdf(file1)
		lon2, lat2, val2 = comp.input_netcdf(file2)
		lon3, lat3, val3 = comp.input_netcdf(file3)
		lon4, lat4, val4 = comp.input_netcdf(file4)
		lon5, lat5, val5 = comp.input_netcdf(file5)

		print("delaunay range: %.2f %.2f %.2f %.2f " % ( min(lon1), max(lon1), min(lat1), max(lat1) ) );
		print("spline range: %.2f %.2f %.2f %.2f " % ( min(lon2), max(lon2), min(lat2), max(lat2) ) );
		print("visr range: %.2f %.2f %.2f %.2f " % ( min(lon3), max(lon3), min(lat3), max(lat3) ) );
		print("nd interp range: %.2f %.2f %.2f %.2f " % ( min(lon4), max(lon4), min(lat4), max(lat4) ) );
		print("tape range: %.2f %.2f %.2f %.2f " % ( min(lon5), max(lon5), min(lat5), max(lat5) ) );

		lons1, lats1, val1 = comp.confine_to_grid(lon1, lat1, val1, -124.38, -121.2, 37.2, 42, 0.04)
		lons2, lats2, val2 = comp.confine_to_grid(lon2, lat2, val2, -124.38, -121.2, 37.2, 42, 0.04)
		lons3, lats3, val3 = comp.confine_to_grid(lon3, lat3, val3, -124.38, -121.2, 37.2, 42, 0.04)
		lons4, lats4, val4 = comp.confine_to_grid(lon4, lat4, val4, -124.38, -121.2, 37.2, 42, 0.04)
		lons5, lats5, val5 = comp.confine_to_grid(lon5, lat5, val5, -124.38, -121.2, 37.2, 42, 0.04)

		comp.check_coregistration(val1, val2, val3, val4, val5);
		my_means = comp.angle_means(lons2, lats2, val1, val2, val3, val4, val5)
		my_sds = comp.angle_sds(lons2, lats2, val1, val2, val3, val4, val5)
		comp.output_nc(lons2, lats2, my_means, "means", "azimuth")
		comp.output_nc(lons2, lats2, my_sds, "deviations", 'azimuth')

	elif component == "I2nd":
		print("Comparing %s across all methods" % component)

		file1 = "Results/Results_Hammond/I2nd.nc"
		file2 = "Results/Results_Numpy_Spline/I2nd.nc"
		file3 = "Results/Results_Visr/I2nd.nc"
		file4 = "Results/Results_ND_interp/I2nd.nc"
		file5 = "Results/Results_Tape/I2nd.nc"

		lon1, lat1, val1 = comp.input_netcdf(file1)
		lon2, lat2, val2 = comp.input_netcdf(file2)
		lon3, lat3, val3 = comp.input_netcdf(file3)
		lon4, lat4, val4 = comp.input_netcdf(file4)
		lon5, lat5, val5 = comp.input_netcdf(file5)

		print("delaunay range: %.2f %.2f %.2f %.2f " % ( min(lon1), max(lon1), min(lat1), max(lat1) ) );
		print("spline range: %.2f %.2f %.2f %.2f " % ( min(lon2), max(lon2), min(lat2), max(lat2) ) );
		print("visr range: %.2f %.2f %.2f %.2f " % ( min(lon3), max(lon3), min(lat3), max(lat3) ) );
		print("nd interp range: %.2f %.2f %.2f %.2f " % ( min(lon4), max(lon4), min(lat4), max(lat4) ) );
		print("tape range: %.2f %.2f %.2f %.2f " % ( min(lon5), max(lon5), min(lat5), max(lat5) ) );

		lons1, lats1, val1 = comp.confine_to_grid(lon1, lat1, val1, -124.38, -121.2, 37.2, 42, 0.04)
		lons2, lats2, val2 = comp.confine_to_grid(lon2, lat2, val2, -124.38, -121.2, 37.2, 42, 0.04)
		lons3, lats3, val3 = comp.confine_to_grid(lon3, lat3, val3, -124.38, -121.2, 37.2, 42, 0.04)
		lons4, lats4, val4 = comp.confine_to_grid(lon4, lat4, val4, -124.38, -121.2, 37.2, 42, 0.04)
		lons5, lats5, val5 = comp.confine_to_grid(lon5, lat5, val5, -124.38, -121.2, 37.2, 42, 0.04)

		comp.check_coregistration(val1, val2, val3, val4, val5);
		my_means = comp.grid_means_log(lons2, lats2, val1, val2, val3, val4, val5)
		my_sds = comp.grid_sds_log(lons2, lats2, val1, val2, val3, val4, val5)

		comp.output_nc(lons2, lats2, my_means, "means", "I2nd")
		comp.output_nc(lons2, lats2, my_sds, "deviations", "I2nd")


	elif component == "dilatation":
		print("Comparing %s across all methods" % component)

		file1 = "Results/Results_Hammond/dila.nc"
		file2 = "Results/Results_Numpy_Spline/dila.nc"
		file3 = "Results/Results_Visr/dila.nc"
		file4 = "Results/Results_ND_interp/dila.nc"
		file5 = "Results/Results_Tape/dila.nc"

		lon1, lat1, val1 = comp.input_netcdf(file1)
		lon2, lat2, val2 = comp.input_netcdf(file2)
		lon3, lat3, val3 = comp.input_netcdf(file3)
		lon4, lat4, val4 = comp.input_netcdf(file4)
		lon5, lat5, val5 = comp.input_netcdf(file5)

		print("delaunay range: %.2f %.2f %.2f %.2f " % ( min(lon1), max(lon1), min(lat1), max(lat1) ) );
		print("spline range: %.2f %.2f %.2f %.2f " % ( min(lon2), max(lon2), min(lat2), max(lat2) ) );
		print("visr range: %.2f %.2f %.2f %.2f " % ( min(lon3), max(lon3), min(lat3), max(lat3) ) );
		print("nd interp range: %.2f %.2f %.2f %.2f " % ( min(lon4), max(lon4), min(lat4), max(lat4) ) );
		print("tape range: %.2f %.2f %.2f %.2f " % ( min(lon5), max(lon5), min(lat5), max(lat5) ) );

		lons1, lats1, val1 = comp.confine_to_grid(lon1, lat1, val1, -124.38, -121.2, 37.2, 42, 0.04)
		lons2, lats2, val2 = comp.confine_to_grid(lon2, lat2, val2, -124.38, -121.2, 37.2, 42, 0.04)
		lons3, lats3, val3 = comp.confine_to_grid(lon3, lat3, val3, -124.38, -121.2, 37.2, 42, 0.04)
		lons4, lats4, val4 = comp.confine_to_grid(lon4, lat4, val4, -124.38, -121.2, 37.2, 42, 0.04)
		lons5, lats5, val5 = comp.confine_to_grid(lon5, lat5, val5, -124.38, -121.2, 37.2, 42, 0.04)

		comp.check_coregistration(val1, val2, val3, val4, val5);
		my_means = comp.grid_means(lons2, lats2, val1, val2, val3, val4, val5)
		my_sds = comp.grid_sds(lons2, lats2, val1, val2, val3, val4, val5)

		comp.output_nc(lons2, lats2, my_means, "means", "dilatation")
		comp.output_nc(lons2, lats2, my_sds, "deviations", "dilatation")

	elif component == "max_shear":
		print("Comparing %s across all methods" % component)

		file1 = "Results/Results_Hammond/max_shear.nc"
		file2 = "Results/Results_Numpy_Spline/max_shear.nc"
		file3 = "Results/Results_Visr/max_shear.nc"
		file4 = "Results/Results_ND_interp/max_shear.nc"
		file5 = "Results/Results_Tape/max_shear.nc"

		lon1, lat1, val1 = comp.input_netcdf(file1)
		lon2, lat2, val2 = comp.input_netcdf(file2)
		lon3, lat3, val3 = comp.input_netcdf(file3)
		lon4, lat4, val4 = comp.input_netcdf(file4)
		lon5, lat5, val5 = comp.input_netcdf(file5)

		print("delaunay range: %.2f %.2f %.2f %.2f " % ( min(lon1), max(lon1), min(lat1), max(lat1) ) );
		print("spline range: %.2f %.2f %.2f %.2f " % ( min(lon2), max(lon2), min(lat2), max(lat2) ) );
		print("visr range: %.2f %.2f %.2f %.2f " % ( min(lon3), max(lon3), min(lat3), max(lat3) ) );
		print("nd interp range: %.2f %.2f %.2f %.2f " % ( min(lon4), max(lon4), min(lat4), max(lat4) ) );
		print("tape range: %.2f %.2f %.2f %.2f " % ( min(lon5), max(lon5), min(lat5), max(lat5) ) );

		lons1, lats1, val1 = comp.confine_to_grid(lon1, lat1, val1, -124.38, -121.2, 37.2, 42, 0.04)
		lons2, lats2, val2 = comp.confine_to_grid(lon2, lat2, val2, -124.38, -121.2, 37.2, 42, 0.04)
		lons3, lats3, val3 = comp.confine_to_grid(lon3, lat3, val3, -124.38, -121.2, 37.2, 42, 0.04)
		lons4, lats4, val4 = comp.confine_to_grid(lon4, lat4, val4, -124.38, -121.2, 37.2, 42, 0.04)
		lons5, lats5, val5 = comp.confine_to_grid(lon5, lat5, val5, -124.38, -121.2, 37.2, 42, 0.04)

		comp.check_coregistration(val1, val2, val3, val4, val5);
		my_means = comp.grid_means(lons2, lats2, val1, val2, val3, val4, val5)
		my_sds = comp.grid_sds(lons2, lats2, val1, val2, val3, val4, val5)

		comp.output_nc(lons2, lats2, my_means, "means", "max_shear")
		comp.output_nc(lons2, lats2, my_sds, "deviations", "max_shear")

	else:
		print("%s is not a recognized strain component" % component)

	return;