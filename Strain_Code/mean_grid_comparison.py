import compare_methods as comp 



def drive(component):

	velfield="midas"
	outdir="Results/"+velfield+"/Means/"
	component_name=component;
	if component=="dilatation":
		component_name="dila";

	print("Comparing %s across all methods" % component_name)
	if component == "I2nd" or component =="dilatation" or component=="max_shear" or component=="azimuth":
		# directories to netcdfs for each method
		file1 = "Results/"+velfield+"/Results_hammond/"+component_name+".nc"
		file2 = "Results/"+velfield+"/Results_spline/"+component_name+".nc"
		file3 = "Results/"+velfield+"/Results_visr/"+component_name+".nc"
		file4 = "Results/"+velfield+"/Results_ND_interp/"+component_name+".nc"
		file5 = "Results/"+velfield+"/Results_tape/"+component_name+".nc"

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
		comp.output_nc(lons2, lats2, my_means, outdir, "means", component_name)
		comp.output_nc(lons2, lats2, my_sds, outdir, "deviations", component_name)

	elif component == "uplift": # Rarely used. 
		file1 = "Results/"+velfield+"/Results_hammond/"+component_name+".txt"
		file2 = "Results/"+velfield+"/Results_spline/"+component_name+".txt"
		file3 = "Results/"+velfield+"/Results_visr/"+component_name+".txt"
		file4 = "Results/"+velfield+"/Results_ND_interp/"+component_name+".txt"
		file5 = "Results/"+velfield+"/Results_tape/"+component_name+".txt"

		lon1, lat1, val1 = comp.input_txt(file1)
		lon2, lat2, val2 = comp.input_txt(file2)
		lon3, lat3, val3 = comp.input_txt(file3)
		lon4, lat4, val4 = comp.input_txt(file4)
		lon5, lat5, val5 = comp.input_txt(file5)

		my_means = comp.array_means(lon1, lat1, val1, val2, val4)
		comp.output_txt(lon1, lat1, my_means, outdir, "means", component_name);
		print("Max uplift: %f mm/yr" % max(my_means))
		print("Min uplift: %f mm/yr" % min(my_means))

	else:
		print("%s is not a recognized strain component" % component)

	return;
