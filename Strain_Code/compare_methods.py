import numpy as np
import netCDF4
import netcdf_functions

# takes in netcdfs from each method of strain calculation, and produces netcdfs of the means and standard deviations

# directories to netcdfs for each method
file1 = "../Results/Results_Delaunay/I2nd.nc"
file2 = "../Results/Results_Numpy_Spline/I2nd.nc"
file3 = "../Results/Results_Visr/I2nd.nc"
file4 = "../Results/Results_ND_interp/I2nd.nc"
file5 = "../Results/Results_Tape/I2nd.nc"

# inputs data from netcdfs. netcdfs must have uniform grid size, and generally cover northern california
# outputs an n-dim array of lons, and m-d array of lats, and an m by n array of values
def input_netcdf(nc):
	file = netCDF4.Dataset(nc, 'r', Format='NETCDF4')
	lon = file.variables['x'][:]
	lat = file.variables['y'][:]
	values = file.variables['z']
	return lon, lat, values

# uses uniform data points from net cdfs with different coord boxes, and confines it to a uniform coord box
# all points must be coregistered even if boxes are different sizes
# i.e. the starting values must be different by multiples of the incriment
# this is controlled by configure_functions and produce_gridded, depending on the method.
def confine_to_grid(x, y, values, xmin, xmax, ymin, ymax, inc):
	print("registering %s to grid [%.2f %.2f %.2f %.2f %.2f] " % ("strain", xmin, xmax, ymin, ymax, inc) );

	new_lon = np.arange(xmin, xmax, inc)
	new_lat = np.arange(ymin, ymax, inc)
	new_vals = []

	x = np.around(x, 5)
	y = np.around(y, 5)
	tol = 0.002

	for i in range(len(x)):
		for j in range(len(y)):
			if (xmin < x[i] < xmax) or (abs((x[i] - xmin)) < tol) or (abs((x[i] - xmax)) < tol):
				if (ymin < y[j] < ymax) or  (abs((y[j] - ymin)) < tol) or (abs((y[j] - ymax)) < tol):
					new_vals.append(values[j][i])

	final_vals = []
	for i in np.arange(0, len(new_vals), len(new_lat)):
		final_vals.append(new_vals[i:i+len(new_lat)])
	final_vals = np.transpose(final_vals)

	return new_lon, new_lat, final_vals

# this function makes sure arrays are of the same dimensions before attempting to produce any statistics
def check_coregistration(v1, v2, v3, v4, v5):
	if np.shape(v1) != np.shape(v2):
		print("\n   Oops! The shape of method 1 vals is %d by %d \n" % (np.shape(v1)[0], np.shape(v1)[1] ) );
		print("   But shape of method 2 vals is %d by %d " % (np.shape(v2)[0], np.shape(v2)[1] ) );
		print("   so not all value arrays are coregistered! \n")
	elif np.shape(v1) != np.shape(v3):
		print("\n   Oops! The shape of method 1 vals is %d by %d \n" % (np.shape(v1)[0], np.shape(v1)[1] ) );
		print("   But shape of method 3 vals is %d by %d \n" % (np.shape(v3)[0], np.shape(v3)[1] ) );
		print("   so not all value arrays are coregistered! \n")
	elif np.shape(v1) != np.shape(v4):
		print("\n   Oops! The shape of method 1 vals is %d by %d \n" % (np.shape(v1)[0], np.shape(v1)[1] ) );
		print("   But shape of method 4 vals is %d by %d \n" % (np.shape(v4)[0], np.shape(v4)[1] ) );
		print("   so not all value arrays are coregistered! \n")
	elif np.shape(v1) != np.shape(v5):
		print("\n   Oops! The shape of method 1 vals is %d by %d \n" % (np.shape(v1)[0], np.shape(v1)[1] ) );
		print("   But shape of method 5 vals is %d by %d " % (np.shape(v5)[0], np.shape(v5)[1] ) );
		print("   so not all value arrays are coregistered! \n")
	else:
		print("All methods are coregistered!")
	return

# gridwise, calculates means and standard deviations, and returns them as arrays with dimension latitude by longitude
def grid_avg_std(x, y, vals1, vals2, vals3, vals4, vals5):
	mean_vals = np.nan * np.ones([len(y), len(x)])
	sd_vals = np.nan * np.ones([len(y), len(x)])
	for j in range(len(y)):
		for i in range(len(x)):
			val1 = vals1[j][i]
			val2 = vals2[j][i]
			val3 = vals3[j][i]
			val4 = vals4[j][i]
			val5 = vals5[j][i]
			mean_val = np.nanmean([val1, val2, val3, val4, val5])
			sd_val = np.nanstd([val1, val2, val3, val4, val5])
			if mean_val != float("-inf"):
				mean_vals[j][i] = mean_val
			sd_vals[j][i] = sd_val
	return mean_vals, sd_vals

# writes the uniform latitude, longitude, and whichever statistical values are desired.
# outputs to result directory for means as a netcdf which can then be manipulated further with gmt.
def output_nc(lon, lat, vals, file):
	netcdf_functions.produce_output_netcdf(lon, lat, vals, 'per yr', "../Results/Results_means/"+file+".nc");
	return


# calling all functions on all methods:

lon1, lat1, val1 = input_netcdf(file1)
lon2, lat2, val2 = input_netcdf(file2)
lon3, lat3, val3 = input_netcdf(file3)
lon4, lat4, val4 = input_netcdf(file4)
lon5, lat5, val5 = input_netcdf(file5)

print("delaunay range: %.2f %.2f %.2f %.2f " % ( min(lon1), max(lon1), min(lat1), max(lat1) ) );

print("spline range: %.2f %.2f %.2f %.2f " % ( min(lon2), max(lon2), min(lat2), max(lat2) ) );

print("visr range: %.2f %.2f %.2f %.2f " % ( min(lon3), max(lon3), min(lat3), max(lat3) ) );

print("nd interp range: %.2f %.2f %.2f %.2f " % ( min(lon4), max(lon4), min(lat4), max(lat4) ) );

print("tape range: %.2f %.2f %.2f %.2f " % ( min(lon5), max(lon5), min(lat5), max(lat5) ) );

lons1, lats1, val1 = confine_to_grid(lon1, lat1, val1, -124.3, -121.2, 37.2, 42, 0.04)
lons2, lats2, val2 = confine_to_grid(lon2, lat2, val2, -124.3, -121.2, 37.2, 42, 0.04)
lons3, lats3, val3 = confine_to_grid(lon3, lat3, val3, -124.3, -121.2, 37.2, 42, 0.04)
lons4, lats4, val4 = confine_to_grid(lon4, lat4, val4, -124.3, -121.2, 37.2, 42, 0.04)
lons5, lats5, val5 = confine_to_grid(lon5, lat5, val5, -124.3, -121.2, 37.2, 42, 0.04)

check_coregistration(val1, val2, val3, val4, val5);
my_means, my_sds = grid_avg_std(lons2, lats2, val1, val3, val4, val4, val5)

output_nc(lons2, lats2, my_means, "means")
output_nc(lons2, lats2, my_sds, "deviations")
