import numpy as np
import netCDF4
import netcdf_functions

# takes in netcdfs from each method of strain calculation, and produces netcdfs of the means and standard deviations

# directories to netcdfs for each method
file1 = "../Results/Results_Delaunay/I2nd.nc"
file2 = "../Results/Results_Numpy_Spline/I2nd.nc"
file3 = "../Results/Results_Visr/I2nd.nc"
file4 = "../Results/Results_ND_interp/I2nd.nc"

# inputs data from netcdfs. netcdfs must have uniform grid size, and generally cover northern california
# outputs an n-dim array of lons, and m-d array of lats, and an m by n array of values
def input_netcdf(nc):
	file = netCDF4.Dataset(nc, 'r', Format='NETCDF4')
	lon = file.variables['x'][:]
	lat = file.variables['y'][:]
	values = file.variables['z']
	return lon, lat, values

# uses uniform data points from net cdfs with different coord boxes, and confines it to a uniform coord box
def confine_to_grid(x, y, values, xmin, xmax, ymin, ymax, inc):
	new_lon = np.arange(xmin, xmax, inc)
	new_lat = np.arange(ymin, ymax, inc)
	new_vals = []
	for i in range(len(x)):
		for j in range(len(y)):
			if x[i] >= xmin and x[i] <= xmax and y[j] >= ymin and y[j] <= ymax:
				new_vals.append(values[j][i])
	print(len(new_vals))
	final_vals = []
	for i in np.arange(0, len(new_vals), len(new_lat)):
		final_vals.append(new_vals[i:i+len(new_lat)])
	final_vals = np.transpose(final_vals)
	return new_lon, new_lat, final_vals

# gridwise, calculates means and standard deviations, and returns them as arrays with dimension latitude by longitude
def grid_avg_std(x, y, vals1, vals2, vals3, vals4):
	mean_vals = np.nan * np.ones([len(y), len(x)])
	sd_vals = np.nan * np.ones([len(y), len(x)])
	for j in range(len(y)):
		for i in range(len(x)):
			val1 = vals1[j][i]
			val2 = vals2[j][i]
			val3 = vals3[j][i]
			val4 = vals4[j][i]
			mean_val = np.nanmean([val1, val2, val3, val4])
			sd_val = np.nanstd([val1, val2, val3, val4])
			if mean_val != float("-inf"):
				mean_vals[j][i] = mean_val
			sd_vals[j][i] = sd_val
	return mean_vals, sd_vals

# writes a netcdf from the uniform lat, lon, and the means or stardard deviations
def output_nc(lon, lat, vals, file):
	netcdf_functions.produce_output_netcdf(lon, lat, vals, 'per yr', "../Results/Results_means/"+file+".nc");
	return

	

lon1, lat1, val1 = input_netcdf(file1)
lon2, lat2, val2 = input_netcdf(file2)
lon3, lat3, val3 = input_netcdf(file3)
lon4, lat4, val4 = input_netcdf(file4)

lons1, lats1, val1 = confine_to_grid(lon1, lat1, val1, -124.5, -121.5, 39, 42, 0.04)
lons2, lats2, val2 = confine_to_grid(lon2, lat2, val2, -124.5, -121.5, 39, 42, 0.04)
lons3, lats3, val3 = confine_to_grid(lon3, lat3, val3, -124.5, -121.5, 39, 42, 0.04)
lons4, lats4, val4 = confine_to_grid(lon4, lat4, val4, -124.5, -121.5, 39, 42, 0.04)

my_means, my_sds = grid_avg_std(lons2, lats2, val1, val2, val4, val4)

output_nc(lons2, lats2, my_means, "means")
output_nc(lons2, lats2, my_sds, "deviations")
