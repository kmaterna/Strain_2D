# A set of code that reads multiple grid files and produces mean and variance statistics
# In grid form. 

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write


def defensive_programming(MyParams, strain_values_dict):
    for method in strain_values_dict.keys():
        print("%s range: %.2f %.2f %.2f %.2f " % (method, np.min(strain_values_dict[method][0]),
                                                  np.max(strain_values_dict[method][0]),
                                                  np.min(strain_values_dict[method][1]),
                                                  np.max(strain_values_dict[method][1])) );
    check_grids(MyParams, strain_values_dict);
    check_shapes(strain_values_dict);
    return;


# this function makes sure arrays are of the same dimensions before attempting to produce any statistics
def check_shapes(strain_values_dict):
    method1 = list(strain_values_dict.keys())[0];
    for method2 in strain_values_dict.keys():
        if method2 == method1:
            continue;
        v1 = strain_values_dict[method1][2];  # strain values array
        v2 = strain_values_dict[method2][2];  # second strain values array
        if np.shape(v1) != np.shape(v2):
            print("\n   Oops! The shape of %s is %d by %d " % (method1, np.shape(v1)[0], np.shape(v1)[1] ) );
            print("   But shape of %s is %d by %d " % (method2, np.shape(v2)[0], np.shape(v2)[1] ) );
            print("   so not all value arrays are coregistered! \n")
            raise Exception("Error! Not all arrays are coregistered!  Cannot compare.")
    print("All methods are have the same shape!")
    return;


def check_grids(MyParams, strain_values_dict):
    range_strain = np.round(MyParams.range_strain, 6);
    inc = np.round(MyParams.inc, 6);
    for method in strain_values_dict.keys():
        if np.round(np.min(strain_values_dict[method][0]), 6) != range_strain[0]:
            raise Exception("Strain %s doesn't match range_strain longitude" % method);
        if np.round(np.max(strain_values_dict[method][0]), 6) != range_strain[1]:
            raise Exception("Strain %s doesn't match range_strain longitude" % method);
        if np.round(np.min(strain_values_dict[method][1]), 6) != range_strain[2]:
            raise Exception("Strain %s doesn't match range_strain latitude" % method);
        if np.round(np.max(strain_values_dict[method][1]), 6) != range_strain[3]:
            raise Exception("Strain %s doesn't match range_strain latitude" % method);
        if np.round(strain_values_dict[method][0][1] - strain_values_dict[method][0][0], 6) != inc[0]:
            print("Inc %s : %f " % (method, np.round(strain_values_dict[method][0][1] - strain_values_dict[method][0][0], 6)) );
            print("Inc config: %f" % inc[0]);
            raise Exception("Increment %s doesn't match inc longitude" % method);
        if np.round(strain_values_dict[method][1][1] - strain_values_dict[method][1][0], 6) != inc[1]:
            print("Inc %s : %f " % (method, np.round(strain_values_dict[method][1][1] - strain_values_dict[method][1][0], 6)) );
            print("Inc config: %f" % inc[1]);
            raise Exception("Increment %s doesn't match inc latitude" % method);
    print("All methods have the same range/increment!")
    return;


def grid_means_stds(strain_values_dict):
    # calculate grid-wise mean and standard deviation, returning arrays with dimension latitude by longitude
    first_key = list(strain_values_dict.keys())[0]
    x = strain_values_dict[first_key][0];
    y = strain_values_dict[first_key][1];
    mean_vals = np.nan * np.ones([len(y), len(x)])
    sd_vals = np.nan * np.ones([len(y), len(x)])
    for j in range(len(y)):
        for i in range(len(x)):
            val_list = [];
            for method in strain_values_dict.keys():
                val_list.append(strain_values_dict[method][2][j][i]);
            mean_val = np.nanmean(val_list)
            sd_val = np.nanstd(val_list)
            if mean_val != float("-inf"):
                mean_vals[j][i] = mean_val
            sd_vals[j][i] = sd_val
    return mean_vals, sd_vals


def grid_means_log(strain_values_dict):
    # calculate grid-wise mean on a log quantity, returning array with dimension latitude by longitude
    first_key = list(strain_values_dict.keys())[0]
    x = strain_values_dict[first_key][0];
    y = strain_values_dict[first_key][1];
    mean_vals = np.nan * np.ones([len(y), len(x)])
    sd_vals = np.nan * np.ones([len(y), len(x)])
    for j in range(len(y)):
        for i in range(len(x)):
            val_list = [];
            for method in strain_values_dict.keys():
                val_list.append(10**strain_values_dict[method][2][j][i]);
            mean_val = np.nanmean(val_list);
            sd_val = np.nanstd(val_list)
            if mean_val != float("-inf"):
                mean_vals[j][i] = np.log10(mean_val)
            sd_vals[j][i] = np.log10(sd_val)
    return mean_vals, sd_vals


def angle_means(strain_values_dict):
    # Implementing the angular mean formulas
    first_key = list(strain_values_dict.keys())[0]
    x = strain_values_dict[first_key][0];
    y = strain_values_dict[first_key][1];
    mean_vals = np.nan * np.ones([len(y), len(x)])
    sd_vals = np.nan * np.ones([len(y), len(x)])
    for j in range(len(y)):
        for i in range(len(x)):
            azimuth_values = [];
            for method in strain_values_dict.keys():
                azimuth_values.append(strain_values_dict[method][2][j][i]);
            theta, sd = angle_mean_math(azimuth_values);
            if theta != float("-inf"):
                mean_vals[j][i] = theta
            if sd != float("inf"):
                sd_vals[j][i] = sd
    return mean_vals, sd_vals


def angle_mean_math(azimuth_values):
    # Angles in degrees
    # separated out so we can unit-test this math
    sin_list, cos_list = [], [];
    for phi in azimuth_values:
        sin_list.append(np.sin(2 * np.radians(90 - phi)));
        cos_list.append(np.cos(2 * np.radians(90 - phi)));
    s = np.nanmean(sin_list);
    c = np.nanmean(cos_list);
    R = ((s ** 2 + c ** 2) ** .5)
    V = 1 - R
    sd = np.degrees((-2 * np.log(R)) ** .5) / 2
    # sd = np.degrees((2*V)**.5)
    # t = np.arctan2(s, c)
    # strike = R*math.e**(math.i*t)
    strike = np.arctan2(s, c) / 2
    theta = 90 - np.degrees(strike)
    if theta < 0:
        theta = 180 + theta
    elif theta > 180:
        theta = theta - 180
    return theta, sd;


def mask_by_value(outdir, grid1, grid2, cutoff_value):
    # grid1 = usually azimuth deviations
    # grid2 = usually I2nd
    lon1, lat1, val1 = netcdf_read_write.read_any_grd(outdir+"/deviations_"+grid1+".nc");
    lon2, lat2, val2 = netcdf_read_write.read_any_grd(outdir+"/means_"+grid2+".nc");
    masked_vals = np.zeros(np.shape(val2));
    for i in range(len(lon1)):
        for j in range(len(lat1)):
            if abs(val2[j][i]) > cutoff_value:
                masked_vals[j][i] = val1[j][i];
            else:
                masked_vals[j][i] = np.nan;
    netcdf_read_write.produce_output_netcdf(lon1, lat1, masked_vals, 'per yr', outdir+"/deviations_"+grid1+".nc");
    return;
