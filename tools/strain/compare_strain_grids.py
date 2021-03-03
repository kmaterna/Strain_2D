import os
from . import utilities, strain_tensor_toolbox
from Tectonic_Utils.read_write import netcdf_read_write
import numpy as np


def drive(MyParams):
    """
    A driver for taking statistics of several strain computations
    """
    compare_grid_means(MyParams, "I2nd.nc", log_means_statistics);
    compare_grid_means(MyParams, "max_shear.nc", simple_means_statistics);
    compare_grid_means(MyParams, "dila.nc", simple_means_statistics);
    compare_grid_means(MyParams, "rot.nc", simple_means_statistics);
    compare_grid_means(MyParams, "azimuth.nc", angular_means_statistics, mask=[MyParams.outdir+'I2nd.grd', 3]);
    return;


def compare_grid_means(MyParams, filename, statistics_function, mask=None):
    """
    A driver for taking the mean of several grid quantities
    The function for taking the mean/std is passed in
    `mask` has format [filename, cutoff_value] if you want to mask based on a particular computation result.
    """
    strain_values_dict = read_strain_files(MyParams, filename);
    lons, lats, my_means, my_stds = compute_grid_statistics(strain_values_dict, statistics_function);
    if mask:
        [_, _, masking_values] = netcdf_read_write.read_any_grd(mask[0]);
        my_means = utilities.mask_by_value(my_means, masking_values, mask[1]);
    netcdf_read_write.produce_output_netcdf(lons, lats, my_means, 'per year', MyParams.outdir+"/means_"+filename);
    netcdf_read_write.produce_output_netcdf(lons, lats, my_stds, 'per year', MyParams.outdir+"/deviations_"+filename);
    return;


# --------- READ FUNCTION ----------- #

def read_strain_files(MyParams, filename):
    """
    Read strain quantities of `filename` into a dictionary
    Each dictionary key is a strain method
    Each dictionary value is a data structure: [lon, lat, value]
    lon : list of floats
    lat: list of floats
    value: 2D array of floats
    We also guarantee the mutual co-registration of the dictionary elements
    """
    strain_values_dict = {};
    for method in MyParams.strain_dict.keys():
        specific_filename = MyParams.strain_dict[method]+"/"+filename
        if os.path.isfile(specific_filename):
            [lon, lat, val] = netcdf_read_write.read_any_grd(specific_filename);
            strain_values_dict[method] = [lon, lat, val];
        else:
            raise Exception("Error! Can't find file %s " % specific_filename);
    utilities.check_grids(MyParams, strain_values_dict);
    utilities.check_shapes(strain_values_dict);
    return strain_values_dict;


# --------- COMPUTE FUNCTION ----------- #

def compute_grid_statistics(strain_values_dict, statistic_function):
    """
    A function that takes statistics on several mutually co-registered grids.
    This function basically runs a loop.
    The inner function must return a mean-like value and a standard-deviation-like value
    The strain_values_dict has values that look like [lon, lat, data]
    """
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
            mean_val, sd_val = statistic_function(val_list);
            mean_vals[j][i] = mean_val
            sd_vals[j][i] = sd_val;
    return x, y, mean_vals, sd_vals;


def simple_means_statistics(value_list):
    """
    Take simple mean and standard deviation of a list of values
    """
    mean_val = np.nanmean(value_list)
    sd_val = np.nanstd(value_list)
    if mean_val == float("-inf"):
        mean_val = np.nan;
    return mean_val, sd_val;


def log_means_statistics(value_list):
    """
    Take mean and standard deviation of a list of values that are log quantities
    """
    value_list = [10 ** x for x in value_list];
    mean_val = np.nanmean(value_list);
    sd_val = np.nanstd(value_list)
    if mean_val != float("-inf"):
        mean_val = np.log10(mean_val)
    else:
        mean_val = np.nan
    sd_val = np.log10(sd_val)
    return mean_val, sd_val;


def angular_means_statistics(value_list):
    """
    Take mean and standard deviation of a list of values that are azimuths
    """
    mean_val, sd_val = np.nan, np.nan;
    theta, sd = strain_tensor_toolbox.angle_mean_math(value_list);
    if theta != float("-inf"):
        mean_val = theta
    if sd != float("inf"):
        sd_val = sd
    return mean_val, sd_val;
