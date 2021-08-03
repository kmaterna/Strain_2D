import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

from . import utilities, strain_tensor_toolbox, velocity_io, pygmt_plots


def drive(MyParams):
    """
    A driver for taking statistics of several strain computations
    """
    mean_ds = xr.Dataset()
    mean_ds['max_shear'] = compare_grid_means(MyParams, "max_shear.nc", simple_means_statistics)
    mean_ds['dilatation'] = compare_grid_means(MyParams, "dilatation.nc", simple_means_statistics)
    mean_ds['I2'] = compare_grid_means(MyParams, "I2.nc", log_means_statistics)
    mean_ds['rotation'] = compare_grid_means(MyParams, "rotation.nc", simple_means_statistics)
    mean_ds['azimuth'] = compare_grid_means(MyParams, "azimuth.nc", angular_means_statistics, mask=[MyParams.outdir+'/means_I2.nc', 3])
    visualize_grid_means(MyParams, mean_ds)


def compare_grid_means(MyParams, filename, statistics_function, mask=None):
    """
    A driver for taking the mean of several grid quantities
    The function for taking the mean/std is passed in
    `mask` has format [filename, cutoff_value] if you want to mask based on a particular computation result.
    """
    strain_values_ds = velocity_io.read_multiple_strain_files(MyParams, filename.split('.')[0]);
    mean_ds = strain_values_ds.to_array(dim='new').reduce(np.nanmean, dim='new')
    std_ds = strain_values_ds.to_array(dim='new').reduce(np.nanstd, dim='new')
    mean_ds.to_netcdf(os.path.join(MyParams.outdir, "means_"+filename))
    std_ds.to_netcdf(os.path.join(MyParams.outdir, "devations_"+filename))
    if "dila" in filename or "max_shear" in filename:
        pygmt_plots.plot_method_differences(
            strain_values_ds,
            mean_ds, 
            MyParams.range_strain, 
            MyParams.outdir,
            MyParams.outdir+"/separate_plots_"+filename.split('.')[0]+'.png'
        )
    return mean_ds


def visualize_grid_means(MyParams, ds):
    """ Make pygmt plots of the means of all quantities """
    pygmt_plots.plot_I2nd(ds['I2'], MyParams.range_strain, MyParams.outdir, [], [],
                          MyParams.outdir + "/means_I2nd.png");
    pygmt_plots.plot_dilatation(ds['dilatation'], MyParams.range_strain, MyParams.outdir, [], [],
                                MyParams.outdir + "/means_dila.png");
    pygmt_plots.plot_maxshear(ds['max_shear'], MyParams.range_strain, MyParams.outdir, [], [],
                              MyParams.outdir + "/means_max_shear.png");
    pygmt_plots.plot_azimuth(ds['azimuth'], MyParams.range_strain, MyParams.outdir, [], [],
                             MyParams.outdir + "/means_azimuth.png");
    pygmt_plots.plot_rotation(ds['rotation'], [], MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir + "/means_rot.png");
    plt.close('all') # clear the memory cache


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
