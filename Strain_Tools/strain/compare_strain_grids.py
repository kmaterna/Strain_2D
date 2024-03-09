import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

from . import utilities, strain_tensor_toolbox, velocity_io, pygmt_plots


def drive(MyParams):
    """
    A driver for taking statistics of several strain computations
    """
    mean_dss = xr.Dataset()
    mean_dss['max_shear'] = compare_grid_means(MyParams, "max_shear", simple_means_statistics)
    mean_dss['dilatation'] = compare_grid_means(MyParams, "dilatation", simple_means_statistics)
    mean_dss['I2'] = compare_grid_means(MyParams, "I2", simple_means_statistics)
    mean_dss['rotation'] = compare_grid_means(MyParams, "rotation", simple_means_statistics)
    mean_dss['azimuth'] = compare_grid_means(MyParams, "azimuth", angular_means_statistics,
                                             mask=[MyParams.outdir+'/means_I2.nc', 3])
    visualize_grid_means(MyParams, mean_dss)


def compare_grid_means(MyParams, plot_type, statistics_function, mask=None):
    """
    A driver for comparing strain rate maps

    Parameters
    ----------
    MyParams: dict            - Parameter dictionary
    plot_type: str            - Type of strain quantity to compare 
    statistics_function: func - standard numpy-compatible reducing function (e.g. mean, median, nanmedian)
    mask:                     - length-2 list of [filename, cutoff_value] used for thresholding the plot_type

    Returns
    -------
    mean_stds_ds: xarray Dataset   - Dataset containing the mean and standard deviation of each variable

    Writes
    ------
    mean_ds, std_ds: xarray Dataset - writes these to NETCDF
    """
    # here we extract each grid of plot_type into an xarray.Dataset
    strain_values_ds = velocity_io.read_multiple_strain_netcdfs(MyParams, plot_type)

    # here we compute mean and standard deviation
    mean_stds_ds = compute_grid_statistics(strain_values_ds, statistics_function)

    mean_stds_ds.to_netcdf(os.path.join(MyParams.outdir, "means_stds_"+str(plot_type)+".nc"))
    if "dila" in plot_type or "max_shear" in plot_type:
        pygmt_plots.plot_method_differences(
            strain_values_ds,
            mean_stds_ds['mean'],
            MyParams.range_strain, 
            MyParams.outdir,
            MyParams.outdir+"/separate_plots_"+plot_type.split('.')[0]+'.png'
        )
    return mean_stds_ds['mean']


def visualize_grid_means(MyParams, ds):
    """ Make pygmt plots of the means of all quantities """
    pygmt_plots.plot_I2nd(ds['I2'], (), MyParams.range_strain, MyParams.outdir, MyParams.outdir + "/means_I2nd.png")
    pygmt_plots.plot_dilatation(ds['dilatation'], (), MyParams.range_strain, MyParams.outdir,
                                MyParams.outdir + "/means_dila.png")
    pygmt_plots.plot_maxshear(ds['max_shear'], (), MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir + "/means_max_shear.png")
    pygmt_plots.plot_azimuth(ds['azimuth'], (), MyParams.range_strain, MyParams.outdir,
                             MyParams.outdir + "/means_azimuth.png")
    pygmt_plots.plot_rotation(ds['rotation'], [], MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir + "/means_rot.png")
    plt.close('all')  # clear the memory cache


# --------- COMPUTE FUNCTION ----------- #

def compute_grid_statistics(strain_values_ds, statistic_function):
    """
    A function that takes statistics on several mutually co-registered grids in an xarray.DataSet.
    This function basically runs a loop.
    The inner function must return a mean-like value and a standard-deviation-like value
    Returns a dataset with two layers, mean and standard deviation
    """

    x = np.array(strain_values_ds['x'])
    y = np.array(strain_values_ds['y'])
    num_grids = len(strain_values_ds.data_vars.items())

    # Unpacking into 3D numpy array
    comparative_strain_values = np.zeros((len(y), len(x), num_grids))
    for i, (varname, da) in enumerate(strain_values_ds.data_vars.items()):
        comparative_strain_values[:, :, i] = np.array(da)

    mean_vals = np.nan * np.ones([len(y), len(x)])
    sd_vals = np.nan * np.ones([len(y), len(x)])
    for j in range(len(y)):
        for i in range(len(x)):
            mean_vals[j][i], sd_vals[j][i] = statistic_function(comparative_strain_values[j][i][:])

    # Repacking result into DS
    mean_stds_ds = xr.Dataset(
        {
            "mean": (("y", "x"), mean_vals),
            "stds": (("y", "x"), sd_vals),
        },
        coords={
            "x": ('x', x),
            "y": ('y', y),
        },
    )
    return mean_stds_ds


def simple_means_statistics(value_list):
    """
    Take simple mean and standard deviation of a list of values
    """
    mean_val = np.nanmean(value_list)
    sd_val = np.nanstd(value_list)
    if mean_val == float("-inf"):
        mean_val = np.nan
    return mean_val, sd_val


def log_means_statistics(value_list):
    """
    Take mean and standard deviation of a list of values that are log quantities
    """
    value_list = [10 ** x for x in value_list]
    mean_val = np.nanmean(value_list)
    sd_val = np.nanstd(value_list)
    if mean_val != float("-inf"):
        mean_val = np.log10(mean_val)
    else:
        mean_val = np.nan
    sd_val = np.log10(sd_val)
    return mean_val, sd_val


def angular_means_statistics(value_list):
    """
    Take mean and standard deviation of a list of values that are azimuths
    """
    mean_val, sd_val = np.nan, np.nan
    theta, sd = strain_tensor_toolbox.angle_mean_math(value_list)
    if theta != float("-inf"):
        mean_val = theta
    if sd != float("inf"):
        sd_val = sd
    return mean_val, sd_val
