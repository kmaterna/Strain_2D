import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr

from . import utilities, strain_tensor_toolbox, velocity_io, pygmt_plots


def drive(MyParams):
    """
    A driver for taking statistics of several strain computations
    """
    mean_dss = compare_grid_means(MyParams)
    visualize_grid_means(MyParams, mean_dss)


def compare_grid_means(MyParams, statistics_function=np.nanmean):
    """
    A driver for comparing strain rate maps

    Parameters
    ----------
    MyParams: dict            - Parameter dictionary
    statistics_function: func - standard numpy-compatible reducing function (e.g. mean, 
                                median, nanmedian)

    Returns
    -------
    mean_dss: xarray Dataset   - Dataset containing the mean and standard deviation of 
                                 each variable

    Writes
    ------
    mean_dss: xarray Dataset   - Dataset containing the mean and standard deviation of 
                                 each variable
    """
    ds_list = []
    for key, value in MyParams.strain_dict.items():
        specific_filename = glob.glob(value + os.sep + '*' + "_strain.nc")[0]
        ds = xr.load_dataset(specific_filename)
        ds_list.append(ds)

    # here we compute mean and standard deviation
    mean_dss = compute_grid_statistics(ds_list, statistics_function)
    mean_dss.to_netcdf(os.path.join(MyParams.outdir, "mean_std.nc"))

    # plot differences between methods for dilatation and max shear strain rates
    pygmt_plots.plot_method_differences(
        velocity_io.read_multiple_strain_netcdfs(MyParams, 'dilatation'),
        mean_dss['dilatation'],
        MyParams.range_strain, 
        MyParams.outdir,
        MyParams.outdir+"/separate_plots_"+"dilatation"+'.png',
    )
    pygmt_plots.plot_method_differences(
        velocity_io.read_multiple_strain_netcdfs(MyParams, 'max_shear'),
        mean_dss['max_shear'],
        MyParams.range_strain, 
        MyParams.outdir,
        MyParams.outdir+"/separate_plots_"+"max_shear"+'.png'
    )
    return mean_dss


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

def compute_grid_statistics(strain_grid_list, stat=np.nanmean):
    """
    A function that takes statistics on several mutually co-registered 
    grids in an xarray.DataSet. This function basically runs a loop.
    The inner function must return a mean-like value and a standard-deviation-like value
    Returns a dataset with two layers, mean and standard deviation
    """
    exx = stat([ds['exx'].data for ds in strain_grid_list], axis=0)
    eyy = stat([ds['eyy'].data for ds in strain_grid_list], axis=0)
    exy = stat([ds['exy'].data for ds in strain_grid_list], axis=0)
    rot = stat([ds['rotation'].data for ds in strain_grid_list], axis=0)
    Ve =  stat([ds['Ve'].data for ds in strain_grid_list], axis=0)
    Vn =  stat([ds['Vn'].data for ds in strain_grid_list], axis=0)

    exx_std = np.nanstd([ds['exx'].data for ds in strain_grid_list], axis=0)
    eyy_std = np.nanstd([ds['eyy'].data for ds in strain_grid_list], axis=0)
    exy_std = np.nanstd([ds['exy'].data for ds in strain_grid_list], axis=0)

    ms_std = np.nanstd([ds['max_shear'].data for ds in strain_grid_list], axis=0)
    dil_std = np.nanstd([ds['dilatation'].data for ds in strain_grid_list], axis=0)

    try:
        Se = stat([ds['Se'].data for ds in strain_grid_list], axis=0)
        Sn = stat([ds['Sn'].data for ds in strain_grid_list], axis=0)
    except KeyError:
        Se = np.empty(Ve.shape)
        Sn = np.empty(Vn.shape)

    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy)

    # Repacking result into DS
    mean_stds_ds = xr.Dataset(
        {
            "Ve":  (("y", "x"), Ve),
            "Vn":  (("y", "x"), Vn),
            "Se":  (("y", "x"), Se),
            "Sn":  (("y", "x"), Sn),
            "exx": (("y", "x"), exx),
            "eyy": (("y", "x"), eyy),
            "exy": (("y", "x"), exy),
            "rotation": (("y", "x"), rot),
            "max_shear": (("y", "x"), max_shear),
            "max_shear_sigma": (("y", "x"), ms_std),
            "dilatation": (("y", "x"), dilatation),
            "dilatation_sigma": (("y", "x"), dil_std),
            "I2": (("y", "x"), I2nd),
            "azimuth": (("y", "x"), azimuth),
            "exx_std": (("y", "x"), exx_std),
            "eyy_std": (("y", "x"), eyy_std),
            "exy_std": (("y", "x"), exy_std),
        },
        coords={
            "x": ('x', strain_grid_list[0].x.data),
            "y": ('y', strain_grid_list[0].y.data),
        },
    )

    
    
    return mean_stds_ds


