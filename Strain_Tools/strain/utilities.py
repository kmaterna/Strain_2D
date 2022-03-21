# A set of utility functions used throughout the Strain_2D library
import subprocess
import numpy as np
import xarray as xr

from . import velocity_io


def get_float_range(string_range):
    """
    :param string_range: format "-125/-121/32/35"
    :type string_range: string
    :returns: list of floats
    :rtype: list
    """
    number_strings = string_range.split('/')
    float_range = [float(number_strings[0]), float(number_strings[1]),
                   float(number_strings[2]), float(number_strings[3])];
    if float_range[1] <= float_range[0]:
        raise ValueError("Error! Given range is invalid", float_range);
    if float_range[3] <= float_range[2]:
        raise ValueError("Error! Given range is invalid", float_range);
    return float_range;


def get_string_range(float_range, x_buffer=0, y_buffer=0):
    """
    Buffer is for the possible interface between pixel-node-registered and gridline-node-registered files

    :param float_range: list, [w, e, s, n]
    :type float_range: list
    :param x_buffer: possible interface between pixel-node-registered etc.
    :type x_buffer: float
    :param y_buffer: possible interface between pixel-node-registered etc.
    :type y_buffer: float
    :returns: string range
    :rtype: string
    """
    string_range = str(float_range[0] - x_buffer) + '/' + str(float_range[1] + x_buffer) + '/' + \
                   str(float_range[2] - y_buffer) + '/' + str(float_range[3] + y_buffer);
    return string_range;


def get_float_inc(string_inc):
    """
    :param string_inc: string, e.g., '0.04/0.04'
    :type string_inc: string
    :returns: list of floats
    :rtype: list
    """
    number_incs = string_inc.split('/')
    float_inc = [float(number_incs[0]), float(number_incs[1])];
    return float_inc;


def get_string_inc(float_inc):
    """
    :type float_inc: list
    :returns: string separated by slash, e.g., '0.04/0.04'
    :rtype: string
    """
    string_inc = str(float_inc[0]) + '/' + str(float_inc[1]);
    return string_inc;


def get_gmt_range_inc(lons, lats):
    """
    Take lons and lats associated with pixel-node-registered files read into Python using xarray.ds

    :param lons: list of pixel centers lons
    :type lons: np.array
    :param lats: list of pixel centers lats
    :type lats: np.array
    :returns: string range, string inc
    :rtype: string, string
    """
    lon_inc = np.round(lons[1] - lons[0], 6)
    edge_of_west_pixel = np.round(lons[0] - lon_inc/2, 5);
    edge_of_east_pixel = np.round(lons[-1] + lon_inc/2, 5);
    lat_inc = np.round(lats[1] - lats[0], 6);
    edge_of_south_pixel = np.round(lats[0] - lat_inc/2, 5);
    edge_of_north_pixel = np.round(lats[-1] + lat_inc/2, 5);
    gmt_range_string = str(edge_of_west_pixel) + '/' + str(edge_of_east_pixel) + '/' + str(edge_of_south_pixel) + \
                       '/' + str(edge_of_north_pixel);
    gmt_inc_string = str(lon_inc) + '/' + str(lat_inc);
    return gmt_range_string, gmt_inc_string;


# --------- DEFENSIVE PROGRAMMING FOR COMPARING MULTIPLE GRIDS ------------------ #

def check_coregistered_shapes(strain_values_ds):
    """
    Make sure arrays are of the same dimensions before attempting to produce any statistics

    :param strain_values_ds: xarray.DataSet of strain values from different calculations
    :returns: None
    """
    x = np.array(strain_values_ds['x']);
    y = np.array(strain_values_ds['y']);
    for varname, da in strain_values_ds.data_vars.items():
        nparray = np.array(strain_values_ds[varname]);
        arrayshape = np.shape(nparray);
        assert (arrayshape == (len(y), len(x)), ValueError(
            "Error! Not all arrays have the same shape!  Cannot compare."));
    print("All methods have the same shape.");
    return;


# --------- GRID AND VELFIELD NAMED TUPLE UTILITIES ------------------ #

def make_grid(coordbox, inc):
    """
    Assumption is a pixel-node-registered grid.
    :param coordbox: [float, float, float, float] corresponding to [W, E, S, N]
    :type coordbox: list
    :param inc: [float, float] corresponding to [xinc, yinc]
    :type inc: list
    :returns: 1d array of lons, 1d array of lats, 2d array of zeros
    """
    lonmin, lonmax = coordbox[0], coordbox[1]
    latmin, latmax = coordbox[2], coordbox[3]
    lons = np.arange(lonmin, lonmax+0.00001, inc[0])
    lats = np.arange(latmin, latmax+0.00001, inc[1])
    grid = np.zeros((len(lats), len(lons)));
    return lons, lats, grid


def getVels(velField):
    """Extract velocities from a list of NamedTuples"""
    lon = np.array([x.elon for x in velField]);
    lat = np.array([x.nlat for x in velField]);
    e = np.array([x.e for x in velField]);
    n = np.array([x.n for x in velField]);
    se = np.array([x.se for x in velField]);
    sn = np.array([x.se for x in velField]);
    return lon, lat, e, n, se, sn;


def get_index_of_nearest_point(xvals, target_val):
    idx = (np.abs(xvals - target_val)).argmin()
    return idx;


def subtract_two_velfields(obsfield, modelfield):
    # Perform residual subtraction
    residual_velfield = [];
    for obs, model in zip(obsfield, modelfield):
        newVelPoint = velocity_io.StationVel(elon=obs.elon, nlat=obs.nlat, e=obs.e - model.e,
                                             n=obs.n - model.n, u=0, se=0, sn=0, su=0, name=obs.name);
        residual_velfield.append(newVelPoint);
    return residual_velfield;


def create_model_velfield(xdata, ydata, Ve, Vn, myVelfield):
    """xdata, ydata are 1D arrays of lon and lat.  Ve, Vn are 2d arrays of interpolated velocities
    Returns MODEL velocities at points in myVelfield. """
    model_velfield = [];
    for item in myVelfield:
        target_lon = item.elon;
        target_lat = item.nlat;
        lon_idx = get_index_of_nearest_point(xdata, target_lon);
        lat_idx = get_index_of_nearest_point(ydata, target_lat);
        modelPoint = velocity_io.StationVel(elon=item.elon, nlat=item.nlat, e=Ve[lat_idx, lon_idx],
                                            n=Vn[lat_idx, lon_idx], u=0, se=0, sn=0, su=0, name=item.name);
        model_velfield.append(modelPoint);
    return model_velfield;


def filter_by_bounding_box(velfield, bbox):
    filtered_velfield = [];
    for item in velfield:
        if bbox[0] <= item.elon <= bbox[1]:
            if bbox[2] <= item.nlat <= bbox[3]:
                filtered_velfield.append(item);
    return filtered_velfield;


# --------- GRD/NETCDF UTILITIES ------------------ #

def mask_by_value(grid1, grid_maskingbasis, cutoff_value):
    """
    Implement NAN-mask for one grid in all places where values are smaller than cutoff value in corresponding grid.

    :param grid1: usually azimuth deviations
    :type grid1: 2D array
    :param grid_maskingbasis: usually I2nd
    :type grid_maskingbasis: 2D array
    :param cutoff_value: cutoff for nans
    :type cutoff_value: float
    :returns: masked grid
    :rtype: 2D array
    """
    (y, x) = np.shape(grid1);
    masked_vals = np.zeros(np.shape(grid1));
    for i in range(x):
        for j in range(y):
            if abs(grid_maskingbasis[j][i]) > cutoff_value:
                masked_vals[j][i] = grid1[j][i];
            else:
                masked_vals[j][i] = np.nan;
    return masked_vals;


def read_basic_fields_from_netcdf(netcdf_name):
    """
    Reads x, y, and strain from formatted NetCDF output for Strain project.
    Return type is class 'xarray.core.dataarray.DataArray' for all 1D and 2D arrays
    """
    ds = xr.open_dataset(netcdf_name);
    lons = ds["x"]
    lats = ds["y"]
    exx = np.reshape(ds['exx'], (len(lats), len(lons)));
    exy = np.reshape(ds['exy'], (len(lats), len(lons)));
    eyy = np.reshape(ds['eyy'], (len(lats), len(lons)));
    return lons, lats, exx, exy, eyy;


def make_gmt_landmask(lons, lats, grd_filename):
    """
    Use GMT to construct a landmask for calculating total moment accumulation from strain rate.

    param lons: np.array of centers of pixels
    param lats: np.array of centers of pixels
    param filename: string
    returns landmask array
    """
    gmt_range_string, gmt_inc_string = get_gmt_range_inc(np.array(lons), np.array(lats));
    subprocess.call(['gmt', 'grdlandmask', '-G'+grd_filename, '-R'+gmt_range_string, '-I'+gmt_inc_string, '-r'],
                    shell=False);  # guarantee pixel node registration
    print('gmt grdlandmask -G'+grd_filename+' -R'+gmt_range_string+'-I'+gmt_inc_string+' -r');
    landmask_array = read_landmask(grd_filename);
    return landmask_array;


def read_landmask(netcdf_name):
    """Read the pixel node grd file created by gmt grdlandmask"""
    print("Reading landmask %s " % netcdf_name);
    ds = xr.open_dataset(netcdf_name);
    landmask = ds["z"];
    return landmask;
