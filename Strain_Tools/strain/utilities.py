# A set of utility functions used throughout the Strain_2D library
import numpy as np


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


# --------- GRID AND NAMED TUPLE UTILITIES ------------------ #

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
    """Read velocities from a NamedTuple"""
    lon, lat, e, n, se, sn = [], [], [], [], [], [];
    for item in velField:
        lon.append(item.elon)
        lat.append(item.nlat)
        e.append(item.e)
        n.append(item.n)
        se.append(item.se)
        sn.append(item.sn)
    return np.array(lon), np.array(lat), np.array(e), np.array(n), np.array(se), np.array(sn)
