# A set of utility functions used throughout the Strain_2D library
import numpy as np


def get_float_range(string_range):
    """
    string range: format "-125/-121/32/35"
    float range: array of floats
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
    """Buffer is for the possible interface between pixel-node-registered and gridline-node-registered files"""
    string_range = str(float_range[0] - x_buffer) + '/' + str(float_range[1] + x_buffer) + '/' + \
                   str(float_range[2] - y_buffer) + '/' + str(float_range[3] + y_buffer);
    return string_range;


def get_float_inc(string_inc):
    """
    string_inc: like '0.04/0.04'
    float_inc: array of floats
    """
    number_incs = string_inc.split('/')
    float_inc = [float(number_incs[0]), float(number_incs[1])];
    return float_inc;


def get_string_inc(float_inc):
    string_inc = str(float_inc[0]) + '/' + str(float_inc[1]);
    return string_inc;


def mask_by_value(grid1, grid_maskingbasis, cutoff_value):
    """
    Will NAN mask one grid in all places where values are smaller than a cutoff value in a corresponding grid.
    grid1 = usually azimuth deviations
    grid2 = usually I2nd
    Returns: masked grid.
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

def check_coregistered_shapes(strain_values_dict):
    """Make sure arrays are of the same dimensions before attempting to produce any statistics"""
    method1 = list(strain_values_dict.keys())[0];
    for method2 in strain_values_dict.keys():
        v1 = strain_values_dict[method1][2];  # strain values array
        v2 = strain_values_dict[method2][2];  # second strain values array
        assert(np.shape(v1) == np.shape(v2)), ValueError("Error! Not all arrays have the same shape!  Cannot compare.");
    print("All methods have the same shape!");
    return;


def check_coregistered_grids(range_strain, inc, strain_values_dict):
    """
    Check a number of grids for ranges and increment consistent with the parameter ranges and increments
    Within a certain number of decimal places
    """
    range_strain = np.round(range_strain, 6);
    inc = np.round(inc, 6);
    for method in strain_values_dict.keys():
        assert (range_strain[0] == np.round(strain_values_dict[method][0][0], 6)), ValueError(
            "Lon of " + method + " doesn't match specs");
        assert (range_strain[1] == np.round(strain_values_dict[method][0][-1], 6)), ValueError(
            "Lon of " + method + " doesn't match specs");
        assert (range_strain[2] == np.round(strain_values_dict[method][1][0], 6)), ValueError(
            "Lat of " + method + " doesn't match specs");
        assert (range_strain[3] == np.round(strain_values_dict[method][1][-1], 6)), ValueError(
            "Lat of " + method + " doesn't match specs");
        xinc = np.round(strain_values_dict[method][0][1] - strain_values_dict[method][0][0], 6);
        yinc = np.round(strain_values_dict[method][1][1] - strain_values_dict[method][1][0], 6);
        assert (xinc == inc[0]), ValueError("xinc of " + method + " doesn't match specs");
        assert (yinc == inc[1]), ValueError("yinc of " + method + " doesn't match specs");
    return;
