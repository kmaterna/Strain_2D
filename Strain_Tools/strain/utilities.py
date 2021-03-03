# A set of utility functions used throughout the Strain_2D library
import numpy as np


def get_float_range(string_range):
    # string range: format "-125/-121/32/35"
    # float range: array of floats
    number_strings = string_range.split('/')
    float_range = [float(number_strings[0]), float(number_strings[1]),
                   float(number_strings[2]), float(number_strings[3])];
    if float_range[1] <= float_range[0]:
        raise ValueError("Error! Given range is invalid", float_range);
    if float_range[3] <= float_range[2]:
        raise ValueError("Error! Given range is invalid", float_range);
    return float_range;


def get_string_range(float_range, x_buffer=0, y_buffer=0):
    # Buffer is for the possible interface between pixel-node-registered and gridline-node-registered files
    string_range = str(float_range[0]-x_buffer)+'/'+str(float_range[1]+x_buffer)+'/' +\
                   str(float_range[2]-y_buffer)+'/'+str(float_range[3]+y_buffer);
    return string_range;


def get_float_inc(string_inc):
    # string_inc: like '0.04/0.04'
    # float_inc: array of floats
    number_incs = string_inc.split('/')
    float_inc = [float(number_incs[0]), float(number_incs[1])];
    return float_inc;


def get_string_inc(float_inc):
    string_inc = str(float_inc[0])+'/'+str(float_inc[1]);
    return string_inc;


def mask_by_value(grid1, grid_maskingbasis, cutoff_value):
    # Will NAN mask one grid in all places where values are smaller than a cutoff value in a corresponding grid.
    # grid1 = usually azimuth deviations
    # grid2 = usually I2nd
    # Returns: masked grid.
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
    for method in strain_values_dict.keys():
        print("%s range: %.2f %.2f %.2f %.2f " % (method, np.min(strain_values_dict[method][0]),
                                                  np.max(strain_values_dict[method][0]),
                                                  np.min(strain_values_dict[method][1]),
                                                  np.max(strain_values_dict[method][1])) );
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
