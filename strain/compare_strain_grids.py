import os
from . import compare_grd_functions as comp
from Tectonic_Utils.read_write import netcdf_read_write


def drive(MyParams, filename):
    print("Comparing %s across all methods" % filename)
    # Some defensive programming (check all the same grid size)
    strain_values_dict = read_strain_files(MyParams, filename);
    comp.defensive_programming(MyParams, strain_values_dict);

    method1 = list(strain_values_dict.keys())[0];
    lons = strain_values_dict[method1][0];
    lats = strain_values_dict[method1][1];

    if filename == "azimuth.nc":
        my_means, my_stds = comp.angle_means(strain_values_dict)
    elif filename == "I2nd.nc":
        my_means, my_stds = comp.grid_means_log(strain_values_dict)
    elif filename == "dila.nc":
        my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    elif filename == "rot.nc":
        my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    elif filename == "max_shear.nc":
        my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    else:
        raise Exception("Cannot find filename %s " % filename);

    netcdf_read_write.produce_output_netcdf(lons, lats, my_means, 'per year', MyParams.outdir+"/means_"+filename);
    netcdf_read_write.produce_output_netcdf(lons, lats, my_stds, 'per year', MyParams.outdir+"/deviations_"+filename);

    # if component == "azimuth":
    #     # A special code to mask out values of azimuth where the magnitude of strain is really low.
    #     comp.mask_by_value(outdir, "azimuth", "I2nd", 3);
    # if component == "dilatation":
    #     # A special code to mask out values of azimuth where the magnitude of strain is really low.
    #     comp.mask_by_value(outdir, "dila", "dila", 15);
    return;


def read_strain_files(MyParams, filename):
    strain_values_dict = {};
    for method in MyParams.strain_dict.keys():
        specific_filename = MyParams.strain_dict[method]+"/"+filename
        if os.path.isfile(specific_filename):
            [lon, lat, val] = netcdf_read_write.read_any_grd(specific_filename);
            strain_values_dict[method] = [lon, lat, val];
        else:
            raise("Error! Can't find file %s " % specific_filename);
    return strain_values_dict;
