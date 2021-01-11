import os
from . import compare_grd_functions as comp
from Tectonic_Utils.read_write import netcdf_read_write


def drive(MyParams, filename):
    # Some defensive programming (check all the same grid size)
    strain_values_dict = read_strain_files(MyParams, filename);
    comp.defensive_programming(MyParams, strain_values_dict);

    print("Comparing %s across all methods" % filename)

    # if component == "azimuth":
    #     my_means = comp.angle_means(lons2, lats2, val1, val2, val3, val4, val5)
    #     my_sds = comp.angle_sds(lons2, lats2, val1, val2, val3, val4, val5)
    # if component == "I2nd":
    #     my_means, my_sds = comp.grid_means_log(strain_values_dict)
    # if component == "dilatation":
    #     my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    # if component == "max_shear":
    #     my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    #
    # netcdf_read_write.produce_output_netcdf(lons2, lats2, my_means, 'per year', outdir+"means"+component_name+".nc");
    # netcdf_read_write.produce_output_netcdf(lons2, lats2, my_sds, 'per year', outdir+"deviations"+component_name+".nc");
    #
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
