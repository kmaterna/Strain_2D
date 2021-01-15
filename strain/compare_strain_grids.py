import os
from . import compare_grd_functions as comp
from Tectonic_Utils.read_write import netcdf_read_write


def drive(MyParams):
    print("Comparing across all strain methods");
    compare_second_invariants(MyParams, "I2nd.nc");
    compare_max_shear(MyParams, "max_shear.nc");
    compare_dilatation(MyParams, "dila.nc");
    compare_rotation(MyParams, "rot.nc");
    compare_azimuth(MyParams, "azimuth.nc");
    return;


def compare_second_invariants(MyParams, filename):
    strain_values_dict, lons, lats = read_strain_files(MyParams, filename);
    my_means, my_stds = comp.grid_means_log(strain_values_dict)
    write_means_stds(lons, lats, my_means, my_stds, MyParams.outdir, filename);
    return;


def compare_max_shear(MyParams, filename):
    strain_values_dict, lons, lats = read_strain_files(MyParams, filename);
    my_means, my_stds = comp.grid_means_stds(strain_values_dict)
    write_means_stds(lons, lats, my_means, my_stds, MyParams.outdir, filename);
    return;


def compare_dilatation(MyParams, filename):
    strain_values_dict, lons, lats = read_strain_files(MyParams, filename);
    my_means, my_stds = comp.grid_means_stds(strain_values_dict);
    #     comp.mask_by_value(outdir, "dila", "dila", 15);
    write_means_stds(lons, lats, my_means, my_stds, MyParams.outdir, filename);
    return;


def compare_rotation(MyParams, filename):
    strain_values_dict, lons, lats = read_strain_files(MyParams, filename);
    my_means, my_stds = comp.grid_means_stds(strain_values_dict);
    write_means_stds(lons, lats, my_means, my_stds, MyParams.outdir, filename);
    return;


def compare_azimuth(MyParams, filename):
    strain_values_dict, lons, lats = read_strain_files(MyParams, filename);
    my_means, my_stds = comp.angle_means(strain_values_dict)
    write_means_stds(lons, lats, my_means, my_stds, MyParams.outdir, filename);
    comp.mask_by_value(MyParams.outdir, "azimuth", "I2nd", 3);
    return;


def read_strain_files(MyParams, filename):
    # Read strain quantities and guarantee their co-registration
    strain_values_dict = {};
    for method in MyParams.strain_dict.keys():
        specific_filename = MyParams.strain_dict[method]+"/"+filename
        if os.path.isfile(specific_filename):
            [lon, lat, val] = netcdf_read_write.read_any_grd(specific_filename);
            strain_values_dict[method] = [lon, lat, val];
        else:
            raise("Error! Can't find file %s " % specific_filename);
    comp.defensive_programming(MyParams, strain_values_dict);
    method1 = list(strain_values_dict.keys())[0];
    lons = strain_values_dict[method1][0];
    lats = strain_values_dict[method1][1];
    return strain_values_dict, lons, lats;


def write_means_stds(lons, lats, means, stds, outdir, filename):
    netcdf_read_write.produce_output_netcdf(lons, lats, means, 'per year', outdir+"/means_"+filename);
    netcdf_read_write.produce_output_netcdf(lons, lats, stds, 'per year', outdir+"/deviations_"+filename);
    return;
