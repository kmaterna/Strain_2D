import sys
import os
import collections
import configparser
import shutil
from . import utilities

Params = collections.namedtuple("Params", ['strain_method', 'input_file', 'range_strain', 'range_data', 'inc',
                                           'xdata', 'ydata', 'outdir', 'method_specific', 'write_metrics'])
Comps_Params = collections.namedtuple("Comps_Params", ['range_strain', 'inc', 'strain_dict', 'outdir'])

avail_modules = "  delaunay\n  delaunay_flat\n  geostats\n  gpsgridder\n  loc_avg_grad\n  wavelets\n  visr\n  simple_visr\n   velmap\n"
help_message = "  Welcome to a geodetic strain-rate calculator.\n\n" \
               "  USAGE 1: strain_rate_compute.py config.txt      <-- for running a strain calculation\n" \
               "  USAGE 2: strain_rate_compute.py --help          <-- for printing help message\n" \
               "  USAGE 3: strain_rate_compute.py --print_config  <-- for writing example config file\n"

comps_help_message = "  Welcome to a geodetic strain-rate comparison tool.\n\n" \
                     "  USAGE 1: strain_rate_comparison.py config.txt      <-- for running a strain comparison\n" \
                     "  USAGE 2: strain_rate_comparison.py --help        <-- for printing help message\n"


def strain_cmd_parser(cmdargs):
    """The configfile is passed as arg. This function will create the output directory if it does not exist. """
    if len(cmdargs) < 2:
        print(help_message)
        sys.exit(0)
    else:
        if cmdargs[1] == '--help':    # help message
            print(help_message)
            print("Available modules are:")
            print(avail_modules)
            sys.exit(0)
        elif cmdargs[1] == '--print_config':    # print example config file
            print("Writing example config file.")
            write_example_strain_config("example_strain_config.txt")
            sys.exit(0)
        else:                           # run the main program
            configfile = cmdargs[1]
            MyParams = read_strain_config(configfile)
            print_Params(MyParams)
            os.makedirs(str(MyParams.outdir), exist_ok=True)
            shutil.copyfile(configfile, os.path.join(str(MyParams.outdir), configfile))
            return MyParams


def read_strain_config(configfile, desired_method=None):
    """
    Build a Params structure from the configfile and set up the output directory.
    Will only read the configuration parameters of the chosen method.
    An optional method override is provided for testing purposes.

    :param configfile: string, filename
    :param desired_method: optional, string, can override the chosen method in the config file to read other values
    :returns: a Params named tuple
    """
    assert os.path.isfile(configfile), FileNotFoundError("Error! config file "+configfile+" was not found.")

    config = configparser.ConfigParser()
    config.read(configfile)
    if desired_method is None:
        strain_method = config.get('general', 'method')
    else:
        strain_method = desired_method
    output_dir = config.get('general', 'output_dir')
    input_file = config.get('general', 'input_vel_file')
    write_metrics = config.getint('general', 'write_metrics') if config.has_option('general', 'write_metrics') else 0
    range_strain = config.get('strain', 'range_strain')
    range_data = config.get('strain', 'range_data') if config.has_option('strain', 'range_data') else range_strain
    inc = config.get('strain', 'inc')
    if range_data == '':
        range_data = range_strain

    # Reading the method-specific stuff
    specific_keys = [item for item in config[strain_method].keys()]
    method_specific = {}
    for item in specific_keys:
        method_specific[item] = config.get(strain_method, item)

    # Cleanup and Grid Specification
    output_dir = os.path.join(output_dir, strain_method, '')
    range_strain = utilities.get_float_range(range_strain)
    range_data = utilities.get_float_range(range_data)
    inc = utilities.get_float_inc(inc)
    xdata, ydata, _ = utilities.make_grid(range_strain, inc)
    MyParams = Params(strain_method=strain_method, input_file=input_file, range_strain=range_strain,
                      range_data=range_data, inc=inc, xdata=xdata, ydata=ydata,
                      outdir=output_dir, method_specific=method_specific, write_metrics=write_metrics)
    return MyParams


def print_Params(MyParams):
    print("\n------------------------------")
    print("Hello! We are...")
    print("   Computing strain using : %s " % MyParams.strain_method)
    print("   Input data from        : %s" % MyParams.input_file)
    print("   Calculation range      : %s" % MyParams.range_strain)
    print("   Putting the outputs    : %s \n" % MyParams.outdir)
    return


def write_example_strain_config(outfile):
    """Write an example strain config file, useful for tutorial."""
    configobj = configparser.ConfigParser()
    configobj["general"] = {}
    configobj["strain"] = {}
    configobj["delaunay"] = {}
    configobj["delaunay_flat"] = {}
    configobj["visr"] = {}
    configobj["simple_visr"] = {}
    configobj["gpsgridder"] = {}
    configobj["loc_avg_grad"] = {}
    configobj["wavelets"] = {}
    configobj["geostats"] = {}
    configobj["strain-comparison"] = {}
    configobj["velmap"] = {}
    genconfig = configobj["general"]
    genconfig["method"] = "delaunay"
    genconfig["output_dir"] = "Output"
    genconfig["input_vel_file"] = "../test/testing_data/NorCal_stationvels.txt"
    genconfig["write_metrics"] = "0"
    strainconfig = configobj["strain"]
    strainconfig["range_strain"] = "-125/-120/38/42"
    strainconfig["range_data"] = "-125/-119/37.5/42.5"
    strainconfig["inc"] = "0.04/0.04"
    d1 = configobj["visr"]
    d1["distance_weighting"] = "gaussian"
    d1["spatial_weighting"] = "voronoi"
    d1["min_max_inc_smooth"] = "1/100/1"
    d1["weighting_threshold"] = "2"
    d1["uncertainty_threshold"] = "0.05"
    d1["num_creeping_faults"] = "0"
    d1["creep_file"] = "crp.dat"
    d1["executable"] = "../contrib/visr/visr.exe"
    d1s = configobj["simple_visr"]
    d1s["weighting_threshold"] = "2"
    d1s["distance_method"] = "gaussian"
    d1s["coverage_method"] = "voronoi"
    d2 = configobj["gpsgridder"]
    d2["poisson"] = "0.5"
    d2["fd"] = "0.01"
    d2["eigenvalue"] = "0.0005"
    d3 = configobj["loc_avg_grad"]
    d3["EstimateRadiusKm"] = "80"
    d3["nstations"] = "8"
    d4 = configobj["wavelets"]
    d4["code_dir"] = ""
    d4["qmin"] = "4"
    d4["qmax"] = "7"
    d4["qsec"] = "7"
    d4["output_tag"] = ""
    d5 = configobj["geostats"]
    d5["model_type"] = "Gaussian"
    d5["sill_east"] = "20"
    d5["range_east"] = "0.2613"  # This is 29 km / 111 km /deg
    d5["nugget_east"] = "3"
    d5["sill_north"] = "20"
    d5["range_north"] = "0.342"  # This is 38 km / 111 km / deg
    d5["nugget_north"] = "6"
    d5["trend"] = "0"
    d6 = configobj["velmap"]
    d6["smoothing_constant"] = "1e-2"
    dcomps = configobj["strain-comparison"]
    dcomps["output_dir"] = "Output/_strain_comparison"
    dcomps["input_dirs"] = "Output/delaunay:Output/gpsgridder:Output/loc_avg_grad:Output/visr"

    with open(outfile, 'w') as configfile:
        configobj.write(configfile)
    print("Writing file %s " % outfile)
    return


def comparison_cmd_parser(args):
    """The configfile is passed as arg"""
    if len(args) < 2:
        args = ("", "--help", "-h")   # help message
    if args[1] == '--help':
        print(comps_help_message)
        sys.exit(0)
    else:
        configfile = args[1]
        MyParams = read_comparison_config(configfile)
        os.makedirs(str(MyParams.outdir), exist_ok=True)
        shutil.copyfile(configfile, str(MyParams.outdir)+'/'+configfile)
        return MyParams


def read_comparison_config(configfile):
    """Build a comparison Params structure from the configfile"""
    assert (os.path.isfile(configfile)), FileNotFoundError("Error! config file " + configfile + " was not found.")
    config = configparser.ConfigParser()
    config.read(configfile)
    range_strain = config.get('strain', 'range_strain')
    range_strain = utilities.get_float_range(range_strain)
    inc = config.get('strain', 'inc')
    inc = utilities.get_float_inc(inc)
    # Reading the methods and their associated strain data
    input_dirs = config.get('strain-comparison', 'input_dirs')
    strain_dict = {}
    for item in input_dirs.split(':'):
        strain_method = item  # get the name of the technique
        strain_dict[strain_method] = item
    output_dir = config.get('strain-comparison', 'output_dir')
    MyParams = Comps_Params(inc=inc, range_strain=range_strain, outdir=output_dir, strain_dict=strain_dict)
    print("\n------------------------------")
    print("Hello! We are comparing strain calculations...")
    print("   Comparing strain from calculations : \n  %s " % strain_dict.values())
    print("   Putting the outputs    : %s \n" % MyParams.outdir)
    return MyParams
