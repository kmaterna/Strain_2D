import subprocess, sys, os
import collections
import configparser
from . import utilities

Params = collections.namedtuple("Params", ['strain_method', 'input_file', 'range_strain', 'range_data',
                                           'inc', 'outdir', 'method_specific']);
Comps_Params = collections.namedtuple("Comps_Params", ['range_strain', 'inc', 'strain_dict', 'outdir']);

avail_modules = "  delaunay\n  delaunay_flat\n  geostats\n  gpsgridder\n  huang\n  tape\n  visr\n "
help_message = "  Welcome to a geodetic strain-rate calculator.\n\n" \
               "  USAGE 1: strain_rate_compute.py config.txt      <-- for running a strain calculation\n" \
               "  USAGE 2: strain_rate_compute.py --help          <-- for printing help message\n" \
               "  USAGE 3: strain_rate_compute.py --print_config  <-- for writing example config file\n"

comps_help_message = "  Welcome to a geodetic strain-rate comparison tool.\n\n" \
                     "  USAGE 1: strain_rate_comparison.py config.txt      <-- for running a strain comparison\n" \
                     "  USAGE 2: strain_rate_comparison.py --help        <-- for printing help message\n"


def strain_cmd_parser(cmdargs):
    """The configfile is passed as arg"""
    if len(cmdargs) < 2:
        print(help_message);
        sys.exit(0);
    else:
        if cmdargs[1] == '--help':    # help message
            print(help_message);
            print("Available modules are:")
            print(avail_modules);
            sys.exit(0);
        elif cmdargs[1] == '--print_config':    # print example config file
            print("Writing example config file.");
            write_example_strain_config("example_strain_config.txt");
            sys.exit(0);
        else:                           # run the main program
            configfile = cmdargs[1];
            MyParams = read_strain_config(configfile);
            subprocess.call(['mkdir', '-p', MyParams.outdir], shell=False);
            subprocess.call(['cp', configfile, MyParams.outdir], shell=False);
            return MyParams;


def read_strain_config(configfile):
    """Build a Params structure from the configfile and set up the output directory"""
    assert(os.path.isfile(configfile)), FileNotFoundError("Error! config file "+configfile+" was not found.");

    config = configparser.ConfigParser()
    config.read(configfile)
    strain_method = config.get('general', 'method');
    output_dir = config.get('general', 'output_dir');
    input_file = config.get('general', 'input_vel_file');
    range_strain = config.get('strain', 'range_strain');
    range_data = config.get('strain', 'range_data') if config.has_option('strain', 'range_data') else range_strain;
    inc = config.get('strain', 'inc');
    if range_data == '':
        range_data = range_strain;

    # Reading the method-specific stuff
    specific_keys = [item for item in config[strain_method].keys()];
    method_specific = {};
    for item in specific_keys:
        method_specific[item] = config.get(strain_method, item);

    # Cleanup
    output_dir = output_dir + '/' + strain_method + '/'
    range_strain = utilities.get_float_range(range_strain);
    range_data = utilities.get_float_range(range_data);
    inc = utilities.get_float_inc(inc);
    MyParams = Params(strain_method=strain_method, input_file=input_file, range_strain=range_strain,
                      range_data=range_data, inc=inc, outdir=output_dir, method_specific=method_specific);

    print("\n------------------------------");
    print("Hello! We are...");
    print("   Computing strain using : %s " % MyParams.strain_method);
    print("   Input data from        : %s" % MyParams.input_file);
    print("   Calculation range      : %s" % MyParams.range_strain);
    print("   Putting the outputs    : %s \n" % MyParams.outdir);
    return MyParams;


def write_example_strain_config(outfile):
    """Write an example strain config file, useful for tutorial."""
    configobj = configparser.ConfigParser()
    configobj["general"] = {};
    configobj["strain"] = {};
    configobj["delaunay"] = {}
    configobj["delaunay_flat"] = {}
    configobj["visr"] = {}
    configobj["gpsgridder"] = {}
    configobj["huang"] = {}
    configobj["tape"] = {}
    configobj["geostats"] = {}
    configobj["strain-comparison"] = {}
    genconfig = configobj["general"];
    genconfig["method"] = "delaunay"
    genconfig["output_dir"] = "Output"
    genconfig["input_vel_file"] = "../test/testing_data/NorCal_stationvels.txt";
    strainconfig = configobj["strain"];
    strainconfig["range_strain"] = "-125/-120/38/42"
    strainconfig["range_data"] = "-125/-119/37.5/42.5"
    strainconfig["inc"] = "0.04/0.04"
    d1 = configobj["visr"];
    d1["distance_weighting"] = "gaussian";
    d1["spatial_weighting"] = "voronoi";
    d1["min_max_inc_smooth"] = "1/100/1";
    d1["executable"] = "../contrib/visr/visr.exe";
    d2 = configobj["gpsgridder"];
    d2["poisson"] = "0.5";
    d2["fd"] = "0.01";
    d2["eigenvalue"] = "0.0005";
    d3 = configobj["huang"];
    d3["EstimateRadiusKm"] = "80";
    d3["nstations"] = "8";
    d4 = configobj["tape"];
    d4["code_dir"] = "";
    d4["qmin"] = "4";
    d4["qmax"] = "7";
    d4["qsec"] = "7";
    d5 = configobj["geostats"];
    d5["model_type"] = "Gaussian";
    d5["sill_east"] = "30";
    d5["range_east"] = "1";
    d5["nugget_east"] = "1";
    d5["sill_north"] = "30";
    d5["range_north"] = "1";
    d5["nugget_north"] = "1";
    d5["trend"] = "0";
    dcomps = configobj["strain-comparison"];
    dcomps["output_dir"] = "Output/_strain_comparison"
    dcomps["input_dirs"] = "Output/delaunay:Output/gpsgridder:Output/huang:Output/visr"
    with open(outfile, 'w') as configfile:
        configobj.write(configfile)
    print("Writing file %s " % outfile);
    return;


def comparison_cmd_parser(args):
    """The configfile is passed as arg"""
    if len(args) < 2:
        args = ("", "--help");   # help message
    if args[1] == '--help':
        print(comps_help_message);
        sys.exit(0);
    else:
        configfile = args[1];
        MyParams = read_comparison_config(configfile);
        subprocess.call(['mkdir', '-p', MyParams.outdir], shell=False);
        subprocess.call(['cp', configfile, MyParams.outdir], shell=False);
        return MyParams;


def read_comparison_config(configfile):
    """Build a comparison Params structure from the configfile"""
    assert (os.path.isfile(configfile)), FileNotFoundError("Error! config file " + configfile + " was not found.");
    config = configparser.ConfigParser();
    config.read(configfile);
    range_strain = config.get('strain', 'range_strain');
    range_strain = utilities.get_float_range(range_strain);
    inc = config.get('strain', 'inc');
    inc = utilities.get_float_inc(inc);
    # Reading the methods and their associated strain data
    input_dirs = config.get('strain-comparison', 'input_dirs');
    strain_dict = {};
    for item in input_dirs.split(':'):
        strain_method = item.split('/')[-1]  # get the name of the technique
        strain_dict[strain_method] = item;
    output_dir = config.get('strain-comparison', 'output_dir');
    MyParams = Comps_Params(inc=inc, range_strain=range_strain, outdir=output_dir, strain_dict=strain_dict);
    print("\n------------------------------");
    print("Hello! We are comparing strain calculations...");
    print("   Comparing strain from calculations : \n  %s " % MyParams.strain_dict);
    print("   Putting the outputs    : %s \n" % MyParams.outdir);
    return MyParams;
