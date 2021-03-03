import subprocess, sys, os
import collections
import configparser
from . import utilities

Params = collections.namedtuple("Params", ['strain_method', 'input_file', 'range_strain', 'range_data',
                                           'inc', 'outdir', 'method_specific']);
Comps_Params = collections.namedtuple("Comps_Params", ['range_strain', 'inc', 'strain_dict', 'outdir']);

help_message = "  Welcome to a geodetic strain calculator.\n" \
               "  USAGE: strain_driver.py config.txt\n" \
               "  See repository source for an example config file.\n"
comps_help_message = "  Welcome to a geodetic strain-rate comparison tool.\n" \
                     "  USAGE: compare_driver config.txt\n" \
                     "  See repository source for an example config file.\n"


def strain_config_parser(cmdargs=None, configfile=None):
    # The configfile can be passed as argv (if bash API) or passed as argument (if python API)
    if not configfile:
        if len(cmdargs) < 2:
            print(help_message);
            sys.exit(0);
        else:
            if cmdargs[1] == '--help':
                print(help_message);
                sys.exit(0);
            configfile = cmdargs[1];

    MyParams = parse_config_file_into_Params(configfile);
    subprocess.call(['mkdir', '-p', MyParams.outdir], shell=False);
    subprocess.call(['cp', configfile, MyParams.outdir], shell=False);

    print("\n------------------------------");
    print("Hello! We are...");
    print("   Computing strain using : %s " % MyParams.strain_method);
    print("   Input data from        : %s" % MyParams.input_file);
    print("   Calculation range      : %s" % MyParams.range_strain);
    print("   Putting the outputs    : %s \n" % MyParams.outdir);
    return MyParams;


def comparison_config_parser(args=None, configfile=None):
    # The configfile can be passed as argv (if bash API) or passed as argument (if python API)
    if not configfile:
        if len(args) < 2:
            print(comps_help_message);
            sys.exit(0);
        else:
            if args[1] == '--help':
                print(comps_help_message);
                sys.exit(0);
            configfile = args[1];
    MyParams = parse_comparison_config_into_Params(configfile);
    subprocess.call(['mkdir', '-p', MyParams.outdir], shell=False);
    subprocess.call(['cp', configfile, MyParams.outdir], shell=False);
    print("\n------------------------------");
    print("Hello! We are comparing strain calculations...");
    print("   Comparing strain from calculations : \n  %s " % MyParams.strain_dict);
    print("   Putting the outputs    : %s \n" % MyParams.outdir);
    return MyParams;


def parse_config_file_into_Params(configfile):
    """ Dedicated function to building a valid Params structure from the configfile """
    if not os.path.isfile(configfile):
        print("config file =  %s" % configfile);
        raise Exception("Error! config file was not found.");

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
    return MyParams;


def parse_comparison_config_into_Params(configfile):
    # Dedicated file to building a valid Params structure from the comps configfile
    if not os.path.isfile(configfile):
        print("config file =  %s" % configfile);
        raise Exception("Error! config file was not found.");
    config = configparser.ConfigParser();
    config.read(configfile);
    output_dir = config.get('general', 'output_dir');
    range_strain = config.get('strain', 'range_strain');
    range_strain = utilities.get_float_range(range_strain);
    inc = config.get('strain', 'inc');
    inc = utilities.get_float_inc(inc);
    # Reading the methods and their associated strain data
    specific_keys = [item for item in config["inputs"].keys()];
    strain_dict = {};
    for item in specific_keys:
        strain_dict[item] = config.get("inputs", item);
    MyParams = Comps_Params(inc=inc, range_strain=range_strain, outdir=output_dir, strain_dict=strain_dict);
    return MyParams;
