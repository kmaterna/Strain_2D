import subprocess, sys, os
import collections
import configparser

Params = collections.namedtuple("Params", ['strain_method',
                                           'input_file', 'range_strain', 'range_data',
                                           'num_years', 'max_sigma', 'inc', 'outdir', 'blacklist_file',
                                           'method_specific']);
available_methods = ['delaunay',
                     'delaunay_flat',
                     'visr',
                     'gps_gridder',
                     'tape',
                     'huang'];

help_message = "  Welcome to a geodetic strain calculator.\n" \
               "  USAGE: strain_driver config.txt\n" \
               "  See repository source for an example config file.\n"


def config_parser(args=None, configfile=None):
    # The configfile can be passed as argv (if bash API) or passed as argument (if python API)
    if not configfile:
        if len(args) < 2:
            print(help_message);
            sys.exit(0);
        else:
            configfile = args[1];

    MyParams = parse_config_file_into_Params(configfile);
    make_output_dir(MyParams.outdir);
    subprocess.call(['cp', configfile, MyParams.outdir], shell=False);

    print("\n------------------------------");
    print("Hello! We are...");
    print("   Computing strain using : %s " % MyParams.strain_method);
    print("   Input data from        : %s" % MyParams.input_file);
    print("   Calculation range      : %s" % MyParams.range_strain);
    print("   Putting the outputs    : %s \n" % MyParams.outdir);
    return MyParams;


def parse_config_file_into_Params(configfile):
    # Dedicated file to building a valid Params structure from the configfile
    if not os.path.isfile(configfile):
        print("config file =  %s" % configfile);
        raise Exception("Error! config file was not found.");

    config = configparser.ConfigParser()
    config.read(configfile)
    strain_method = config.get('general', 'method');
    output_dir = config.get('general', 'output_dir');
    input_file = config.get('inputs', 'vel_file');
    blacklist_file = config.get('inputs', 'blacklist') if config.has_option('inputs', 'blacklist') else '';
    num_years = config.getfloat('inputs', 'num_years') if config.has_option('inputs', 'num_years') else 0;
    max_sigma = config.getfloat('inputs', 'max_sigma') if config.has_option('inputs', 'max_sigma') else 100;
    range_strain = config.get('strain', 'range_strain');
    range_data = config.get('strain', 'range_data') if config.has_option('strain', 'range_data') else range_strain;
    inc = config.get('strain', 'inc');
    if range_data == '':
        range_data = range_strain;
    method_specific = {};   # will write later

    # Cleanup
    output_dir = output_dir + '/' + strain_method + '/'
    range_strain = get_float_range(range_strain);
    range_data = get_float_range(range_data);
    inc = get_float_inc(inc);
    MyParams = Params(strain_method=strain_method, input_file=input_file, range_strain=range_strain,
                      range_data=range_data, num_years=num_years, max_sigma=max_sigma, inc=inc, outdir=output_dir,
                      blacklist_file=blacklist_file, method_specific=method_specific);
    sanity_check_inputs(MyParams)
    return MyParams;


def make_output_dir(outer_name):
    subprocess.call(['mkdir', '-p', outer_name], shell=False);
    return;


def get_float_range(string_range):
    # string range: format "-125/-121/32/35"
    # float range: array of floats
    number_strings = string_range.split('/')
    float_range = [float(number_strings[0]), float(number_strings[1]),
                   float(number_strings[2]), float(number_strings[3])];
    return float_range;


def get_float_inc(string_inc):
    # string_inc: like '0.04/0.04'
    # float_inc: array of floats
    number_incs = string_inc.split('/')
    float_inc = [float(number_incs[0]), float(number_incs[1])];
    return float_inc;


def sanity_check_inputs(MyParams):
    # For options that change based on strain method,
    # Check that the right ones exist.
    # Specific logic here.
    if MyParams.strain_method not in available_methods:
        raise Exception("%s is not a known strain method. Exiting.\n" % MyParams.strain_method);
    if MyParams.strain_method == "gps_gridder":
        if 'distance' not in MyParams.method_specific.keys():
            raise Exception("\ngps_gridder requires distance buffer. Please add to method_specific config. Exiting.\n");
    elif MyParams.strain_method == "visr":
        if 'distance' not in MyParams.method_specific.keys():
            raise Exception("\nvisr requires distance buffer. Please add to method_specific config. Exiting.\n");
    return;
