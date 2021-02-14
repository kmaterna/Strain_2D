import subprocess, sys, os
import collections
import configparser

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
    ''' Dedicated function to building a valid Params structure from the configfile '''
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
    range_strain = get_float_range(range_strain);
    range_data = get_float_range(range_data);
    inc = get_float_inc(inc);
    MyParams = Params(strain_method=strain_method, input_file=input_file, range_strain=range_strain,
                      range_data=range_data, inc=inc, outdir=output_dir, method_specific=method_specific);
    sanity_check_inputs(MyParams)
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
    range_strain = get_float_range(range_strain);
    inc = config.get('strain', 'inc');
    inc = get_float_inc(inc);
    # Reading the methods
    specific_keys = [item for item in config["inputs"].keys()];
    strain_dict = {};
    for item in specific_keys:
        strain_dict[item] = config.get("inputs", item);
    MyParams = Comps_Params(inc=inc, range_strain=range_strain, outdir=output_dir, strain_dict=strain_dict);
    return MyParams;


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


def sanity_check_inputs(MyParams):
    # For options that change based on strain method,
    # Check that the right ones exist.
    # Specific logic here.
    if MyParams.strain_method == "gps_gridder":
        if 'poisson' not in MyParams.method_specific.keys():
            raise ValueError("\ngps_gridder requires poisson's ratio. Please add to method_specific config. Exiting.\n");
        if 'fd' not in MyParams.method_specific.keys():
            raise ValueError("\ngps_gridder requires fudge factor fd. Please add to method_specific config. Exiting.\n");
        if 'eigenvalue' not in MyParams.method_specific.keys():
            raise ValueError("\ngps_gridder requires eigenvalue. Please add to method_specific config. Exiting.\n");
    elif MyParams.strain_method == "visr":
        if 'distance_weighting' not in MyParams.method_specific.keys():
            raise ValueError("\nvisr requires distance weighting. Please add to method_specific config. Exiting.\n");
        if 'spatial_weighting' not in MyParams.method_specific.keys():
            raise ValueError("\nvisr requires spatial weighting. Please add to method_specific config. Exiting.\n");
        if 'min_max_inc_smooth' not in MyParams.method_specific.keys():
            raise ValueError("\nvisr requires smoothing information. Please add to method_specific config. Exiting.\n");
        if 'executable' not in MyParams.method_specific.keys():
            raise ValueError("\nvisr requires path to executable. Please add to method_specific config. Exiting.\n");
    elif MyParams.strain_method == "huang":
        if 'estimateradiuskm' not in MyParams.method_specific.keys():
            raise ValueError("\nmethod requires estimateradiuskm. Please add to method_specific config. Exiting.\n");
        if 'nstations' not in MyParams.method_specific.keys():
            raise ValueError("\nmethod requires nstations. Please add to method_specific config. Exiting.\n");
    return;
