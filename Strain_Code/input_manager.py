# The input manager for GPS Strain analysis. 

import gps_vel_functions
import gps_io_functions


# ----------------- INPUTS -------------------------
def inputs(MyParams):
    print("------------------------------");
    # Purpose: generate input velocity field.
    if 'PBO' in MyParams.input_file or 'pbo' in MyParams.input_file:
        [myVelfield] = gps_io_functions.read_pbo_vel_file(MyParams.input_file);
    elif 'MAGNET' in MyParams.input_file or 'unr' in MyParams.input_file or 'midas' in MyParams.input_file:
        raise Exception("MAGNET files not yet supported.");
        # [myVelfield] = gps_io_functions.read_unr_vel_file(MyParams.input_file, 'coord_file.txt');  # COME BACK LATER
    elif 'SCEC' in MyParams.input_file:
        [myVelfield] = gps_io_functions.read_gamit_velfile(MyParams.input_file);
    else:
        raise Exception("Error! Cannot read %s " % MyParams.input_file);

    print("%d stations before applying cleaning." % (len(myVelfield)));
    blacklist = gps_io_functions.read_blacklist(MyParams.blacklist_file);
    myVelfield = gps_vel_functions.remove_blacklist_vels(myVelfield, blacklist);
    myVelfield = gps_vel_functions.clean_velfield(myVelfield, num_years=MyParams.num_years,
                                                  max_horiz_sigma=MyParams.max_sigma, coord_box=MyParams.range_data);
    myVelfield = gps_vel_functions.remove_duplicates(myVelfield);
    print("%d stations after selection criteria.\n" % (len(myVelfield)));
    return myVelfield;
