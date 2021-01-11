"""
Driver program for empirical strain statistics.
Configure which grids are being compared in the config file passed into this library.
"""

import sys
from Strain_2D.strain import configure_functions
from Strain_2D.strain import compare_strain_grids


if __name__ == "__main__":
    MyParams = configure_functions.comparison_config_parser(args=sys.argv);
    compare_strain_grids.drive(MyParams, "I2nd.nc");
    # compare_strain_grids.drive("azimuth");
    # compare_strain_grids.drive("dilatation");
    # compare_strain_grids.drive("max_shear");
