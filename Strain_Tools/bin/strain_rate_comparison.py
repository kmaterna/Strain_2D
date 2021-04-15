#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Driver program for empirical strain statistics on strain computations.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import sys
from Strain_2D.Strain_Tools.strain import configure_functions
from Strain_2D.Strain_Tools.strain import compare_strain_grids


if __name__ == "__main__":
    MyParams = configure_functions.comparison_config_parser(args=sys.argv);
    compare_strain_grids.drive(MyParams);
