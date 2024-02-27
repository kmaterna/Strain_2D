#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Driver program for empirical strain statistics on strain computations.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import sys
from strain import configure_functions, compare_strain_grids

if __name__ == "__main__":
    MyParams = configure_functions.comparison_cmd_parser(args=sys.argv)
    compare_strain_grids.drive(MyParams)
