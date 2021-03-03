#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Driver program for strain calculation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import sys
from Strain_2D.Tools.strain import internal_coordinator
from Strain_2D.Tools.strain import configure_functions


if __name__ == "__main__":
    MyParams = configure_functions.strain_config_parser(cmdargs=sys.argv);
    internal_coordinator.strain_coordinator(MyParams);
