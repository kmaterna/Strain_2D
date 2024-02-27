#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Driver program for computing moment via Savage and Simpson, 1997
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import sys
from strain import moment_functions

if __name__ == "__main__":
    MyParams = moment_functions.cmd_parser(cmdargs=sys.argv)
    moment_functions.moment_coordinator(MyParams)
