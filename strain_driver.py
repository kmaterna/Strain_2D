"""
Driver program for strain calculation
Options:
1. Delaunay (on sphere, or delaunay_flat)
2. gpsgridder
3. visr
4. tape
5. huang
"""


import sys
from Strain_2D.strain import configure_functions
from Strain_2D.strain import internal_coordinator


if __name__ == "__main__":
    MyParams = configure_functions.strain_config_parser(args=sys.argv);
    internal_coordinator.strain_coordinator(MyParams);
