"""
Driver program for strain calculation
Options:
1. Delaunay (on sphere)
2. gpsgridder
3. visr
4. tape
5. huang
"""

import input_manager
import strain_delaunay
import strain_gpsgridder
import strain_visr


compute_dict = {
    "delaunay": strain_delaunay.compute,
    "gpsgridder": strain_gpsgridder.compute,
    "visr": strain_visr.compute };


def strain_coordinator(MyParams):
    Inputs = input_manager.inputs(MyParams);
    compute_dict[MyParams.strain_method](Inputs, MyParams);
    return;
