"""
Driver program for strain calculation
Options:
1. Delaunay (on sphere, or delaunay_flat)
2. gpsgridder
3. visr
4. tape
5. huang
"""

from Strain_2D.strain import input_manager
from Strain_2D.strain import strain_delaunay
from Strain_2D.strain import strain_delaunay_flatearth
from Strain_2D.strain import strain_gpsgridder
from Strain_2D.strain import strain_visr
from Strain_2D.strain import strain_huang
from Strain_2D.strain import output_manager


compute_dict = {
    "delaunay": strain_delaunay.compute,
    "delaunay_flat": strain_delaunay_flatearth.compute,
    "gps_gridder": strain_gpsgridder.compute,
    "visr": strain_visr.compute,
    "huang": strain_huang.compute };


def strain_coordinator(MyParams):
    Inputs = input_manager.inputs(MyParams);
    [lons, lats, rot, exx, exy, eyy] = compute_dict[MyParams.strain_method](Inputs, MyParams);
    output_manager.outputs_2d(lons, lats, rot, exx, exy, eyy, MyParams, Inputs);  # constant 2D output format
    return;
