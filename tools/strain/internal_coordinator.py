"""
Driver program for strain calculation
Options:
1. Delaunay (on sphere, or delaunay_flat)
2. gpsgridder
3. visr
4. tape
5. huang
"""

from strain import (
    input_manager, 
    strain_delaunay, 
    strain_delaunay_flatearth, 
    strain_gpsgridder,
    strain_visr,
    strain_huang,
    output_manager,
)


compute_dict = {
        "delaunay": strain_delaunay.compute,
        "delaunay_flat": strain_delaunay_flatearth.compute,
        "gps_gridder": strain_gpsgridder.compute,
        "visr": strain_visr.compute,
        "huang": strain_huang.compute 
    };


def strain_coordinator(MyParams):
    Inputs = input_manager.inputs(MyParams);
    strain_fun = compute_dict[MyParams.strain_method]
    [lons, lats, rot, exx, exy, eyy] = strain_fun(Inputs, MyParams);
    output_manager.outputs_2d(
            lons, 
            lats, 
            rot, 
            exx, 
            exy, 
            eyy, 
            MyParams, 
            Inputs
        );  # constant 2D output format
    return
