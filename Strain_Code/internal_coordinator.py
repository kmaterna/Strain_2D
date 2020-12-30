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
import strain_delaunay_flatearth
import strain_gpsgridder
import strain_visr
import strain_tensor_toolbox
import output_manager


compute_dict = {
    "delaunay": strain_delaunay.compute,
    "delaunay_flat": strain_delaunay_flatearth.compute,
    "gpsgridder": strain_gpsgridder.compute,
    "visr": strain_visr.compute };


def strain_coordinator(MyParams):
    Inputs = input_manager.inputs(MyParams);
    # Here we will make a constant output format with 2d grids
    [xcentroid, ycentroid, triangle_vertices, rot, e1, e2, v00, v01, v10, v11] = compute_dict[MyParams.strain_method](Inputs, MyParams);
    # For 2D grid outputs, eventually
    # [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(e1, e2, v00, v01, v10, v11);
    # output_manager.outputs_1d(xcentroid, ycentroid, triangle_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation, azimuth, Inputs, MyParams);
    return;
