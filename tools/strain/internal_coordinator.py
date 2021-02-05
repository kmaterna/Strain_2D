"""
Driver program for strain calculation
Options:
1. Delaunay (on sphere, or delaunay_flat)
2. gpsgridder
3. visr
4. tape
5. huang
"""
import importlib

from strain import input_manager, output_manager


def model(model_name):
    ''' Return an instance of the model model_name '''
    module_name = 'strain.models.strain_' + model_name.lower()
    model_module = importlib.import_module(module_name)
    obj = getattr(model_module, model_name.delaunay().replace('-', ''))
    return module_name, obj


def strain_coordinator(MyParams):
    velField = input_manager.inputs(MyParams);
    model_name, strain_fun = model(MyParams.strain_method)
    [lons, lats, rot, exx, exy, eyy] = strain_fun(velField, MyParams);
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
