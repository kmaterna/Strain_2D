"""
Driver program for strain calculation
"""
import importlib
from . import input_manager, output_manager


def model(model_name):
    """ Return an instance of the model model_name.
        Dynamic module discovery. """
    module_name = 'Strain_2D.tools.strain.models.strain_' + model_name.lower()
    model_module = importlib.import_module(module_name)
    obj = getattr(model_module, model_name)
    return module_name, obj()


def strain_coordinator(MyParams):
    velField = input_manager.inputs(MyParams);
    model_name, strain_model = model(MyParams.strain_method);
    [lons, lats, rot, exx, exy, eyy] = strain_model.compute(velField, MyParams);
    output_manager.outputs_2d(lons, lats, rot, exx, exy, eyy, MyParams, velField);  # 2D grid output format
    return
