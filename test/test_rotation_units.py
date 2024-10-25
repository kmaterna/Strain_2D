import unittest
from strain import velocity_io, configure_functions, output_manager
from strain.internal_coordinator import get_model
import os


class Tests(unittest.TestCase):

    def test_euler_pole_rotation(self):
        """
        Test a synthetic case in which a triangle of points is rotated about an Euler Pole by a known rate.
        The Euler Pole is ep = [-75, 65, 1]  # longitude, latitude, degrees/Ma
        Computed by: vels1 = euler_pole.point_rotation_by_Euler_Pole(point1, ep)
        The correct rotation rate is 17.45 radians per Ga.
        Test: Is that what is output?
        """
        print("Testing Euler Pole rotation for rotation units.")
        # Create the hard-coded setup to the test problem
        config_file = "test/testing_data/00_example_strain_config.txt"
        myvelfield = velocity_io.read_stationvels("test/testing_data/euler_pole_rotation_vels.txt")

        # Read the config parameters
        delaunay_params = configure_functions.read_strain_config(config_file, desired_method='delaunay')
        delaunay_flat_params = configure_functions.read_strain_config(config_file, desired_method='delaunay_flat')
        gpsgridder_params = configure_functions.read_strain_config(config_file, desired_method='gpsgridder')
        visr_params = configure_functions.read_strain_config(config_file, desired_method='visr')
        lag_params = configure_functions.read_strain_config(config_file, desired_method='loc_avg_grad')
        geostats_params = configure_functions.read_strain_config(config_file, desired_method='geostats')
        os.makedirs(name=os.path.join("test", "testing_data", "output"), exist_ok=True)  # make parent dir for outputs
        os.makedirs(str(delaunay_params.outdir), exist_ok=True)
        os.makedirs(str(delaunay_flat_params.outdir), exist_ok=True)
        os.makedirs(str(gpsgridder_params.outdir), exist_ok=True)
        os.makedirs(str(visr_params.outdir), exist_ok=True)
        os.makedirs(str(lag_params.outdir), exist_ok=True)
        os.makedirs(str(geostats_params.outdir), exist_ok=True)

        # Running several example strain methods on the simple dataset
        # Delaunay on Flat Earth
        strain_model = get_model("delaunay_flat")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(delaunay_flat_params)  # calling constructor, build strain model from params
        [Ve, Vn, rot_flat, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_flat, exx, exy, eyy, delaunay_flat_params, vels, resids)  # 2D grid outputs

        # Delaunay
        strain_model = get_model("delaunay")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(delaunay_params)  # calling constructor, building strain model from params
        [Ve, Vn, rot_del, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_del, exx, exy, eyy, delaunay_params, vels, resids)  # 2D grid output format

        # gpsgridder
        strain_model = get_model("gpsgridder")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(gpsgridder_params)  # calling constructor, building strain model from params
        [Ve, Vn, rot_gps, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_gps, exx, exy, eyy, gpsgridder_params, vels, resids)  # 2D grid output format

        # visr
        # Requires the fortran executable to be compiled on your architecture and findable.
        strain_model = get_model("visr")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(visr_params)  # calling constructor, building strain model from params
        [Ve, Vn, rot_visr, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_visr, exx, exy, eyy, visr_params, vels, resids)  # 2D grid output format

        # local_average_gradient
        strain_model = get_model("loc_avg_grad")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(lag_params)  # calling constructor, building strain model from params
        [Ve, Vn, rot_lag, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_lag, exx, exy, eyy, lag_params, vels, resids)  # 2D grid output format

        # geostats
        strain_model = get_model("geostats")  # getting an object of type that inherits from Strain_2d
        constructed_object = strain_model(geostats_params)  # calling constructor, building strain model from params
        [Ve, Vn, rot_geo, exx, exy, eyy, vels, resids] = constructed_object.compute(myvelfield)  # computing strain
        output_manager.outputs_2d(Ve, Vn, rot_geo, exx, exy, eyy, geostats_params, vels, resids)  # 2D grid output format

        print("DELAUNAY_FLAT:", rot_flat[3][3])
        print("DELAUNAY:", rot_del[3][3])
        print("GPSGRIDDER:", rot_gps[3][3])
        print("VISR:", rot_visr[3][3])
        print("LOC_AVG_GRAD:", rot_lag[3][3])
        print("GEOSTATS:", rot_geo[3][3])


if __name__ == "__main__":
    unittest.main()
