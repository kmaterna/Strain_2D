import unittest
import strain_tensor_toolbox
import configure_functions


class Tests(unittest.TestCase):

    def test_solid_body_rotation(self):
        # This example is solid body rotation.
        # Does solid body rotation come out?
        up, vp = 0, -1;
        ur, vr = -1, -1;
        uq, vq = 0, 0;
        xinc, yinc = 1, 1;
        [dudx, dvdx, dudy, dvdy] = strain_tensor_toolbox.compute_displacement_gradients(up, vp, ur, vr, uq, vq, xinc,
                                                                                        yinc);
        [exx, exy, eyy, rotation] = strain_tensor_toolbox.compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy);
        print("Testing solid body rotation.")
        print("exx, exy, eyy: %f %f %f" % (exx, exy, eyy));
        print("rotation: %f" % rotation);
        self.assertEqual(exx, 0);
        self.assertEqual(exy, 0);
        self.assertEqual(eyy, 0);
        self.assertEqual(rotation, 1000);
        return;

    def test_reading_config(self):
        MyParams = configure_functions.config_parser(configfile="test/example_config.txt");
        self.assertTrue(MyParams);
        return;


if __name__ == "__main__":
    unittest.main();
