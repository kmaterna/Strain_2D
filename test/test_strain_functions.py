import unittest
import strain_tensor_toolbox
import configure_functions
import strain_delaunay_flatearth
import strain_delaunay
import gps_io_functions


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
        MyParams = configure_functions.parse_config_file_into_Params(configfile="test/testing_data/example_config.txt");
        self.assertTrue(MyParams);
        return;

    def test_delaunay_signs(self):
        station1 = gps_io_functions.Station_Vel(name="xxxx", elon=-123, nlat=39, e=0, n=0, u=0, se=1, sn=1, su=1,
                                                first_epoch=0, last_epoch=0, refframe='', proccenter='', subnetwork='',
                                                survey=0);
        station2 = gps_io_functions.Station_Vel(name="yyyy", elon=-124, nlat=39, e=-1, n=2, u=0, se=1, sn=1, su=1,
                                                first_epoch=0, last_epoch=0, refframe='', proccenter='', subnetwork='',
                                                survey=0);
        station3 = gps_io_functions.Station_Vel(name="zzzz", elon=-124, nlat=40, e=-3, n=4, u=0, se=1, sn=1, su=1,
                                                first_epoch=0, last_epoch=0, refframe='', proccenter='', subnetwork='',
                                                survey=0);
        station4 = gps_io_functions.Station_Vel(name="wwww", elon=-123, nlat=40, e=0, n=0, u=0, se=1, sn=1, su=1,
                                                first_epoch=0, last_epoch=0, refframe='', proccenter='', subnetwork='',
                                                survey=0);
        myVelfield = [station1, station2, station3, station4];
        [_, _, _, rot1, exx1, exy1, eyy1] = strain_delaunay.compute_with_delaunay_polygons(myVelfield);
        [_, _, _, rot2, exx2, exy2, eyy2] = strain_delaunay_flatearth.compute_with_delaunay_polygons(myVelfield);
        print("delaunay-sphere vs delaunay-flat:")
        print('exx:', exx1, ' vs ', exx2)
        print('exy:', exy1, ' vs ', exy2)
        print('eyy:', eyy1, ' vs ', eyy2);
        self.assertLess(abs(exx1[0]-exx1[0]), abs(exx1[0]*0.05));  # less than 5% difference
        self.assertLess(abs(exy1[0]-exy2[0]), abs(exy1[0]*0.05));  # less than 5% difference
        self.assertLess(abs(eyy1[0]-eyy2[0]), abs(eyy1[0]*0.05));  # less than 5% difference
        self.assertLess(abs(rot1[0]-rot2[0]), abs(rot1[0]*0.05));  # less than 5% difference
        return;


if __name__ == "__main__":
    unittest.main();
