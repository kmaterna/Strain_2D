# Methods for turning a 1D computation into a 2D computation

import produce_gridded as pg
import strain_tensor_toolbox
from Tectonic_Utils.read_write import netcdf_read_write


def drive_tape(myParams):
    # generic steps
    map_range = [float(i) for i in myParams.map_range.split("/")]
    coord_box = myParams.coord_box;
    indir = "../compearth/surfacevel2strain/matlab_output/"
    print("Producing gridded dataset of: ")
    x, y, tt, tp, pp = pg.input_tape(indir, "cascadia_d02_q03_q06_b1_2D_s1_u1_strain.dat",
                                     "cascadia_d02_q03_q06_b1_2D_s1_u1_Dtensor_6entries.dat");
    [I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = pg.compute_tape(tt, tp, pp);
    # second invariant
    newx, newy, newI2nd = pg.nn_interp(x, y, I2nd, coord_box[0], coord_box[1], coord_box[2], coord_box[3],
                                       myParams.grid_inc);
    pg.output_tape(newx, newy, newI2nd, myParams.outdir, "I2nd.nc");
    # dilatation
    newx, newy, newdila = pg.nn_interp(x, y, dilatation, coord_box[0], coord_box[1], coord_box[2], coord_box[3],
                                       myParams.grid_inc);
    pg.output_tape(newx, newy, newdila, myParams.outdir, "dila.nc");
    # max shear
    newx, newy, newmax = pg.nn_interp(x, y, max_shear, coord_box[0], coord_box[1], coord_box[2], coord_box[3],
                                      myParams.grid_inc);
    pg.output_tape(newx, newy, newmax, myParams.outdir, "max_shear.nc");
    # azimuth
    azimuth = strain_tensor_toolbox.max_shortening_azimuth_1d(e1, e2, v00, v01, v10, v11)
    newx, newy, newaz = pg.nn_interp(x, y, azimuth, coord_box[0], coord_box[1], coord_box[2], coord_box[3],
                                     myParams.grid_inc);
    netcdf_read_write.produce_output_netcdf(newx, newy, newaz, 'degrees', myParams.outdir+'azimuth.nc');
    pg.write_tape_eigenvectors(x, y, e1, e2, v00, v01, v10, v11)
    return;
