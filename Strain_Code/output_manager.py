# The output manager for GPS Strain analysis. 
# ----------------- OUTPUTS -------------------------

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
import strain_tensor_toolbox
import gps_io_functions


def outputs_2d(xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation, myVelfield, MyParams):
    print("Writing 2d outputs:");
    outfile = open(MyParams.outdir + "tempgps.txt", 'w');
    for i in range(len(myVelfield.n)):
        outfile.write("%f %f %f %f %f %f 0.0\n" % (
            myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i],
            myVelfield.sn[i]));
    outfile.close();
    azimuth = strain_tensor_toolbox.max_shortening_azimuth(e1, e2, v00, v01, v10, v11)
    netcdf_read_write.produce_output_netcdf(xdata, ydata, azimuth, 'degrees', MyParams.outdir + 'azimuth.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, I2nd, 'per yr', MyParams.outdir + 'I2nd.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, rot, 'per yr', MyParams.outdir + 'rot.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, dilatation, 'per yr', MyParams.outdir + 'dila.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, max_shear, 'per yr', MyParams.outdir + 'max_shear.nc');
    print("Max I2: %f " % (np.amax(I2nd)));
    print("Max rot: %f " % (np.amax(rot)));
    print("Min rot: %f " % (np.amin(rot)));
    write_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11, MyParams);
    gmt_file = open(MyParams.outdir + "run_gmt.gmt", 'w');
    gmt_file.write("../../../" + MyParams.gmtfile + " " + MyParams.map_range + "\n");
    gmt_file.close();
    upfile = open(MyParams.outdir + "uplift.txt", 'w');
    for i in range(len(myVelfield.n)):
        upfile.write("%f %f %f \n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.u[i]));
    upfile.close();
    print("../../../" + MyParams.gmtfile + " " + MyParams.map_range);
    return;


def write_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11, MyParams):
    # Need eigs_interval and outdir from MyParams.
    positive_file = open(MyParams.outdir + "positive_eigs.txt", 'w');
    negative_file = open(MyParams.outdir + "negative_eigs.txt", 'w');
    if MyParams.strain_method == 'visr':
        eigs_dec = 8;
    elif MyParams.strain_method == 'gpsgridder':
        eigs_dec = 12;
    elif MyParams.strain_method == 'spline':
        eigs_dec = 8;
    elif MyParams.strain_method == 'ND_interp':
        eigs_dec = 12;
    else:
        print("Error! strain method not recognized for eigenvector plotting.");

    do_not_print_value = 200;
    overmax_scale = 200;

    for j in range(len(ydata)):
        for k in range(len(xdata)):
            if np.mod(j, eigs_dec) == 0 and np.mod(k, eigs_dec) == 0:
                if w1[j][k] > 0:
                    scale = w1[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    positive_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k] * scale, v10[j][k] * scale));
                    positive_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k] * scale, -v10[j][k] * scale));
                if w1[j][k] < 0:
                    scale = w1[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    negative_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v00[j][k] * scale, v10[j][k] * scale));
                    negative_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v00[j][k] * scale, -v10[j][k] * scale));
                if w2[j][k] > 0:
                    scale = w2[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    positive_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k] * scale, v11[j][k] * scale));
                    positive_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k] * scale, -v11[j][k] * scale));
                if w2[j][k] < 0:
                    scale = w2[j][k];
                    if abs(scale) > do_not_print_value:
                        scale = overmax_scale;
                    negative_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], v01[j][k] * scale, v11[j][k] * scale));
                    negative_file.write(
                        "%s %s %s %s 0 0 0\n" % (xdata[k], ydata[j], -v01[j][k] * scale, -v11[j][k] * scale));
    positive_file.close();
    negative_file.close();

    return;


def outputs_1d(xcentroid, ycentroid, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation,
               azimuth, myVelfield, MyParams):
    print("------------------------------\nWriting 1d outputs:");
    write_multisegment_file(polygon_vertices, rot, MyParams.outdir+"rotation.txt");
    write_multisegment_file(polygon_vertices, I2nd, MyParams.outdir+"I2nd.txt");
    write_multisegment_file(polygon_vertices, dilatation, MyParams.outdir+"Dilatation.txt");
    write_multisegment_file(polygon_vertices, max_shear, MyParams.outdir+"max_shear.txt");
    write_multisegment_file(polygon_vertices, azimuth, MyParams.outdir+"azimuth.txt");
    gps_io_functions.write_humanread_vel_file(myVelfield, MyParams.outdir+"tempgps.txt");

    positive_file = open(MyParams.outdir + "positive_eigs.txt", 'w');
    negative_file = open(MyParams.outdir + "negative_eigs.txt", 'w');
    for i in range(len(I2nd)):
        # Write the eigenvectors and eigenvalues
        write_single_eigenvector(positive_file, negative_file, e1[i], v00[i], v10[i], xcentroid[i], ycentroid[i]);
        write_single_eigenvector(positive_file, negative_file, e2[i], v01[i], v11[i], xcentroid[i], ycentroid[i]);
    positive_file.close();
    negative_file.close();

    print("Max I2: %f " % (max(I2nd)));
    print("Max rot: %f " % (max(rot)));
    print("Min rot: %f " % (min(rot)));
    return;


def write_multisegment_file(polygon_vertices, quantity, filename):
    # Write a quantity for each polygon, in GMT-readable format
    ofile = open(filename, 'w');
    for i in range(len(quantity)):
        # Write the value associated with the triangle
        ofile.write("> -Z" + str(quantity[i]) + "\n");
        ofile.write(str(polygon_vertices[i, 0, 0]) + " " + str(polygon_vertices[i, 0, 1]) + "\n");
        ofile.write(str(polygon_vertices[i, 1, 0]) + " " + str(polygon_vertices[i, 1, 1]) + "\n");
        ofile.write(str(polygon_vertices[i, 2, 0]) + " " + str(polygon_vertices[i, 2, 1]) + "\n");
    ofile.close();
    return;


def write_single_eigenvector(positive_file, negative_file, e, v0, v1, x, y):
    # e = eigenvalue, [v0, v1] = eigenvector.
    # Writes a single eigenvector eigenvalue pair.
    # Also has functionality to saturate eigenvectors so they don't blow up.
    overall_max = 40.0;
    scale = 0.4 * e;

    vx = v0 * scale;
    vy = v1 * scale;
    if np.sqrt(vx * vx + vy * vy) > overall_max:
        scale = scale * (overall_max / np.sqrt(vx * vx + vy * vy))
        vx = v0 * scale;
        vy = v1 * scale;

    if e > 0:
        positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy));
        positive_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy));
    else:
        negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, vx, vy));
        negative_file.write("%s %s %s %s 0 0 0\n" % (x, y, -vx, -vy));
    return;
