# The output manager for GPS Strain analysis. 
# ----------------- OUTPUTS -------------------------

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
import strain_tensor_toolbox
import gps_io_functions


def outputs_2d(xdata, ydata, rot, exx, exy, eyy, MyParams, myVelfield):
    print("------------------------------\nWriting 2d outputs:");
    gps_io_functions.write_humanread_vel_file(myVelfield, MyParams.outdir+"tempgps.txt");
    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy);
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy);
    netcdf_read_write.produce_output_netcdf(xdata, ydata, exx, 'microstrain', MyParams.outdir + 'exx.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, exy, 'microstrain', MyParams.outdir + 'exy.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, eyy, 'microstrain', MyParams.outdir + 'eyy.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, azimuth, 'degrees', MyParams.outdir + 'azimuth.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, I2nd, 'per yr', MyParams.outdir + 'I2nd.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, rot, 'per yr', MyParams.outdir + 'rot.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, dilatation, 'per yr', MyParams.outdir + 'dila.nc');
    netcdf_read_write.produce_output_netcdf(xdata, ydata, max_shear, 'per yr', MyParams.outdir + 'max_shear.nc');
    print("Max I2: %f " % (np.amax(I2nd)));
    print("Max rot: %f " % (np.nanmax(rot)));
    print("Min rot: %f " % (np.nanmin(rot)));
    write_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11, MyParams);
    return;


def write_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11, MyParams):
    # Need eigs_interval and outdir from MyParams.
    positive_file = open(MyParams.outdir + "positive_eigs.txt", 'w');
    negative_file = open(MyParams.outdir + "negative_eigs.txt", 'w');
    if MyParams.strain_method == 'visr':
        eigs_dec = 8;
    elif MyParams.strain_method == 'gpsgridder':
        eigs_dec = 12;
    elif MyParams.strain_method == 'delaunay' or MyParams.strain_method == 'delaunay_flat':
        eigs_dec = 12;
    else:
        raise Exception("Error! strain method not recognized for eigenvector plotting.");

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


def outputs_1d(xcentroid, ycentroid, polygon_vertices, rot, exx, exy, eyy, myVelfield, MyParams):
    print("------------------------------\nWriting 1d outputs:");
    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy);
    write_multisegment_file(polygon_vertices, rot, MyParams.outdir+"rot_polygons.txt");
    write_multisegment_file(polygon_vertices, I2nd, MyParams.outdir+"I2nd_polygons.txt");
    write_multisegment_file(polygon_vertices, dilatation, MyParams.outdir+"Dilatation_polygons.txt");
    write_multisegment_file(polygon_vertices, max_shear, MyParams.outdir+"max_shear_polygons.txt");
    write_multisegment_file(polygon_vertices, azimuth, MyParams.outdir+"azimuth_polygons.txt");
    write_multisegment_file(polygon_vertices, exx, MyParams.outdir + "exx_polygons.txt");
    write_multisegment_file(polygon_vertices, exy, MyParams.outdir + "exy_polygons.txt");
    write_multisegment_file(polygon_vertices, eyy, MyParams.outdir + "eyy_polygons.txt");
    gps_io_functions.write_humanread_vel_file(myVelfield, MyParams.outdir+"tempgps.txt");

    positive_file = open(MyParams.outdir + "positive_eigs_polygons.txt", 'w');
    negative_file = open(MyParams.outdir + "negative_eigs_polygons.txt", 'w');
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy);
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
