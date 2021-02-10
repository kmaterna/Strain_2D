# The output manager for GPS Strain analysis. 
# ----------------- OUTPUTS -------------------------

import numpy as np
import os, subprocess
import pygmt
from Tectonic_Utils.read_write import netcdf_read_write
from Strain_2D.strain import strain_tensor_toolbox
from Strain_2D.strain import velocity_io
from Strain_2D.strain import configure_functions


def outputs_2d(xdata, ydata, rot, exx, exy, eyy, MyParams, myVelfield):
    print("------------------------------\nWriting 2d outputs:");
    velocity_io.write_stationvels(myVelfield, MyParams.outdir+"tempgps.txt");
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
    print("Min/Max rot:   %f,   %f " % (np.amin(rot), np.amax(rot)) );
    write_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11, MyParams);
    positive_eigs = velocity_io.read_horiz_vels(MyParams.outdir+"positive_eigs.txt");
    negative_eigs = velocity_io.read_horiz_vels(MyParams.outdir+"negative_eigs.txt");

    # PYGMT PLOTS
    plot_rotation(MyParams.outdir+'rot.nc', myVelfield, MyParams.range_strain, MyParams.outdir,
                  MyParams.outdir+'rotation.png');
    plot_dilatation(MyParams.outdir+'dila.nc', MyParams.range_strain, MyParams.outdir, positive_eigs, negative_eigs,
                    MyParams.outdir+'dilatation.png');
    plot_I2nd(MyParams.outdir+'I2nd.nc', MyParams.range_strain, MyParams.outdir, positive_eigs, negative_eigs,
              MyParams.outdir+'I2nd.png');
    plot_maxshear(MyParams.outdir+'max_shear.nc', MyParams.range_strain, MyParams.outdir, positive_eigs, negative_eigs,
                  MyParams.outdir+'max_shear.png');
    plot_azimuth(MyParams.outdir+'azimuth.nc', MyParams.range_strain, MyParams.outdir, positive_eigs, negative_eigs,
                 MyParams.outdir+'azimuth.png');
    return;


def write_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11, MyParams):
    # Need eigs_interval and outdir from MyParams.
    positive_file = open(MyParams.outdir + "positive_eigs.txt", 'w');
    negative_file = open(MyParams.outdir + "negative_eigs.txt", 'w');
    
    eigs_dec = 12;
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
    velocity_io.write_stationvels(myVelfield, MyParams.outdir+"tempgps.txt");

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

    # Plot the polygons as additional output (more intuitive)
    curr_directory = os.getcwd();
    dir_path = os.path.dirname(os.path.realpath(__file__));
    range_string = configure_functions.get_string_range(MyParams.range_data);
    os.chdir(MyParams.outdir);
    subprocess.call([dir_path + '/GMT_mapping_codes/delaunay_gmt.gmt', range_string], shell=False);
    os.chdir(curr_directory);
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


def plot_rotation(filename, station_vels, region, outdir, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="magma", T="0/300/1", G="0.3/1.0", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Rotation\"");
    fig.grdimage(filename, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    for n in station_vels:
        fig.plot(x=n.elon, y=n.nlat, S='c0.04i', G='black', W='0.4p,white');  # station locations
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
                 direction=[[n.e], [n.n]]);  # displacement vectors
    fig.plot(x=region[0] + 0.9, y=region[2] + 0.1, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
             direction=[[20], [0]]);  # scale vector
    fig.text(x=region[0] + 0.5, y=region[2] + 0.1, text="20 mm/yr", font='10p,Helvetica,black')
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="0/300", B=["x50", "y+L\"Rad/Ka\""]);
    print("Saving rotation figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_dilatation(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="polar", T="-200/200/2", I=True, D="o", H=outdir+"mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Dilatation\"");
    fig.grdimage(filename, region=region, C=outdir+"mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    for n in positive_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[[n.e], [n.n]]);  # vectors
    for n in negative_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[[n.e], [n.n]]);  # vectors
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"mycpt.cpt", G="-200/200", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving dilatation figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_I2nd(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="batlow", T="-1/5/0.1", D="o", H=outdir+"mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Second Invariant\"");
    fig.grdimage(filename, region=region, C=outdir+"mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    for n in positive_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[[n.e], [n.n]]);  # vectors
    for n in negative_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[[n.e], [n.n]]);  # vectors
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"mycpt.cpt", G="-1/5", B=["x1", "y+L\"Log(I2)\""]);
    print("Saving I2nd figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_maxshear(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="polar", T="0/300/2", G="0/1.0", D="o", H=outdir+"mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Maximum Shear\"");
    fig.grdimage(filename, region=region, C=outdir+"mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    for n in positive_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[[n.e], [n.n]]);  # vectors
    for n in negative_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[[n.e], [n.n]]);  # vectors
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"mycpt.cpt", G="0/300", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving MaxShear figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_azimuth(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="rainbow", T="0/180/1", D="o", H=outdir+"mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Azimuth of Max Shortening\"");
    fig.grdimage(filename, region=region, C=outdir+"mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    for n in positive_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[[n.e], [n.n]]);  # vectors
    for n in negative_eigs:
        fig.plot(x=n.elon, y=n.nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[[n.e], [n.n]]);  # vectors
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"mycpt.cpt", G="0/180", B=["x30", "y+L\"Deg from North\""]);
    print("Saving azimuth figure as %s." % outfile)
    fig.savefig(outfile);
    return;
