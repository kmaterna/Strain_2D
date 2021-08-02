# The output manager for Strain analysis.
# ----------------- OUTPUTS -------------------------
import numpy as np
import os
from xarray import Dataset

from . import strain_tensor_toolbox, velocity_io, pygmt_plots


def outputs_2d(xdata, ydata, rot, exx, exy, eyy, MyParams, myVelfield):
    print("------------------------------\nWriting 2d outputs:");
    velocity_io.write_stationvels(myVelfield, MyParams.outdir+"tempgps.txt");
    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy);
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy);

    # First create an xarray data multi-cube to write
    ds = Dataset(
        {
            "exx": (("y", "x"), exx),
            "eyy": (("y", "x"), eyy),
            "exy": (("y", "x"), exy),
            "azimuth":  (("y", "x"), azimuth),
            "rotation":  (("y", "x"), rot),
            "I2":  (("y", "x"), I2nd),
            "dilatation":  (("y", "x"), dilatation),
            "max_shear":  (("y", "x"), max_shear),
        },
        coords={
            "x": ('x', xdata),
            "y": ('y', ydata),
        },
    )
    ds.to_netcdf(os.path.join(MyParams.outdir, '{}_strain.nc'.format(MyParams.strain_method)))

    print("Max I2: %f " % (np.amax(I2nd)));
    print("Min/Max rot:   %f,   %f " % (np.nanmin(rot), np.nanmax(rot)) );

    # get grid eigenvectors for plotting
    [positive_eigs, negative_eigs] = get_grid_eigenvectors(xdata, ydata, e1, e2, v00, v01, v10, v11);
    velocity_io.write_gmt_format(positive_eigs, MyParams.outdir + 'positive_eigs.txt');
    velocity_io.write_gmt_format(negative_eigs, MyParams.outdir + 'negative_eigs.txt');

    # PYGMT PLOTS
    pygmt_plots.plot_rotation(ds['rotation'], myVelfield, MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir+'rotation.png');
    pygmt_plots.plot_dilatation(ds['dilatation'], MyParams.range_strain, MyParams.outdir, positive_eigs,
                                negative_eigs, MyParams.outdir+'dilatation.png');
    pygmt_plots.plot_I2nd(ds['I2'], MyParams.range_strain, MyParams.outdir, positive_eigs,
                          negative_eigs, MyParams.outdir+'I2nd.png');
    pygmt_plots.plot_maxshear(ds['max_shear'], MyParams.range_strain, MyParams.outdir, positive_eigs,
                              negative_eigs, MyParams.outdir+'max_shear.png');
    pygmt_plots.plot_azimuth(ds['azimuth'], MyParams.range_strain, MyParams.outdir, positive_eigs,
                             negative_eigs, MyParams.outdir+'azimuth.png');


def outputs_1d(xcentroid, ycentroid, polygon_vertices, rot, exx, exy, eyy, range_strain, myVelfield, outdir):
    print("------------------------------\nWriting 1d outputs:");
    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy);
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy);
    [positive_eigs, negative_eigs] = get_list_eigenvectors(xcentroid, ycentroid, e1, e2, v00, v01, v10, v11);

    velocity_io.write_multisegment_file(polygon_vertices, rot, outdir + "rot_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, I2nd, outdir + "I2nd_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, dilatation, outdir + "Dilatation_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, max_shear, outdir + "max_shear_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, azimuth, outdir + "azimuth_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, exx, outdir + "exx_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, exy, outdir + "exy_polygons.txt");
    velocity_io.write_multisegment_file(polygon_vertices, eyy, outdir + "eyy_polygons.txt");
    velocity_io.write_stationvels(myVelfield, outdir+"tempgps.txt");

    # Write the eigenvectors and eigenvalues
    velocity_io.write_gmt_format(positive_eigs, outdir + 'positive_eigs_polygons.txt');
    velocity_io.write_gmt_format(negative_eigs, outdir + 'negative_eigs_polygons.txt');

    print("Max I2: %f " % (max(I2nd)));
    print("Min/Max rot:   %f,   %f " % (np.amin(rot), np.amax(rot)) );

    # Plot the polygons as additional output (more intuitive)
    pygmt_plots.plot_dilatation_1D(range_strain, polygon_vertices, dilatation, outdir, positive_eigs,
                                   negative_eigs, outdir+'polygon_dilatation.eps');
    pygmt_plots.plot_I2nd_1D(range_strain, polygon_vertices, I2nd, outdir, positive_eigs,
                             negative_eigs, outdir+'polygon_I2nd.eps');


def get_grid_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11):
    """
    Resamples eigenvectors on regular grid, with maximum eigenvalue imposed
    Returns two lists of "stationvels" objects for plotting vectors

    :param xdata: 1d array of floats (lons)
    :param ydata: 1d array of floats (lats)
    :param w1: 2d arrays of floats (values)
    :param w2: 2d arrays of floats (values)
    :param v00: 2d arrays of floats
    :param v01: 2d arrays of floats
    :param v10: 2d arrays of floats
    :param v11: 2d arrays of floats
    """
    eigs_dec = 12;
    do_not_print_value = 200;
    overmax_scale = 200;
    positive_eigs, negative_eigs = [], [];
    for j in range(len(ydata)):
        for k in range(len(xdata)):
            if np.mod(j, eigs_dec) == 0 and np.mod(k, eigs_dec) == 0:
                if np.isnan(w1[j][k]) or np.isnan(w2[j][k]) or np.isnan(v00[j][k]) or np.isnan(v11[j][k]):
                    continue;
                # Write the first eigenvector pair
                scale = w1[j][k];
                if abs(scale) > do_not_print_value:
                    scale = overmax_scale;
                if w1[j][k] > 0:
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v00[j][k]*scale,
                                                                n=v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v00[j][k]*scale,
                                                                n=-v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                else:
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v00[j][k]*scale,
                                                                n=v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v00[j][k]*scale,
                                                                n=-v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                # Write the second eigenvector pair
                scale = w2[j][k];
                if abs(scale) > do_not_print_value:
                    scale = overmax_scale;
                if w2[j][k] > 0:
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v01[j][k]*scale,
                                                                n=v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v01[j][k]*scale,
                                                                n=-v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                else:
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v01[j][k]*scale,
                                                                n=v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v01[j][k]*scale,
                                                                n=-v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0));
    return positive_eigs, negative_eigs;


def get_list_eigenvectors(xdata, ydata, w1, w2, v00, v01, v10, v11):
    """
    Returns vectors with maximum eigenvalue imposed
    Returns two lists of "stationvels" objects for plotting vectors

    :param xdata: 1d array of floats (lons)
    :param ydata: 1d array of floats (lats)
    :param w1: 1d array of floats (values)
    :param w2: 1d array of floats (values)
    :param v00: 1d array of floats (values)
    :param v01: 1d array of floats (values)
    :param v10: 1d array of floats (values)
    :param v11: 1d array of floats (values)
    """
    positive_eigs, negative_eigs = [], [];
    overall_max = 40.0;
    for i in range(len(xdata)):
        scale = 0.4 * w1[i];
        if abs(w1[i]) > overall_max:
            scale = overall_max;
        if np.isnan(w1[i]) or np.isnan(w2[i]) or np.isnan(v00[i]) or np.isnan(v11[i]):
            continue;
        if w1[i] > 0:
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v00[i] * scale,
                                                        n=v10[i] * scale, u=0, se=0, sn=0, su=0, name=0));
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v00[i] * scale,
                                                        n=-v10[i] * scale, u=0, se=0, sn=0, su=0, name=0));
        else:
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v00[i] * scale,
                                                        n=v10[i] * scale, u=0, se=0, sn=0, su=0, name=0));
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v00[i] * scale,
                                                        n=-v10[i] * scale, u=0, se=0, sn=0, su=0, name=0));
        scale = 0.4 * w2[i];
        if abs(w2[i]) > overall_max:
            scale = overall_max;
        if w2[i] > 0:
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v01[i] * scale,
                                                        n=v11[i] * scale, u=0, se=0, sn=0, su=0, name=0));
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v01[i] * scale,
                                                        n=-v11[i] * scale, u=0, se=0, sn=0, su=0, name=0));
        else:
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v01[i] * scale,
                                                        n=v11[i] * scale, u=0, se=0, sn=0, su=0, name=0));
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v01[i] * scale,
                                                        n=-v11[i] * scale, u=0, se=0, sn=0, su=0, name=0));
    return positive_eigs, negative_eigs;
