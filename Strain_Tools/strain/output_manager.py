# The output manager for Strain analysis.
# ----------------- OUTPUTS -------------------------
import numpy as np
import os
from xarray import Dataset
from . import strain_tensor_toolbox, velocity_io, pygmt_plots, moment_functions, data_misfits


def outputs_2d(Ve, Vn, rot, exx, exy, eyy, MyParams, myVelfield, residfield):
    """Every strain method goes through this function at the end of its output stage"""
    print("------------------------------\nWriting 2d outputs:")

    # Write residual velocities.  Filter observations by range_strain bounding box.
    # TODO: Residual field should be identical to observed except for VE/VN/VU. Currently SE/SN/SU are not the same?
    # TODO: Alternatively, just everything to one file and have two extra columns with the residuals appended.
    velocity_io.write_stationvels(myVelfield, MyParams.outdir + 'obs_vels.txt', header='Obs Velocity.')
    velocity_io.write_stationvels(residfield, MyParams.outdir + 'residual_vels.txt', header='Obs-minus-model.')

    if len(myVelfield) != len(residfield):
        raise ValueError("Error! Velocity field and residual field have different lengths "
                         "("+str(len(myVelfield))+" vs "+str(len(residfield))+").")

    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy)
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy)

    # get grid eigenvectors for plotting
    [positive_eigs, negative_eigs] = get_grid_eigenvectors(MyParams.xdata, MyParams.ydata, e1, e2, v00, v01, v10, v11)
    velocity_io.write_gmt_format(positive_eigs, MyParams.outdir + 'positive_eigs.txt')
    velocity_io.write_gmt_format(negative_eigs, MyParams.outdir + 'negative_eigs.txt')

    # First create an xarray data multi-cube to write
    ds = Dataset(
        {
            "Ve": (("y", "x"), Ve),
            "Vn": (("y", "x"), Vn),
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
            "x": ('x', MyParams.xdata),
            "y": ('y', MyParams.ydata),
        },
    )

    output_filename = os.path.join(MyParams.outdir, '{}_strain.nc'.format(MyParams.strain_method))
    print("Writing file %s " % output_filename)
    ds.to_netcdf(output_filename)

    print("Max I2: %f " % (np.nanmax(I2nd)))
    print("Min/Max rot:   %f,   %f " % (np.nanmin(rot), np.nanmax(rot)))

    # PYGMT PLOTS
    pygmt_plots.plot_rotation(ds['rotation'], myVelfield, MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir+'rotation.png')
    pygmt_plots.plot_dilatation(ds['dilatation'], myVelfield, MyParams.range_strain, MyParams.outdir,
                                MyParams.outdir + 'dilatation.png', positive_eigs, negative_eigs)
    pygmt_plots.plot_I2nd(ds['I2'], myVelfield, MyParams.range_strain, MyParams.outdir, MyParams.outdir + 'I2nd.png',
                          positive_eigs, negative_eigs)
    pygmt_plots.plot_maxshear(ds['max_shear'], myVelfield, MyParams.range_strain, MyParams.outdir,
                              MyParams.outdir + 'max_shear.png', positive_eigs, negative_eigs)
    pygmt_plots.plot_azimuth(ds['azimuth'], myVelfield, MyParams.range_strain, MyParams.outdir,
                             MyParams.outdir + 'azimuth.png', positive_eigs, negative_eigs)

    if MyParams.write_metrics:  # optional: we can automatically compute a metric of the mag. of strain field
        output_params = {"outdir": MyParams.outdir,
                         "netcdf": output_filename,
                         "landmask": MyParams.outdir+"landmask.grd",
                         "mu": 30,
                         "depth": 11,
                         "outfile": MyParams.outdir+"strain_metrics.txt",
                         "obs_velfile": MyParams.outdir + 'obs_vels.txt',
                         "resid_velfile": MyParams.outdir + 'residual_vels.txt',
                         "use_landmask": 1}
        moment_functions.moment_coordinator(output_params)
        data_misfits.misfits_coordinator(output_params)  # append-mode

    return


def outputs_1d(xcentroid, ycentroid, polygon_vertices, rot, exx, exy, eyy, range_strain, myVelfield, outdir):
    print("------------------------------\nWriting 1d outputs:")
    exx, exy, eyy = np.array(exx), np.array(exy), np.array(eyy)
    [I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(exx, exy, eyy)
    [e1, e2, v00, v01, v10, v11] = strain_tensor_toolbox.compute_eigenvectors(exx, exy, eyy)
    [positive_eigs, negative_eigs] = get_list_eigenvectors(xcentroid, ycentroid, e1, e2, v00, v01, v10, v11)
    I2nd = np.log10(np.abs(I2nd))  # for convenient plotting

    dilatation_polygon_outfile = outdir + "Dilatation_polygons.txt"
    second_inv_polygon_outfile = outdir + "I2nd_polygons.txt"
    velocity_io.write_multisegment_file(polygon_vertices, rot, outdir + "rot_polygons.txt")
    velocity_io.write_multisegment_file(polygon_vertices, I2nd, second_inv_polygon_outfile)
    velocity_io.write_multisegment_file(polygon_vertices, dilatation, dilatation_polygon_outfile)
    velocity_io.write_multisegment_file(polygon_vertices, max_shear, outdir + "max_shear_polygons.txt")
    velocity_io.write_multisegment_file(polygon_vertices, azimuth, outdir + "azimuth_polygons.txt")
    velocity_io.write_multisegment_file(polygon_vertices, exx, outdir + "exx_polygons.txt")
    velocity_io.write_multisegment_file(polygon_vertices, exy, outdir + "exy_polygons.txt")
    velocity_io.write_multisegment_file(polygon_vertices, eyy, outdir + "eyy_polygons.txt")
    velocity_io.write_stationvels(myVelfield, outdir+"tempgps.txt")

    # Write the eigenvectors and eigenvalues
    velocity_io.write_gmt_format(positive_eigs, outdir + 'positive_eigs_polygons.txt')
    velocity_io.write_gmt_format(negative_eigs, outdir + 'negative_eigs_polygons.txt')

    print("Max I2: %f " % (max(I2nd)))
    print("Min/Max rot:   %f,   %f " % (np.amin(rot), np.amax(rot)))

    # Plot the polygons as additional output (more intuitive)
    pygmt_plots.plot_dilatation_1D(range_strain, dilatation_polygon_outfile, outdir,
                                   outdir + 'polygon_dilatation.png', positive_eigs, negative_eigs)
    pygmt_plots.plot_I2nd_1D(range_strain, second_inv_polygon_outfile, outdir, outdir + 'polygon_I2nd.png',
                             positive_eigs, negative_eigs)
    return


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
    eigs_dec = 12
    do_not_print_value = 200
    overmax_scale = 200
    positive_eigs, negative_eigs = [], []
    for j in range(len(ydata)):
        for k in range(len(xdata)):
            if np.mod(j, eigs_dec) == 0 and np.mod(k, eigs_dec) == 0:
                if np.isnan(w1[j][k]) or np.isnan(w2[j][k]) or np.isnan(v00[j][k]) or np.isnan(v11[j][k]):
                    continue
                # Write the first eigenvector pair
                scale = w1[j][k]
                if abs(scale) > do_not_print_value:
                    scale = overmax_scale
                if w1[j][k] > 0:
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v00[j][k]*scale,
                                                                n=v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v00[j][k]*scale,
                                                                n=-v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                else:
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v00[j][k]*scale,
                                                                n=v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v00[j][k]*scale,
                                                                n=-v10[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                # Write the second eigenvector pair
                scale = w2[j][k]
                if abs(scale) > do_not_print_value:
                    scale = overmax_scale
                if w2[j][k] > 0:
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v01[j][k]*scale,
                                                                n=v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                    positive_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v01[j][k]*scale,
                                                                n=-v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                else:
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=v01[j][k]*scale,
                                                                n=v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
                    negative_eigs.append(velocity_io.StationVel(elon=xdata[k], nlat=ydata[j], e=-v01[j][k]*scale,
                                                                n=-v11[j][k]*scale, u=0, se=0, sn=0, su=0, name=0))
    return positive_eigs, negative_eigs


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
    positive_eigs, negative_eigs = [], []
    overall_max = 40.0
    for i in range(len(xdata)):
        scale = 0.4 * w1[i]
        if abs(w1[i]) > overall_max:
            scale = overall_max
        if np.isnan(w1[i]) or np.isnan(w2[i]) or np.isnan(v00[i]) or np.isnan(v11[i]):
            continue
        if w1[i] > 0:
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v00[i] * scale,
                                                        n=v10[i] * scale, u=0, se=0, sn=0, su=0, name=0))
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v00[i] * scale,
                                                        n=-v10[i] * scale, u=0, se=0, sn=0, su=0, name=0))
        else:
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v00[i] * scale,
                                                        n=v10[i] * scale, u=0, se=0, sn=0, su=0, name=0))
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v00[i] * scale,
                                                        n=-v10[i] * scale, u=0, se=0, sn=0, su=0, name=0))
        scale = 0.4 * w2[i]
        if abs(w2[i]) > overall_max:
            scale = overall_max
        if w2[i] > 0:
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v01[i] * scale,
                                                        n=v11[i] * scale, u=0, se=0, sn=0, su=0, name=0))
            positive_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v01[i] * scale,
                                                        n=-v11[i] * scale, u=0, se=0, sn=0, su=0, name=0))
        else:
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=v01[i] * scale,
                                                        n=v11[i] * scale, u=0, se=0, sn=0, su=0, name=0))
            negative_eigs.append(velocity_io.StationVel(elon=xdata[i], nlat=ydata[i], e=-v01[i] * scale,
                                                        n=-v11[i] * scale, u=0, se=0, sn=0, su=0, name=0))
    return positive_eigs, negative_eigs
