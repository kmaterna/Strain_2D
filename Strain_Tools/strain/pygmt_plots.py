import pygmt
import numpy as np


def station_vels_to_arrays(vectors):
    """ Unpack a list of station_vels vectors into arrays for pygmt plotting """
    elon, nlat, e, n = [], [], [], []
    for item in vectors:
        elon.append(item.elon)
        nlat.append(item.nlat)
        e.append(item.e)
        n.append(item.n)
    return np.array(elon), np.array(nlat), np.array(e), np.array(n)


def filter_vectors_to_land_only(region, elon, nlat, e, n):
    if region[1] - region[0] < 0.2:
        maskfile = pygmt.grdlandmask(region=region, spacing='1s', resolution='h')
    else:
        maskfile = pygmt.grdlandmask(region=region, spacing='10m', resolution='i')
    points, newelon, newnlat, newe, newn = [], [], [], [], []
    for i in range(len(elon)):
        points.append(np.array([elon[i], nlat[i]]))   # build np array of points
    points = np.array(points)
    if len(points) == 0:
        return [], [], [], []
    values = pygmt.grdtrack(grid=maskfile, points=points)  # values is a pandas DataFrame
    for i, item in enumerate(values[2]):   # column 2 is the grdtrack output
        if item > 0.8:  # if grdtrack gives something on land
            newelon.append(elon[i])
            newnlat.append(nlat[i])
            newe.append(e[i])
            newn.append(n[i])
    return newelon, newnlat, newe, newn


def plot_rotation(rotation_array, station_vels, region, outdir, outfile):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="magma", series="0/300/1", truncate="0.3/1.0", background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Rotation\"")
    fig.grdimage(rotation_array, region=region, cmap=outdir+"/mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50", frame="1.0")
    if station_vels:
        elon, nlat, e, n = station_vels_to_arrays(station_vels)
        fig.plot(x=elon, y=nlat, style='c0.04i', fill='black', pen='0.4p,white')  # station locations
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
                 direction=[e, n])  # displacement vectors
        fig.plot(x=region[0], y=region[2], style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
                 direction=[[20], [0]], offset="0.9i/0.1i")  # scale vector
        fig.text(x=region[0], y=region[2], text="20 mm/yr", font='10p,Helvetica,black', offset='0.4i/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="0/300",
                 frame=["x50", "y+L\"Rad/Ka\""])
    print("Saving rotation figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_dilatation(dila_array, station_vels, region, outdir, outfile, positive_eigs=(), negative_eigs=()):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="polar", series="-200/200/2", reverse=True, background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Dilatation\"")
    fig.grdimage(dila_array, region=region, cmap=outdir+"/mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50", frame="1.0")
    if station_vels:
        elon, nlat, e, n = station_vels_to_arrays(station_vels)
        fig.plot(x=elon, y=nlat, style='c0.02i', fill='black')  # station locations
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black',
             offset='0.4i/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="-200/200",
                 frame=["x50", "y+L\"Nanostr/yr\""])
    print("Saving dilatation figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_I2nd(I2_array, station_vels, region, outdir, outfile, positive_eigs=(), negative_eigs=()):
    plotting_array = np.log10(np.abs(I2_array))  # for plotting the map of second invariant
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="batlow", series="-1/5/0.1", background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Second Invariant\"")
    fig.grdimage(plotting_array, region=region, cmap=outdir+"/mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50", frame="1.0")
    if station_vels:
        elon, nlat, e, n = station_vels_to_arrays(station_vels)
        fig.plot(x=elon, y=nlat, style='c0.02i', fill='black')  # station locations
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black', offset='0.4i/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="-1/5",
                 frame=["x1", "y+L\"Log10(I2)\""])
    print("Saving I2nd figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_maxshear(max_shear_array, station_vels, region, outdir, outfile, positive_eigs=(), negative_eigs=()):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="polar", series="0/300/2", truncate="0/1.0", background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Maximum Shear\"")
    fig.grdimage(max_shear_array, projection=proj, region=region, cmap=outdir+"/mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50", frame="1.0")
    if station_vels:
        elon, nlat, e, n = station_vels_to_arrays(station_vels)
        fig.plot(x=elon, y=nlat, style='c0.02i', fill='black')  # station locations
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black', offset='0.4i/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="0/300",
                 frame=["x50", "y+L\"Nanostr/yr\""])
    print("Saving MaxShear figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_azimuth(azimuth_array, station_vels, region, outdir, outfile, positive_eigs=(), negative_eigs=()):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="rainbow", series="0/180/1", background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Azimuth of Max Shortening\"")
    fig.grdimage(azimuth_array, region=region, cmap=outdir+"/mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.12/0.12+c" + str(region[2]) + "+w50", frame="1.0")
    if station_vels:
        elon, nlat, e, n = station_vels_to_arrays(station_vels)
        fig.plot(x=elon, y=nlat, style='c0.02i', fill='black')  # station locations
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black', offset='0.4i/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="0/180",
                 frame=["x30", "y+L\"Deg from North\""])
    print("Saving azimuth figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_dilatation_1D(region, polygon_outdir_file, outdir, outfile, positive_eigs=(), negative_eigs=()):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="polar", series="-200/200/2", reverse=True, background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Dilatation\"")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue', frame="1.0")

    fig.plot(data=polygon_outdir_file, pen="thinner,black", fill="+z", cmap=outdir+"/mycpt.cpt")
    fig.coast(borders='2', shorelines='1.0p,black', water='lightblue', map_scale="n0.12/0.12+c"+str(region[2])+"+w50")
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black', offset="0.4i/0.1i")
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="-200/200",
                 frame=["x50", "y+L\"Nanostr/yr\""])
    print("Saving dilatation figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_I2nd_1D(region, second_inv_polygon_file, outdir, outfile, positive_eigs=(), negative_eigs=()):
    proj = 'M4i'
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="batlow", series="-1/5/0.1", background="o", output=outdir+"/mycpt.cpt")
    fig.basemap(region=region, projection=proj, frame="+t\"Second Invariant\"")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue', frame="1.0")

    # color by value
    fig.plot(data=second_inv_polygon_file, pen="thinner,black", fill="+z", cmap=outdir + "/mycpt.cpt")
    fig.coast(borders='2', shorelines='1.0p,black', water='lightblue', map_scale="n0.12/0.12+c"+str(region[2])+"+w50")
    if positive_eigs:
        elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue',
                 direction=[e, n])
    if negative_eigs:
        elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
        elon, nlat, e, n = filter_vectors_to_land_only(region, elon, nlat, e, n)
        fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black',
                 direction=[e, n])
    # Scale vector
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]], offset="0.9i/0.1i")
    fig.plot(x=region[0], y=region[2], style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]], offset="0.9i/0.1i")
    fig.text(x=region[0], y=region[2], text="200 ns/yr", font='10p,Helvetica,black', offset="0.4i/0.1i")
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap=outdir+"/mycpt.cpt", truncate="-1/5",
                 frame=["x1", "y+L\"Log(I2nd)\""])
    print("Saving I2nd figure as %s." % outfile)
    fig.savefig(outfile)
    return


def plot_method_differences(strain_values_ds, average_strains, region, outdir, outfile):
    """Useful for dilatation and max shear based on values in the color bar"""
    pygmt.makecpt(cmap="polar", series="-300/300/2", truncate="-1.0/1.0", background="o", output=outdir+"/mycpt.cpt")
    fig = pygmt.Figure()
    proj = 'M2.2i'
    numrows = 2
    numcols = 2
    with fig.subplot(nrows=numrows, ncols=numcols, figsize=("7i", "6i"), frame="lrtb"):
        for i in range(numrows):  # row number starting from 0
            for j in range(numcols):  # column number starting from 0
                index = i * numcols + j  # index number starting from 0
                with fig.set_panel(panel=index):  # sets the current panel
                    fig.basemap(region=region, projection=proj, frame=["WeSn", "2.0"])
                    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue')
                    for counter, (varname, da) in enumerate(strain_values_ds.data_vars.items()):
                        if counter == index:
                            plotting_data = np.subtract(da, average_strains)  # testing
                            fig.grdimage(plotting_data, projection=proj, region=region, cmap=outdir+"/mycpt.cpt")
                            fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black',
                                      water='lightblue')
                            fig.text(position="BL", text=varname+" minus mean", region=region, offset='0/0.1i')
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/-0.5i", cmap=outdir+"/mycpt.cpt", truncate="-300/300",
                 frame=["x50", "y+L\"Nanostr/yr\""])
    print("Saving Method Differences as %s." % outfile)
    fig.savefig(outfile)
