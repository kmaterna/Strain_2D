import pygmt
import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write


def station_vels_to_arrays(station_vels):
    """ Unpack station vels into arrays for pygmt plotting of vectors """
    elon, nlat, e, n = [0], [0], [0], [0];
    for item in station_vels:
        elon.append(item.elon);
        nlat.append(item.nlat);
        e.append(item.e)
        n.append(item.n);
    return np.array(elon), np.array(nlat), np.array(e), np.array(n);


def plot_rotation(filename, station_vels, region, outdir, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="magma", T="0/300/1", G="0.3/1.0", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Rotation\"");
    fig.grdimage(filename, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    elon, nlat, e, n = station_vels_to_arrays(station_vels);
    fig.plot(x=elon, y=nlat, S='c0.04i', G='black', W='0.4p,white');  # station locations
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black',
             direction=[e, n]);  # displacement vectors
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
    pygmt.makecpt(C="polar", T="-200/200/2", I=True, D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Dilatation\"");
    fig.grdimage(filename, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="-200/200", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving dilatation figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_I2nd(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="batlow", T="-1/5/0.1", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Second Invariant\"");
    fig.grdimage(filename, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="-1/5", B=["x1", "y+L\"Log(I2)\""]);
    print("Saving I2nd figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_maxshear(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="polar", T="0/300/2", G="0/1.0", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Maximum Shear\"");
    fig.grdimage(filename, projection=proj, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="0/300", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving MaxShear figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_azimuth(filename, region, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="rainbow", T="0/180/1", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Azimuth of Max Shortening\"");
    fig.grdimage(filename, region=region, C=outdir+"/mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.12/0.12+c" + str(region[2]) + "+w50", B="1.0");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="0/180", B=["x30", "y+L\"Deg from North\""]);
    print("Saving azimuth figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_dilatation_1D(region, polygon_vertices, dilatation, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="polar", T="-200/200/2", I=True, D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Dilatation\"");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', B="1.0");

    # color by value
    for i in range(len(dilatation)):
        lons = [polygon_vertices[i, 0, 0], polygon_vertices[i, 1, 0], polygon_vertices[i, 2, 0]];
        lats = [polygon_vertices[i, 0, 1], polygon_vertices[i, 1, 1], polygon_vertices[i, 2, 1]];
        fig.plot(x=lons, y=lats, Z=str(dilatation[i]), pen="thinner,black", G="+z", C=outdir+"/mycpt.cpt");

    fig.coast(N='2', W='1.0p,black', S='lightblue', L="n0.12/0.12+c" + str(region[2]) + "+w50");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="-200/200", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving dilatation figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_I2nd_1D(region, polygon_vertices, I2nd, outdir, positive_eigs, negative_eigs, outfile):
    proj = 'M4i'
    fig = pygmt.Figure();
    pygmt.makecpt(C="batlow", T="-1/5/0.1", D="o", H=outdir+"/mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Second Invariant\"");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', B="1.0");

    # color by value
    for i in range(len(I2nd)):
        lons = [polygon_vertices[i, 0, 0], polygon_vertices[i, 1, 0], polygon_vertices[i, 2, 0]];
        lats = [polygon_vertices[i, 0, 1], polygon_vertices[i, 1, 1], polygon_vertices[i, 2, 1]];
        fig.plot(x=lons, y=lats, Z=str(I2nd[i]), pen="thinner,black", G="+z", C=outdir+"/mycpt.cpt");

    fig.coast(N='2', W='1.0p,black', S='lightblue', L="n0.12/0.12+c" + str(region[2]) + "+w50");
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[e, n]);
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs);
    fig.plot(x=elon, y=nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[e, n]);
    # Scale vector
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[200], [0]]);
    fig.plot(x=region[0] + 1.1, y=region[2] + 0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3',
             pen='0.6p,black', direction=[[-200], [0]]);
    fig.text(x=region[0] + 0.4, y=region[2] + 0.1, text="200 ns/yr", font='10p,Helvetica,black');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C=outdir+"/mycpt.cpt", G="-1/5", B=["x1", "y+L\"Log(I2nd)\""]);
    print("Saving I2nd figure as %s." % outfile)
    fig.savefig(outfile);
    return;


def plot_method_differences(strain_dictionary, average_strains, region, outdir, outfile):
    """Useful for dilatation and max shear based on values in the color bar"""
    pygmt.makecpt(C="polar", T="-300/300/2", G="-1.0/1.0", D="o", H=outdir+"/mycpt.cpt");
    fig = pygmt.Figure();
    proj = 'M2.2i'
    numrows = 2;
    numcols = 2;
    with fig.subplot(nrows=numrows, ncols=numcols, figsize=("7i", "6i"), frame="lrtb"):
        for i in range(numrows):  # row number starting from 0
            for j in range(numcols):  # column number starting from 0
                index = i * numcols + j  # index number starting from 0
                with fig.set_panel(panel=index):  # sets the current panel
                    fig.basemap(region=region, projection=proj, B=["WeSn", "2.0"]);
                    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue');
                    for counter, name in enumerate(strain_dictionary.keys()):
                        if counter == index:
                            plotting_data = np.subtract(strain_dictionary[name][2], average_strains);
                            netcdf_read_write.produce_output_netcdf(strain_dictionary[name][0],
                                                                    strain_dictionary[name][1], plotting_data,
                                                                    'per year', outdir + "/temp.grd");
                            fig.grdimage(outdir+'/temp.grd', projection=proj, region=region, C=outdir+"/mycpt.cpt");
                            fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue');
                            fig.text(position="BL", text=name+" minus mean", region=region, D='0/0.1i');
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/-0.5i", C=outdir+"/mycpt.cpt", G="-300/300", B=["x50", "y+L\"Nanostr/yr\""]);
    print("Saving Method Differences as %s." % outfile);
    fig.savefig(outfile);
    return;
