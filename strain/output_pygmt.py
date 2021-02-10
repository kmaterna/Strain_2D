import pygmt
import velocity_io

# Configure and input
region = [-125, -120, 38, 42]
proj = "M4i"
station_vels = velocity_io.read_stationvels("tempgps.txt");
positive_eigs = velocity_io.read_horiz_vels("positive_eigs.txt");
negative_eigs = velocity_io.read_horiz_vels("negative_eigs.txt");

# # Outputs for Rotation
# infile = 'rot.nc'
# outfile = 'rotation.png'
# fig = pygmt.Figure();
# pygmt.makecpt(C="magma", T="0/300/1", G="0.3/1.0", D="o", H="mycpt.cpt");
# fig.basemap(region=region, projection=proj, B="+t\"Rotation\"");
# fig.grdimage(infile, region=region, C="mycpt.cpt");
# fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', L="n0.15/0.1+c"+str(region[2])+"+w50",
#           B="1.0")
# for n in station_vels:
#     fig.plot(x=n.elon, y=n.nlat, S='c0.04i',G='black',W='0.4p,white');  # station locations
#     fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black', direction=[[n.e], [n.n]]);  # vectors!
# # Scale Vector
# fig.plot(x=region[0]+0.8, y=region[2]+0.1, style='v0.20+e+a40+gblack+h0+p1p,black+z0.04', pen='0.6p,black', direction=[[20], [0]]);
# fig.text(x=region[0]+0.4, y=region[2]+0.1, text="20 mm/yr",font='11p,Helvetica,black')
# fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="0/300", B=["x50", "y+L\"Rad/Ka\""]);
# print("Saving rotation figure as %s." % outfile)
# fig.savefig(outfile);


# # Outputs for Dilatation
# infile = 'dila.nc'
# outfile = 'dilatation.png'
# fig = pygmt.Figure();
# pygmt.makecpt(C="polar", T="-200/200/2", I=True, D="o", H="mycpt.cpt");
# fig.basemap(region=region, projection=proj, B="+t\"Dilatation\"");
# fig.grdimage(infile, region=region, C="mycpt.cpt");
# fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', L="n0.12/0.12+c"+str(region[2])+"+w50",
#           B="1.0")
# for n in positive_eigs:
#     fig.plot(x=n.elon, y=n.nlat, style='v0.20+e+a40+gblue+h0.5+p0.3p,blue+z0.003+n0.3', pen='0.6p,blue', direction=[[n.e], [n.n]]);  # vectors!
# for n in negative_eigs:
#     fig.plot(x=n.elon, y=n.nlat, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[[n.e], [n.n]]);  # vectors!
# # Scale vector
# fig.plot(x=region[0]+1.1, y=region[2]+0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[[200], [0]]);
# fig.plot(x=region[0]+1.1, y=region[2]+0.1, style='v0.20+b+a40+gred+h0.5+p0.3p,black+z0.003+n0.3', pen='0.6p,black', direction=[[-200], [0]]);
# fig.text(x=region[0]+0.4, y=region[2]+0.1, text="200 ns/yr", font='10p,Helvetica,black');
# fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="-200/200", B=["x50", "y+L\"Nanostr/yr\""]);
# print("Saving dilatation figure as %s." % outfile)
# fig.savefig(outfile);

#
# # Outputs for Second Invariant
# infile = 'I2nd.nc'
# outfile = 'I2nd.png'
# fig = pygmt.Figure();
# pygmt.makecpt(C="batlow", T="-1/5/0.1", D="o", H="mycpt.cpt");
# fig.basemap(region=region, projection=proj, B="+t\"Second Invariant\"");
# fig.grdimage(infile, region=region, C="mycpt.cpt");
# fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', L="n0.15/0.1+c"+str(region[2])+"+w50",
#           B="1.0")
# fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="-1/5", B=["x1", "y+L\"Log(I2)\""]);
# print("Saving I2nd figure as %s." % outfile)
# fig.savefig(outfile);
#
#
# # Outputs for Max Shear
# infile = 'max_shear.nc'
# outfile = 'max_shear.png'
# fig = pygmt.Figure();
# pygmt.makecpt(C="polar", T="0/300/2", G="0/1.0", D="o", H="mycpt.cpt");
# fig.basemap(region=region, projection=proj, B="+t\"Maximum Shear\"");
# fig.grdimage(infile, region=region, C="mycpt.cpt");
# fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', L="n0.15/0.1+c"+str(region[2])+"+w50",
#           B="1.0")
# fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="0/300", B=["x50", "y+L\"Nanostr/yr\""]);
# print("Saving MaxShear figure as %s." % outfile)
# fig.savefig(outfile);
#
#
# # Outputs for Azimuth
# infile = 'azimuth.nc'
# outfile = 'azimuth.png'
# fig = pygmt.Figure();
# pygmt.makecpt(C="rainbow", T="0/180/1", D="o", H="mycpt.cpt");
# fig.basemap(region=region, projection=proj, B="+t\"Azimuth of Max Shortening\"");
# fig.grdimage(infile, region=region, C="mycpt.cpt");
# fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue', L="n0.15/0.1+c"+str(region[2])+"+w50",
#           B="1.0")
# fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="0/180", B=["x30", "y+L\"Deg from North\""]);
# print("Saving azimuth figure as %s." % outfile)
# fig.savefig(outfile);
