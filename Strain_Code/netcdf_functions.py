import numpy as np 
import scipy.io.netcdf as netcdf
from subprocess import check_output


def read_grd_z(filename, zname):
	zdata0 = netcdf.netcdf_file(filename,'r').variables[zname][::-1];
	zdata  = zdata0.copy();	
	return [zdata];

def read_grd_xy(filename, xname, yname):
	xdata0 = netcdf.netcdf_file(filename,'r').variables[xname][::-1];
	ydata0 = netcdf.netcdf_file(filename,'r').variables[yname][::-1];
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	return [xdata, ydata]; 

def read_grd_xyz(filename, xname, yname, zname):
	# xdata0 = netcdf.netcdf_file(filename,'r').variables[xname][::-1];
	# ydata0 = netcdf.netcdf_file(filename,'r').variables[yname][::-1];
	# zdata0 = netcdf.netcdf_file(filename,'r').variables[zname][::-1];
	xdata0 = netcdf.netcdf_file(filename,'r').variables[xname][::];
	ydata0 = netcdf.netcdf_file(filename,'r').variables[yname][::];
	zdata0 = netcdf.netcdf_file(filename,'r').variables[zname][::];
	zdata=zdata0.copy();	
	xdata=xdata0.copy();
	ydata=ydata0.copy();
	return [xdata, ydata, zdata]; 





def produce_output_netcdf(xdata, ydata, zdata, zunits, netcdfname):
	# # Write the netcdf velocity grid file.  
	f=netcdf.netcdf_file(netcdfname,'w');
	f.history = 'Created for a test';
	f.createDimension('x',len(xdata));
	f.createDimension('y',len(ydata));
	print(np.shape(zdata));
	x=f.createVariable('x',float,('x',))
	x[:]=xdata;
	x.units = 'longitude';
	y=f.createVariable('y',float,('y',))
	y[:]=ydata;
	y.units = 'latitude';
	z=f.createVariable('z',float,('y','x',));
	z[:,:]=zdata;
	z.units = zunits;
	f.close();
	return;

def flip_if_necessary(filename):
	# IF WE NEED TO FLIP DATA:
	xinc = check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $8}\'',shell=True);  # the x-increment
	yinc = check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $9}\'',shell=True);  # the y-increment
	xinc=float(xinc.split()[0]);
	yinc=float(yinc.split()[0]);

	if xinc < 0:  # FLIP THE X-AXIS
		print("flipping the x-axis");
		[xdata,ydata] = read_grd_xy(filename, 'x', 'y');
		[data] = read_grd_z(filename,'z');
		# This is the key! Flip the x-axis when necessary.  
		#xdata=np.flip(xdata,0);  # This is sometimes necessary and sometimes not!  Not sure why. 
		produce_output_netcdf(xdata, ydata, data, '',filename);
		xinc = check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $8}\'',shell=True);  # the x-increment
		xinc = float(xinc.split()[0]);
		print("New xinc is: %f " % (xinc) );
	if yinc < 0:
		print("flipping the y-axis");
		[xdata,ydata] = read_grd_xy(filename,'x','y');
		[data] = read_grd_z(filename,'z');
		# Flip the y-axis when necessary.  
		# ydata=np.flip(ydata,0);  
		produce_output_netcdf(xdata, ydata, data, '',filename);
		yinc = check_output('gmt grdinfo -M -C '+filename+' | awk \'{print $9}\'',shell=True);  # the x-increment
		yinc = float(yinc.split()[0]);
		print("New yinc is: %f" % (yinc) );
	return;
