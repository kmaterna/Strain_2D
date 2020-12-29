# Use GPS Gridder from GMT to interpolate between GPS stations
# The algorithm is based on the greens functions for elastic sheets with a given Poisson's ratio. 
# From: Sandwell, D. T., and P. Wessel (2016),
# Interpolation of 2-D vector data using constraints from elasticity, Geophys. Res.Lett. 


import numpy as np 
import subprocess
from Tectonic_Utils.read_write import netcdf_read_write
import strain_tensor_toolbox


# ----------------- COMPUTE -------------------------
def compute(myVelfield, MyParams):
	print("Computing strain via gpsgridder method.");
	outfile=open("tempgps.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();
	subprocess.call("gmt gpsgridder tempgps.txt -R"+MyParams.map_range+" -I"+str(MyParams.grid_inc)+" -S0.5 -Fd0.01 -C0.0005 -Emisfitfile.txt -fg -r -Gnc4_%s.nc",shell=True);  # makes a netcdf grid file
	# -R = range. -I = interval. -E prints the model and data fits at the input stations (very useful). 
	# -S = poisson's ratio. -Fd = fudge factor. -C = eigenvalues below this value will be ignored. 
	# -fg = flat earth approximation. -G = output netcdf files (x and y displacements). 
	# You should experiment with Fd and C values to find something that you like (good fit without overfitting). 
	# For Northern California, I like -Fd0.01 -C0.005. -R-125/-121/38/42.2
	# I had Ff0.1 -C0.005 -fg 
	
	# For large domains, GMT netcdf4 files instead of netcdf3. We must turn them all into netcdf3 for python to read them. 
	subprocess.call('nccopy -k classic nc4_u.nc gps_u.nc',shell=True);
	subprocess.call('nccopy -k classic nc4_v.nc gps_v.nc',shell=True);
	subprocess.call(['rm','tempgps.txt'],shell=False);
	subprocess.call(['rm','gmt.history'],shell=False);
	subprocess.call(['mv','misfitfile.txt',MyParams.outdir],shell=False);
	subprocess.call(['mv','nc4_u.nc',MyParams.outdir],shell=False);
	subprocess.call(['mv','nc4_v.nc',MyParams.outdir],shell=False);
	subprocess.call(['mv','gps_u.nc',MyParams.outdir],shell=False);
	subprocess.call(['mv','gps_v.nc',MyParams.outdir],shell=False);

	# Get ready to do strain calculation. 
	file1=MyParams.outdir+"gps_u.nc";
	file2=MyParams.outdir+"gps_v.nc";
	[xdata, ydata, udata] = netcdf_read_write.read_any_grd(file1);
	[_, _, vdata] = netcdf_read_write.read_any_grd(file2);
	xinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $8}\'',shell=True));  # the x-increment
	yinc = float(subprocess.check_output('gmt grdinfo -M -C '+file1+' | awk \'{print $9}\'',shell=True));  # the y-increment	
	typical_lat=float(MyParams.map_range[2]);
	xinc = xinc * 111.000 * np.cos(np.deg2rad(typical_lat)); # in km (not degrees)
	yinc = yinc * 111.000;  # in km (not degrees)
	xdata=np.flipud(xdata); ydata=np.flipud(ydata); udata=np.flipud(udata); vdata=np.flipud(vdata);
	[ydim, xdim] = np.shape(udata)
	rot=np.zeros(np.shape(vdata));  # 2nd invariant of rotation rate tensor
	e1=np.zeros(np.shape(vdata));  # maximum principal strain (array of float)
	e2=np.zeros(np.shape(vdata));  # minimum principal strain (array of float)
	v00=np.zeros(np.shape(vdata));  # more complicated: eigenvectors (array of matrix 2x2)
	v01=np.zeros(np.shape(vdata));  # more complicated: eigenvectors (array of matrix 2x2)
	v10=np.zeros(np.shape(vdata));  # more complicated: eigenvectors (array of matrix 2x2)
	v11=np.zeros(np.shape(vdata));  # more complicated: eigenvectors (array of matrix 2x2)

	# the strain calculation
	for j in range(ydim-1):
		for i in range(xdim-1):
			up=udata[j][i];
			vp=vdata[j][i];
			uq=udata[j][i+1];
			vq=vdata[j][i+1];
			ur=udata[j+1][i];
			vr=vdata[j+1][i];

			[dudx, dvdx, dudy, dvdy] = strain_tensor_toolbox.compute_displacement_gradients(up, vp, ur, vr, uq, vq, xinc, yinc);
			
			# The components that are easily computed
			# Units: nanostrain per year. 
			[exx, exy, eyy, rotation] = strain_tensor_toolbox.compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy);
			rot[j][i]=abs(rotation);

			# Compute a number of values based on tensor properties. 
			[e11, e22, v1] = strain_tensor_toolbox.eigenvector_eigenvalue(exx, exy, eyy);
			e1[j][i]= e11;
			e2[j][i]= e22;
			v00[j][i]=v1[0][0];
			v10[j][i]=v1[1][0];
			v01[j][i]=v1[0][1];
			v11[j][i]=v1[1][1];

	print("Success computing strain via gpsgridder method.\n");

	return [xdata, ydata, rot, e1, e2, v00, v01, v10, v11];







