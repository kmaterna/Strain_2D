# The purpose of this script is to take velocities, interpolate them with a spline function, and then compute the 
# resulting strain field. 

# NOTE: The Numpy interp2d function is a little finicky.  I need to include lots of extra data (all the way to Utah) 
# and use high-order quintic spline to fit the data without messy singularities. 
# If it's still not working, you may need to impose coastlines on the domain. 
# However, I'm trying really hard not to do that, since it's messy. 

# For reference: 
# Params=collections.namedtuple("Params",['strain_method','input_file','map_range','coord_box','coord_box_data',num_years','max_sigma','grid_inc','outdir','gmtfile']);
# Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);


import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.path as path
import sys
from scipy import interpolate
import strain_tensor_toolbox


def compute(myVelfield, MyParams):
	print("Computing strain via Numpy Spline method.");

	interp_kind='quintic';
	
	# Scipy returns a function that you can use on a new set of x,y pairs. 
	f_east = interpolate.interp2d(myVelfield.elon, myVelfield.nlat, myVelfield.e, kind=interp_kind);
	f_north = interpolate.interp2d(myVelfield.elon, myVelfield.nlat, myVelfield.n, kind=interp_kind);


	# The new interpolation grid: a new set of points with some chosen spacing
	xarray=np.arange(MyParams.coord_box[0],MyParams.coord_box[1],MyParams.grid_inc);
	yarray=np.arange(MyParams.coord_box[2],MyParams.coord_box[3],MyParams.grid_inc);
	[X,Y]=np.meshgrid(xarray,yarray);
	

	# Evaluate the linear or cubic interpolation function at new points
	new_east=np.zeros(np.shape(X));
	new_north=np.zeros(np.shape(X));
	for i in range(len(yarray)):
		for j in range(len(xarray)):
			new_east[i][j]=f_east(xarray[j],yarray[i]);  # only want to give the functions one point at a time. 
			new_north[i][j]=f_north(xarray[j],yarray[i]);

	# Grid increments
	lats = np.array([MyParams.map_range[2], MyParams.map_range[3]], dtype=float)
	typical_lat = np.mean(lats);
	xinc = MyParams.grid_inc * 111.000 * np.cos(np.deg2rad(typical_lat)); # in km (not degrees)
	yinc = MyParams.grid_inc * 111.000;  # in km (not degrees)

	# Computing the elements of the strain tensor from the 
	rot=np.zeros(np.shape(X));  # 2nd invariant of rotation rate tensor
	I2nd=np.zeros(np.shape(X));  # 2nd invariant of strain rate tensor
	max_shear=np.zeros(np.shape(X));  # max shear of strain rate tensor
	e1=np.zeros(np.shape(X));  # maximum principal strain (array of float)
	e2=np.zeros(np.shape(X));  # minimum principal strain (array of float)
	v00=np.zeros(np.shape(X));  # more complicated: eigenvectors (array of matrix 2x2)
	v01=np.zeros(np.shape(X));  # more complicated: eigenvectors (array of matrix 2x2)
	v10=np.zeros(np.shape(X));  # more complicated: eigenvectors (array of matrix 2x2)
	v11=np.zeros(np.shape(X));  # more complicated: eigenvectors (array of matrix 2x2)
	dilatation=np.zeros(np.shape(X));

	# the strain calculation
	for j in range(len(yarray)-1):
		for i in range(len(xarray)-1):
			up=new_east[j][i];
			vp=new_north[j][i];
			uq=new_east[j][i+1];
			vq=new_north[j][i+1];
			ur=new_east[j+1][i];
			vr=new_north[j+1][i];

			[dudx, dvdx, dudy, dvdy] = strain_tensor_toolbox.compute_displacement_gradients(up, vp, ur, vr, uq, vq, xinc, yinc);
			
			# The components that are easily computed
			# Units: nanostrain per year. 
			[exx, exy, eyy, rotation] = strain_tensor_toolbox.compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy);
			rot[j][i]=abs(rotation);

			# Compute a number of values based on tensor properties. 
			I2nd[j][i] = np.log10(np.abs(strain_tensor_toolbox.second_invariant(exx, exy, eyy)));
			[e11, e22, v1] = strain_tensor_toolbox.eigenvector_eigenvalue(exx, exy, eyy);
			e1[j][i]= e11;
			e2[j][i]= e22;
			v00[j][i]=v1[0][0];
			v10[j][i]=v1[1][0];
			v01[j][i]=v1[0][1];
			v11[j][i]=v1[1][1];
			max_shear[j][i] = (e11 - e22)/2;
			dilatation[j][i]= e11+e22; 

	print("Success computing strain via Numpy Spline method.\n");

	return [xarray, yarray, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation];


	# For debugging, make a plot of the velocities. 
	# color_levels=np.arange(-40,45,5);
	# plt.figure();
	# f,axarr=plt.subplots(1,2);
	# axarr[0].contourf(X,Y,new_east, color_levels, cmap='jet',vmin=-40, vmax=40);
	# axarr[0].quiver(myVelfield.elon, myVelfield.nlat, myVelfield.e, myVelfield.n);
	# axarr[0].set_title('East');
	# axarr[0].set_xlim([MyParams.coord_box[4],MyParams.coord_box[5]]);
	# axarr[0].set_ylim([MyParams.coord_box[6],MyParams.coord_box[7]]);	
	# h1=axarr[1].contourf(X,Y,new_north, color_levels, cmap='jet',vmin=-40, vmax=40);
	# axarr[1].quiver(myVelfield.elon, myVelfield.nlat, myVelfield.e, myVelfield.n);
	# axarr[1].set_xlim([MyParams.coord_box[4],MyParams.coord_box[5]]);
	# axarr[1].set_ylim([MyParams.coord_box[6],MyParams.coord_box[7]]);
	# axarr[1].set_title('North');
	# cbar=plt.colorbar(h1,ax=axarr[1]); cbar.set_label('mm/yr');
	# plt.savefig('test_velos.jpg');

	# sys.exit(0);
	# return;

