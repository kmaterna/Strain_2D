"""
May 2018
Take a set of velocities, establish delaunay triangles, 
solve a linear inversion problem for the components of the velocity gradient tensor
at the centroid of each triangle. 
The strain rate tensor and the rotation tensor can be readily computed 
from the symmetric and anti-symmetric parts of the velocity gradient tensor. 
Plot the outputs. 

Following a technique learned in Brad Hagar's geodynamics class, and 
modeled off of advice from 2007 Journal of Geodynamcis paper:
ftp://ftp.ingv.it/pub/salvatore.barba/RevEu/Cai_StrainBIFROST_2007.pdf
"""

import numpy as np
from scipy.spatial import Delaunay
from numpy.linalg import inv, det
import strain_tensor_toolbox
import collections
import datetime as dt
import output_manager

Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);


def read_pbo_vel_file(infile):
# Old format
	ifile=open(infile,'r');
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for line in ifile:
		temp=line.split();
		name.append('zzzz');
		nlat.append(float(temp[1]));
		elon_temp=float(temp[0]);
		if elon_temp>180:
			elon_temp=elon_temp-360.0;
		elon.append(elon_temp);
		n.append(float(temp[3]));
		e.append(float(temp[2]));
		u.append(0);
		sn.append(float(temp[4]));
		se.append(float(temp[5]));
		su.append(0);
		first_epoch.append(dt.datetime.strptime("20190202"[0:8],'%Y%m%d'));
		last_epoch.append(dt.datetime.strptime("20300101",'%Y%m%d'));
	ifile.close();
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];


def output_oldstyle(xcentroid, ycentroid, polygon_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation,
				   azimuth, myVelfield, MyParams):
	print("------------------------------\nWriting 1d outputs:");
	output_manager.write_multisegment_file(polygon_vertices, rot, MyParams.outdir + "rotation.txt");
	output_manager.write_multisegment_file(polygon_vertices, I2nd, MyParams.outdir + "I2nd.txt");
	output_manager.write_multisegment_file(polygon_vertices, dilatation, MyParams.outdir + "Dilatation.txt");
	output_manager.write_multisegment_file(polygon_vertices, max_shear, MyParams.outdir + "max_shear.txt");
	output_manager.write_multisegment_file(polygon_vertices, azimuth, MyParams.outdir + "azimuth.txt");
	# gps_io_functions.write_humanread_vel_file(myVelfield, MyParams.outdir + "tempgps.txt");

	positive_file = open(MyParams.outdir + "positive_eigs.txt", 'w');
	negative_file = open(MyParams.outdir + "negative_eigs.txt", 'w');
	for i in range(len(I2nd)):
		# Write the eigenvectors and eigenvalues
		output_manager.write_single_eigenvector(positive_file, negative_file, e1[i], v00[i], v10[i], xcentroid[i], ycentroid[i]);
		output_manager.write_single_eigenvector(positive_file, negative_file, e2[i], v01[i], v11[i], xcentroid[i], ycentroid[i]);
	positive_file.close();
	negative_file.close();

	print("Max I2: %f " % (max(I2nd)));
	print("Max rot: %f " % (max(rot)));
	print("Min rot: %f " % (min(rot)));
	return;


def compute_test(Inputs, MyParams):
	[myVelfield] = read_pbo_vel_file(MyParams.input_file);
	[xcentroid, ycentroid, triangle_vertices, rot, e1, e2, v00, v01, v10, v11] = compute(myVelfield, MyParams);
	[I2nd, max_shear, dilatation, azimuth] = strain_tensor_toolbox.compute_derived_quantities(e1, e2, v00, v01, v10, v11);
	output_oldstyle(xcentroid, ycentroid, triangle_vertices, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation, azimuth, myVelfield, MyParams);
	return [xcentroid, ycentroid, triangle_vertices, rot, e1, e2, v00, v01, v10, v11];


# ----------------- COMPUTE -------------------------

def compute(myVelfield, MyParams):

	print("Computing strain via delaunay method.");
	z = np.array([myVelfield.elon,myVelfield.nlat]);
	z = z.T;
	tri=Delaunay(z);

	triangle_vertices = z[tri.simplices];
	trishape = np.shape(triangle_vertices);  # 516 x 3 x 2, for example

	# We are going to solve for the velocity gradient tensor at the centroid of each triangle. 
	centroids=[];
	for i in range(trishape[0]):
		xcor_mean = np.mean([triangle_vertices[i,0,0],triangle_vertices[i,1,0],triangle_vertices[i,2,0]]);
		ycor_mean = np.mean([triangle_vertices[i,0,1],triangle_vertices[i,1,1],triangle_vertices[i,2,1]]);
		centroids.append([xcor_mean,ycor_mean]);
	xcentroid=[x[0] for x in centroids];
	ycentroid=[x[1] for x in centroids];


	# Initialize arrays.
	rot=[];
	e1=[]; # eigenvalues
	e2=[];
	v00=[];  # eigenvectors
	v01=[];
	v10=[];
	v11=[];


	# for each triangle:
	for i in range(trishape[0]):

		# Get the velocities of each vertex (VE1, VN1, VE2, VN2, VE3, VN3)
		# Get velocities for Vertex 1 (triangle_vertices[i,0,0] and triangle_vertices[i,0,1])
		xindex1 = np.where(myVelfield.elon==triangle_vertices[i,0,0])
		yindex1 = np.where(myVelfield.nlat==triangle_vertices[i,0,1])
		index1=np.intersect1d(xindex1,yindex1);
		xindex2 = np.where(myVelfield.elon==triangle_vertices[i,1,0])
		yindex2 = np.where(myVelfield.nlat==triangle_vertices[i,1,1])
		index2=np.intersect1d(xindex2,yindex2);
		xindex3 = np.where(myVelfield.elon==triangle_vertices[i,2,0])
		yindex3 = np.where(myVelfield.nlat==triangle_vertices[i,2,1])
		index3=np.intersect1d(xindex3,yindex3);

		VE1=myVelfield.e[index1[0]]; VN1=myVelfield.n[index1[0]];
		VE2=myVelfield.e[index2[0]]; VN2=myVelfield.n[index2[0]];
		VE3=myVelfield.e[index3[0]]; VN3=myVelfield.n[index3[0]];
		obs_vel = np.array([[VE1],[VN1],[VE2],[VN2],[VE3],[VN3]]);

		# Get the distance between centroid and vertex (in km)
		dE1 = (triangle_vertices[i,0,0]-xcentroid[i])*111.0*np.cos(np.deg2rad(ycentroid[i]));
		dE2 = (triangle_vertices[i,1,0]-xcentroid[i])*111.0*np.cos(np.deg2rad(ycentroid[i]));
		dE3 = (triangle_vertices[i,2,0]-xcentroid[i])*111.0*np.cos(np.deg2rad(ycentroid[i]));
		dN1 = (triangle_vertices[i,0,1]-ycentroid[i])*111.0;
		dN2 = (triangle_vertices[i,1,1]-ycentroid[i])*111.0;
		dN3 = (triangle_vertices[i,2,1]-ycentroid[i])*111.0;

		Design_Matrix = np.array([[1,0,dE1,dN1,0,0],[0,1,0,0,dE1,dN1],[1,0,dE2,dN2,0,0],[0,1,0,0,dE2,dN2],[1,0,dE3,dN3,0,0],[0,1,0,0,dE3,dN3]]);

		# Invert to get the components of the velocity gradient tensor. 
		DMinv = inv(Design_Matrix);
		vel_grad = np.dot(DMinv, obs_vel);  # this is the money step. 
		VE_centroid=vel_grad[0][0];
		VN_centroid=vel_grad[1][0];
		dVEdE=vel_grad[2][0];
		dVEdN=vel_grad[3][0];
		dVNdE=vel_grad[4][0];
		dVNdN=vel_grad[5][0];


		# The components that are easily computed
		[exx, exy, eyy, rotation] = strain_tensor_toolbox.compute_strain_components_from_dx(dVEdE, dVNdE, dVEdN, dVNdN);

		# # Compute a number of values based on tensor properties. 
		[e11, e22, v] = strain_tensor_toolbox.eigenvector_eigenvalue(exx, exy, eyy);

		e1.append(e11);
		e2.append(e22);
		rot.append(abs(rotation));
		v00.append(v[0][0]);
		v10.append(v[1][0]);
		v01.append(v[0][1]);
		v11.append(v[1][1]);

	print("Success computing strain via delaunay flat-earth method.\n");


	return [xcentroid, ycentroid, triangle_vertices, rot, e1, e2, v00, v01, v10, v11];




