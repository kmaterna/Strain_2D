# A toolbox once we have strain tensor components. 
# We want to compute:
# Rotation vs. Volume Strain
# 2nd invariant
# Dilatation
# Max Shear Strain
# Eigenvectors and Eigenvalues


import numpy as np 


def second_invariant(exx, exy, eyy):
	e2nd = exx*eyy - exy*exy;
	return e2nd;

def eigenvector_eigenvalue(exx, exy, eyy):
	T = np.array([[exx, exy],[exy, eyy]]);  # the tensor
	w,v=np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
	return [w[0], w[1], v]; 

def max_shear_strain(exx, exy, eyy):
	T = np.array([[exx, exy],[exy, eyy]]);  # the tensor
	w,v=np.linalg.eig(T);  # The eigenvectors and eigenvalues (principal strains) of the strain rate tensor
	# w = eigenvalues
	# v = eigenvectors
	max_shear = (w[0]-w[1]) * 0.5;
	return max_shear;


def compute_displacement_gradients(up, vp, ur, vr, uq, vq, dx, dy):
	# up, vp describe velocity at a reference point P
	# R and Q are two other points: Q offset by dx in the x direction, and R offset by dy in the y direction. 
	# In practical usage, these are in mm/yr and km. 
	dudx = (uq-up)/dx;
	dvdx = (vq-vp)/dx;
	dudy = (ur-up)/dy;
	dvdy = (vr-vp)/dy;
	return [dudx, dvdx, dudy, dvdy];


def compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy):
	# Given a displacement tensor, compute the relevant parts of the strain and rotation tensors. 
	# Also converts to nanostrain per year.
	# Rot is the off-diagonal element of the rotation tensor
	# http://www.engr.colostate.edu/~thompson/hPage/CourseMat/Tutorials/Solid_Mechanics/rotations.pdf
	exx=dudx*1000;
	exy=(0.5 * (dvdx+dudy) )*1000;
	eyy=dvdy*1000;
	rot=(0.5*(dvdx-dudy));
	rot=rot*1000.0;
	return [exx, exy, eyy, rot];

	