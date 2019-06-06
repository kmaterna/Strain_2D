

import numpy as np 
import strain_tensor_toolbox


# Testing framework:
# This example is solid body rotation. 
# Does solid body rotation come out? 

def test_rotation_tensor():
	up=0;
	vp=-1;
	ur=-1;
	vr=-1;
	uq=0;
	vq=0;
	xinc=1;
	yinc=1;

	[dudx, dvdx, dudy, dvdy] = strain_tensor_toolbox.compute_displacement_gradients(up, vp, ur, vr, uq, vq, xinc, yinc);
	[exx, exy, eyy, rotation] = strain_tensor_toolbox.compute_strain_components_from_dx(dudx, dvdx, dudy, dvdy);

	print("exx: %f" % exx); 
	print("exy: %f" % exy); 
	print("eyy: %f" % eyy); 
	print("rotation: %f" % rotation); 
	return;

if __name__=="__main__":
	test_rotation_tensor();