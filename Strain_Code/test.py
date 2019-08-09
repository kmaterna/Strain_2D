import numpy as np 
import math as m
import ND_interp_strain
import configure_functions
import common_functions

def max_shortening_azimuth(e1, e2, v00, v01, v10, v11):
	az = np.zeros(e1.shape)
	n = [0, 1]
	for i in range(len(e1)):
		for j in range(len(e1[0])):
			if e1[i][j] > e2[i][j]:
				maxv = np.array([v00[i][j], v10[i][j]])
			else:
				maxv = np.array([v01[i][j], v11[i][j]])
			cos_theta = np.dot(maxv, n)
			degree = m.degrees(np.arccos(cos_theta))
			az[i][j] = degree
	return az


[MyParams] = configure_functions.configure("ND_interp")
[myVelfield] = common_functions.inputs(MyParams)
[xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = ND_interp_strain.compute(myVelfield, MyParams);
my_az = max_shortening_azimuth(e1, e2, v00, v01, v10, v11)

print(my_az[11])