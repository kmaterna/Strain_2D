import numpy as np 
import collections
import read_velo as rv
import scipy.spatial as sp 
import matplotlib.pyplot as plt

filename = "..\\Example_Data\\fake_data.txt"
velo = rv.input(filename)

# makes arrays of coordinates into pairs
def coordinate_pairs(longitude, latitude):
	coords = []
	for i in range(len(latitude)):
		pair = [longitude[i], latitude[i]]
		coords.append(pair)
	return coords

# calls delaunay on coordinate pairs, returns delaunay indices and an array with the coordinates of each triangle  
def delaunay_to_coords(coordinates):
	pure_delaun = sp.Delaunay(coordinates)
	triangle_coords = coordinates[pure_delaun.simplices]
	return pure_delaun, triangle_coords

# finds the center of each triangle in degrees
def find_centroids(delaunay_matrix, coordinate_pairs):
	centroids = []
	pair = []
	for triangle in delaunay_matrix.simplices:
		vertices = coordinate_pairs[triangle]
		x = (vertices[0][0] + vertices[1][0] + vertices[2][0])/3
		y = (vertices[0][1] + vertices[1][1] + vertices[2][1])/3
		centroids.append([x, y])
	return centroids

# finds the distance from the centroid of eahc triangle to its vertices, in kilometers
def dist_to_centroids(triangle_matrix, centroids):
	dx_arr = []
	dy_arr = []
	for i in range(len(centroids)):
		triangle_dx = []
		triangle_dy = []
		for j in range(3):
			dx = abs(centroids[i][0] - triangle_matrix[i][j][0]) * 111.139 * np.cos(np.deg2rad(centroids[i][1]))
			triangle_dx.append(dx)
			dy = abs(centroids[i][1] - triangle_matrix[i][j][1]) * 111.139
			triangle_dy.append(dy)
		dx_arr.append(triangle_dx)
		dy_arr.append(triangle_dy)
	return dx_arr, dy_arr

# creates an array of the "d" arrays: the velocities of the vertices of each triangle, in mm/yr
def configure_d(dEdt_arr, dNdt_arr, pure_delaun):
	triangle_dEdt = dEdt_arr[pure_delaun.simplices]
	triangle_dNdt = dNdt_arr[pure_delaun.simplices]
	one_d = []
	all_d = []
	for i in range(len(triangle_dEdt)):
		one_d = (triangle_dEdt[i][0], triangle_dNdt[i][0], triangle_dEdt[i][1], triangle_dNdt[i][1], triangle_dEdt[i][2], triangle_dNdt[i][2])
		all_d.append(np.array(one_d)*1000)
	return all_d

# creates an array of "G" matrices from the distances from centroid to vertices of each triangle according to technique from Cai et al.
def configure_G(dx_arr, dy_arr):
	all_G = []
	for i in range(len(dx_arr)):
		col0 = (1, 0)*3
		col1 = (0, 1)*3
		col2 = (dx_arr[i][0], 0, dx_arr[i][1], 0, dx_arr[i][2], 0)
		col3 = (dy_arr[i][0], 0, dy_arr[i][1], 0, dy_arr[i][2], 0)
		col4 = (0, dx_arr[i][0], 0, dx_arr[i][1], 0, dx_arr[i][2])
		col5 = (0, dy_arr[i][0], 0, dy_arr[i][1], 0, dy_arr[i][2])
		G = np.column_stack([col0, col1, col2, col3, col4, col5])
		all_G.append(G)
	return all_G

# performs matrix inversion and multiplication to find an "m" array for each triangle, with the velocity gradients
def solve_for_m(all_G, all_d):
	all_m = []
	for i in range(len(all_G)):
		G = all_G[i]
		d = all_d[i]
		# GtG = np.matmul(np.transpose(G), G)
		# Gtd = np.matmul(np.transpose(G), d)
		# m = np.matmul(np.linalg.inv(GtG), Gtd)
		m = np.linalg.solve(G, d)
		all_m.append(m)
	return all_m

# takes the gradients for each triangle and produces a corresponding strain rate tensor in mm/Km*yr
def gradients_to_strain(all_m):
	all_T = []
	for m in all_m:
		row0 = (m[2], .5 * (m[3] + m[4]))
		row1 = (.5 * (m[3] + m[4]) , m[5])
		T = np.stack([row0, row1])
		all_T.append(T)
	return all_T

# takes a strain rate tensor and returns the eigenvectors and eigenvalues for each triangle
def strain_to_eigen(all_T):
	exx = []
	exy = []
	eyy = []
	e1 = []
	e2 = []
	v00 = []
	v10 = []
	v01 = []
	v11 = []
	for T in all_T:
		exx.append(T[0][0])
		exy.append(T[0][1])
		eyy.append(T[1][1])
		w, v = np.linalg.eig(T)
		e1.append(w[0])
		e2.append(w[1])
		v00.append(v[0][0]);
		v10.append(v[1][0]);
		v01.append(v[0][1]);
		v11.append(v[1][1]);
	return exx, exy, eyy, e1, e2, v00, v10, v01, v11

# computes the second invariant from the eigenvalues
# def second_invariant(e1, e2):
# 	all_inv = []
# 	for i in range(len(e1)):
# 		inv = (e1[i] ** 2 + e2[i] ** 2) **.5
# 		if inv == 0:
# 			all_inv.append(0)
# 		else:
# 			all_inv.append(np.log10(np.abs(inv)))
# 	return all_inv

# computes the second invariant from the eigenvectors
def second_invariant(exx, exy, eyy):
	all_I2 = []
	for i in range(len(exx)):
		I2 = (exx[i] * eyy[i]) - (exy[i] ** 2)
		if I2 == 0:
			all_I2.append(0)
		else:
			all_I2.append(np.log10(np.abs(I2)))
	return all_I2

# computes maximum shear strain from the eigenvalues
def max_shear(e1, e2):
	shears = []
	for i in range(len(e1)):
		shear = (e1[i] - e2[i]) / 2
		if shear == 0: 
			shears.append(shear)
		else:
			shears.append(np.log10(np.abs(shear)))
	return shears

# outputs a .txt with corresponding columns to be read by psxy in GMT
# def output_invariants(filename, invars, triangle_matrix):
# 	with open(filename, 'w') as f:
# 		for i in range(len(invars)):
# 			f.write("> -Z"+str(invars[i])+"\n");
# 			f.write(str(triangle_matrix[i,0,0])+" "+str(triangle_matrix[i,0,1])+"\n");
# 			f.write(str(triangle_matrix[i,1,0])+" "+str(triangle_matrix[i,1,1])+"\n");
# 			f.write(str(triangle_matrix[i,2,0])+" "+str(triangle_matrix[i,2,1])+"\n");
# 	return




def compute(myVelfield, myParams):
	coords = np.array(coordinate_pairs(myVelfield.elong, myVelfield.nlat))
	delaun, triangle_vertices = delaunay_to_coords(coords)
	centroids = find_centroids(delaun, coords)
	xcentroid=[x[0] for x in centroids]
	ycentroid=[x[1] for x in centroids]
	mydxs, mydys = dist_to_centroids(triangle_vertices, centroids)
	my_ds = configure_d(myVelfield.e, myVelfield.n, delaun)
	my_Gs = configure_G(mydxs, mydys)
	my_ms = solve_for_m(my_Gs, my_ds)
	my_Ts = gradients_to_strain(my_ms)
	exx, exy, eyy, e1, e2, v00, v10, v01, v11 = strain_to_eigen(my_Ts)
	I2 = second_invariant(exx, exy, eyy)
	max_shear = max_shear(e1, e2)

	return [xcentroid, ycentroid, triangle_vertices, I2, max_shear, e1, e2, v00, v10, v01, v11];

