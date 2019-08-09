# Use the interpolation scheme of Zheng-Kang Shen et al.
# Shen, Z.-K., M. Wang, Y. Zeng, and F. Wang, Strain determination using spatially discrete geodetic data, 
# Bull. Seismol. Soc. Am., 105(4), 2117-2127, doi: 10.1785/0120140247, 2015.
# http://scec.ess.ucla.edu/~zshen/visr/visr.html

# The fortran files must be compiled and linked like this: 
# gfortran -c voronoi_area_version.f90
# gfortran visr.f voronoi_area_version.o -o visr.exe


import numpy as np 
import subprocess
import sys
import strain_tensor_toolbox

# Reference for velocity field objects:
# Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);

def compute(myVelfield, MyParams):

	print("Computing strain via Visr method.");

	strain_config_file='Strain_Code/visr/visr_strain.drv';
	strain_data_file='Strain_Code/visr/vec';  # THIS CAN ONLY BE 20 CHARACTERS LONG!!!
	strain_output_file='Strain_Code/visr/str';  # THIS CAN ONLY BE 20 CHARACTERS LONG!!!
	write_fortran_config_file(strain_config_file, strain_data_file, strain_output_file, MyParams);
	write_fortran_data_file(strain_data_file, myVelfield);
	call_fortran_compute(strain_config_file);
	# We convert that text file into grids, which we will write as GMT grd files. 
	[xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation] = make_output_grids_from_strain_out(strain_output_file);

	print("Success computing strain via Visr method.\n");
	return [xdata, ydata, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation];


def write_fortran_config_file(strain_config_file, strain_data_file, strain_output_file, MyParams):
	# The config file will have the following components. 
	"""
	visr/visr_drive_strain.drv contains: 
	visr/velh.cmm4                             ! Station coordinate and velocity solution file
	visr/strain.out                            ! Strain rate output file
	1                                          ! distance weighting scheme: 1=gaussian, 2=quadratic
	2                                          ! spatial weighting scheme: 1=azimuth, 2=voronoi area
	1 100 1                                    ! minimum, maximum, and incremental spatial smoothing constants (km)
	24                                         ! weighting threshold Wt
	0.5                                        ! uncertainty threshold for reset
	3                                          ! function: 1=velocity compatibility checking; 2=velocity interpolation; 3=strain rate interpolation
	-122.5 -114.0 32.0 37.5 0.04 0.04          ! Lon_min, Lon_max, Lat_min, Lat_max, dLon, dLat
	0                                          ! number of creep faults
	crp.dat                                    ! creep fault data file
	"""
	ofile=open(strain_config_file,'w');
	ofile.write(strain_data_file+'                              ! Station coordinate and velocity solution file\n');
	ofile.write(strain_output_file+'                             ! Strain rate output file\n');
	ofile.write('1                                          ! distance weighting scheme: 1=gaussian, 2=quadratic\n');
	ofile.write('2                                          ! spatial weighting scheme: 1=azimuth, 2=voronoi area\n');
	ofile.write('1 100 1                                    ! minimum, maximum, and incremental spatial smoothing constants (km)\n');
	ofile.write('2                                          ! weighting threshold Wt\n');
	ofile.write('0.05                                       ! uncertainty threshold for reset\n');
	ofile.write('3                                          ! function: 1=velocity compatibility checking; 2=velocity interpolation; 3=strain rate interpolation\n');
	ofile.write(str(MyParams.coord_box[0])+' '+str(MyParams.coord_box[1])+' '+str(MyParams.coord_box[2])+' '+str(MyParams.coord_box[3])+' '+str(MyParams.grid_inc)+' '+str(MyParams.grid_inc)+'                ! Lon_min, Lon_max, Lat_min, Lat_max, dLon, dLat\n');
	ofile.write('0                                          ! number of creep faults\n');
	ofile.write('crp.dat                                    ! creep fault data file\n');
	ofile.close();
	return;


def write_fortran_data_file(data_file, Velfield):
	ofile=open(data_file,'w');
	# 35    format(a8,2f10.4,2(f7.2,f5.2),f7.3)
	# 0102_GPS -119.2642   34.5655 -29.02 0.79  22.96 0.73  0.082     4   7.2  1994.4  
	for i in range(len(Velfield.name)):
		ofile.write(Velfield.name[i]+"_GPS ");
		ofile.write("%9.4f %9.4f " % (Velfield.elon[i], Velfield.nlat[i]) );
		ofile.write("%6.2f %4.2f %6.2f %4.2f %6.3f " % (Velfield.e[i], Velfield.se[i], Velfield.n[i], Velfield.sn[i], 0.001) );
		ofile.write("    5   2.1  2005.0\n");
	ofile.close();
	return;


def call_fortran_compute(config_file):
	# Here we will call the strain compute function, using visr's fortran code. 
	# It will output a large text file. 	
	print("Calling visr.exe fortran code to compute strain. ");
	subprocess.call('Strain_Code/visr/visr.exe < '+config_file, shell=True);
	return;



def make_output_grids_from_strain_out(infile):
	ifile=open(infile,'r');
	x=[]; y=[]; rotation=[]; exx=[]; exy=[]; eyy=[];
	for line in ifile:
		temp=line.split();
		if 'index' in temp or 'longitude' in temp or 'deg' in temp:
			continue;
		else:
			x.append(float(temp[0]));
			y.append(float(temp[1]));
			rotation.append(float(temp[7]));
			exx.append(float(temp[9]));
			exy.append(float(temp[11]));
			eyy.append(float(temp[13]));
	ifile.close();
	ax1=set(x);
	ax2=set(y);
	xlen=len(ax1);
	ylen=len(ax2);

	xaxis=sorted(ax1);
	yaxis=sorted(ax2);

	if xlen==0 and ylen==0:
		print("ERROR! No valid strains have been computed. Try again.")
		sys.exit(0);
	
	# Loop through x and y lists, find the index of the coordinates in the xaxis and yaxis sets, 
	# Then place them into the 2d arrays. 
	# Then go compute I2nd, eigenvectors and eigenvalues. 

	I2nd=np.zeros((ylen, xlen)); max_shear=np.zeros((ylen, xlen)); rot=np.zeros((ylen, xlen)); e1=np.zeros((ylen, xlen)); e2=np.zeros((ylen, xlen)); dilatation=np.zeros((ylen, xlen))
	v00=np.zeros((ylen, xlen)); v01=np.zeros((ylen, xlen)); v10=np.zeros((ylen, xlen)); v11=np.zeros((ylen, xlen)); 
	print(np.shape(xaxis));
	print(np.shape(I2nd))

	for i in range(len(x)):
		xindex=xaxis.index(x[i]);
		yindex=yaxis.index(y[i]);
		rot[yindex][xindex]=rotation[i];
		I2nd[yindex][xindex] = np.log10(np.abs(strain_tensor_toolbox.second_invariant(exx[i], exy[i], eyy[i])));
		[e11, e22, v1] = strain_tensor_toolbox.eigenvector_eigenvalue(exx[i], exy[i], eyy[i]);
		e1[yindex][xindex]= e11;
		e2[yindex][xindex]= e22;
		v00[yindex][xindex]=v1[0][0];
		v10[yindex][xindex]=v1[1][0];
		v01[yindex][xindex]=v1[0][1];
		v11[yindex][xindex]=v1[1][1];
		dilatation[yindex][xindex]=e11+e22;
		max_shear[yindex][xindex] = (e11 - e22)/2;


	return [xaxis, yaxis, I2nd, max_shear, rot, e1, e2, v00, v01, v10, v11, dilatation];


