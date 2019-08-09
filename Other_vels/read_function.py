# Inter-ETS velocity files
# From Noel Bartlow
# Ex: Bartlow_interETSvels.txt
# Run through get_coords(). 
# Want: 
# E = east
# N = north
# U = up
# Ea1 = Annual 1st term.  Annual 2nd term.  Semi-annual 1st term. Semi-annual 2nd term. Ignore all of these. 


def read_interETS_vels(filename):
	names = np.genfromtxt(filename,skip_header=8, usecols=(0), dtype=('unicode') );
	E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2 = np.genfromtxt(filename,skip_header=8, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True );
	return [names, E, N, U, Ea1, Na1, Ua1, Ea2, Na2, Ua2, Es1, Ns1, Us1, Es2, Ns2, Us2];
