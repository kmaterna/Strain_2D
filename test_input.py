import numpy as np
import gps_input_functions
import collections
import datetime as dt 
import sys, os


Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);

def read_ETS_vel_file(infile):
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; 
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		if len(temp)== 0:
			continue;
		elif temp[0]=="#" or temp[0]=="WSLR":
			continue;
		else:
			name.append(temp[0]);
			e.append(float(temp[1]));
			n.append(float(temp[2]));
			u.append(float(temp[3]));
			se = 0.1 * np.ones(len(name));
			sn = 0.1 * np.ones(len(name));
			su = 0.1 * np.ones(len(name));
	ifile.close();

	[elon,nlat] = gps_input_functions.get_coordinates_for_stations(name);
	[first_epoch, last_epoch] = gps_input_functions.get_start_times_for_stations(name);
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);

	return [myVelfield]; 


[myVelfield] = read_ETS_vel_file("Other_vels/Bartlow_interETSvels.txt")

print(len(myVelfield))
print(myVelfield[2])
