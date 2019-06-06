
import numpy as np
import collections
import datetime as dt 
import sys, os


Velfield = collections.namedtuple("Velfield",['name','nlat','elon','n','e','u','sn','se','su','first_epoch','last_epoch']);
Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);


def read_pbo_vel_file(infile):
# Meant for reading velocity files from the PBO/UNAVCO website. 
# Returns a Velfield collections object. 
	start=0;
	ifile=open(infile,'r');
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for line in ifile:
		if start==1:
			temp=line.split();
			name.append(temp[0]);
			nlat.append(float(temp[7]));
			elon_temp=float(temp[8]);
			if elon_temp>180:
				elon_temp=elon_temp-360.0;
			elon.append(elon_temp);
			n.append(float(temp[19])*1000.0);
			e.append(float(temp[20])*1000.0);
			u.append(float(temp[21])*1000.0);
			sn.append(float(temp[22])*1000.0);
			se.append(float(temp[23])*1000.0);
			su.append(float(temp[24])*1000.0);
			t1=temp[-2];
			t2=temp[-1];
			first_epoch.append(dt.datetime.strptime(t1[0:8],'%Y%m%d'));
			last_epoch.append(dt.datetime.strptime(t2[0:8],'%Y%m%d'));
		if "*" in line:
			start=1;
	ifile.close();
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];

def write_humanread_vel_file(myVelfield, outfile):
	ofile=open(outfile,'w');
	ofile.write("Format: lon(deg) lat(deg) e(mm) n(mm) u(mm) Se(mm) Sn(mm) Su(mm) first_date(yyyymmdd) last_date(yyyymmdd) name\n");
	for i in range(len(myVelfield.name)):
		first_epoch = dt.datetime.strftime(myVelfield.first_epoch[i],'%Y%m%d');
		last_epoch = dt.datetime.strftime(myVelfield.last_epoch[i],'%Y%m%d');
		ofile.write("%f %f %f %f %f %f %f %f %s %s %s\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.u[i], myVelfield.se[i], myVelfield.sn[i], myVelfield.su[i],first_epoch, last_epoch, myVelfield.name[i]) );
	ofile.close();
	return;

def write_gmt_velfile(myVelfield, outfile):
	ofile=open(outfile,'w');
	ofile.write("# Format: lon(deg) lat(deg) e(mm) n(mm) Se(mm) Sn(mm) 0 0 1 name\n");
	for i in range(len(myVelfield.name)):
		if myVelfield.sn[i] < 0.2:  # trying to make a clean dataset
			ofile.write("%f %f %f %f %f %f 0 0 1 %s\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i], myVelfield.name[i]) );
	ofile.close();
	return;


def read_unr_vel_file(infile):
# Meant for reading velocity files from the MAGNET/MIDAS website. 
# Returns a Velfield collections object. 	
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; 
	ifile=open(infile,'r');
	for line in ifile:
		temp=line.split();
		if temp[0]=="#":
			continue;
		else:
			name.append(temp[0]);
			e.append(float(temp[8])*1000.0);
			n.append(float(temp[9])*1000.0);
			u.append(float(temp[10])*1000.0);
			se.append(float(temp[11])*1000.0);
			sn.append(float(temp[12])*1000.0);
			su.append(float(temp[13])*1000.0);
	ifile.close();

	[elon,nlat]=get_coordinates_for_stations(name);
	[first_epoch, last_epoch] = get_start_times_for_stations(name);

	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield]; 



def clean_velfield(velfield, num_years, max_sigma, coord_box):
# Take the raw GPS velocities, and clean them up. 
# Remove data that's less than num_years long, 
# has formal uncertainties above max_sigma, 
# or is outside our box of intersest. 
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	for i in range(len(velfield.n)):
		if velfield.sn[i] > max_sigma:  
			continue;
		if velfield.se[i] > max_sigma:
			continue;
		deltatime=velfield.last_epoch[i]-velfield.first_epoch[i];
		if deltatime.days <= num_years*365.24:
			continue;
		if (velfield.elon[i]>coord_box[0] and velfield.elon[i]<coord_box[1] and velfield.nlat[i]>coord_box[2] and velfield.nlat[i]<coord_box[3]):
			#The station is within the box of interest. 
			name.append(velfield.name[i]);
			nlat.append(velfield.nlat[i]);
			elon.append(velfield.elon[i]);
			n.append(velfield.n[i]);
			sn.append(velfield.sn[i]);
			e.append(velfield.e[i]);
			se.append(velfield.se[i]);
			u.append(velfield.u[i]);
			su.append(velfield.su[i]);			
			first_epoch.append(velfield.first_epoch[i]);
			last_epoch.append(velfield.last_epoch[i]);
	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);
	return [myVelfield];


def remove_duplicates(velfield):
	name=[]; nlat=[]; elon=[]; n=[]; e=[]; u=[]; sn=[]; se=[]; su=[]; first_epoch=[]; last_epoch=[];
	
	for i in range(len(velfield.n)):
		is_duplicate = 0;
		for j in range(len(name)):
			if abs(nlat[j]-velfield.nlat[i])<0.0005 and abs(elon[j]-velfield.elon[i])<0.0005:
				# we found a duplicate measurement. 
				is_duplicate = 1;
				# Right now assuming all entries at the same lat/lon have the same velocity values. 

		if is_duplicate == 0:
			name.append(velfield.name[i]);
			nlat.append(velfield.nlat[i]);
			elon.append(velfield.elon[i]);
			n.append(velfield.n[i]);
			sn.append(velfield.sn[i]);
			e.append(velfield.e[i]);
			se.append(velfield.se[i]);
			u.append(velfield.u[i]);
			su.append(velfield.su[i]);			
			first_epoch.append(velfield.first_epoch[i]);
			last_epoch.append(velfield.last_epoch[i]);

	myVelfield = Velfield(name=name, nlat=nlat, elon=elon, n=n, e=e, u=u, sn=sn, se=sn, su=su, first_epoch=first_epoch, last_epoch=last_epoch);	
	return [myVelfield];



def get_coordinates_for_stations(station_names,coordinates_file="../../GPS_POS_DATA/UNR_DATA/UNR_coords_july2018.txt"):
	lon=[];
	lat=[];
	reference_names=[]; reference_lons=[]; reference_lats=[];

	# Read the file
	ifile=open(coordinates_file,'r');
	for line in ifile:
		temp=line.split();
		if temp[0]=="#":
			continue;
		reference_names.append(temp[0]);
		reference_lats.append(float(temp[1]));
		testlon=float(temp[2]);
		if testlon>180:
			testlon=testlon-360;
		reference_lons.append(testlon);
	ifile.close();

	# find the stations
	for i in range(len(station_names)):
		myindex=reference_names.index(station_names[i]);
		lon.append(reference_lons[myindex]);
		lat.append(reference_lats[myindex]);
		if myindex==[]:
			print("Error! Could not find coordinates for station %s " % station_names[i]);
			print("Returning [0,0]. ");
			lon.append(0.0);
			lat.append(0.0);

	return [lon,lat];


def get_start_times_for_stations(station_names,coordinates_file="../../GPS_POS_DATA/UNR_DATA/UNR_coords_july2018.txt"):
	# Meant for UNR stations
	end_time=[];
	start_time=[];
	reference_names=[]; reference_start_time=[]; reference_end_time=[];

	# Read the file
	ifile=open(coordinates_file,'r');
	for line in ifile:
		temp=line.split();
		if temp[0]=="#":
			continue;
		reference_names.append(temp[0]);
		reference_start_time.append(temp[7]);
		reference_end_time.append(temp[8]);
	ifile.close();

	# find the stations
	for i in range(len(station_names)):
		myindex=reference_names.index(station_names[i]);
		start_time.append(dt.datetime.strptime(reference_start_time[myindex],'%Y-%m-%d'));
		end_time.append(dt.datetime.strptime(reference_end_time[myindex],'%Y-%m-%d'));
		if myindex==[]:
			print("Error! Could not find startdate for station %s " % station_names[i]);
			print("Returning [0,0]. ");
			start_time.append(dt.datetime.strptime('2000-01-01','%Y-%m-%d'));
			end_time.append(dt.datetime.strptime('2000-01-01','%Y-%m-%d'));		
	return [start_time,end_time];




