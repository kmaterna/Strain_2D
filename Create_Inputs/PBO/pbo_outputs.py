import gps_input_functions as gpsin

[myVelfield] = gpsin.read_pbo_vel_file("../../Vel_Data/NAM08_pbovelfile_feb2018.vel")
# [myVelfield] = gpsin.blacklist(myVelfield, "../../Vel_Data/blacklist_stations.txt")
# [myVelfield] = gpsin.remove_duplicates(myVelfield)
# [myVelfield] = gpsin.clean_velfield(myVelfield, 2, .95, [-125, -121, 37, 42])

def output_unr(velfield, outdir):
	outfile=open(outdir+"pbo_velo.txt",'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0.0\n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i]) );
	outfile.close();
	return;

output_unr(myVelfield, "../../Vel_Data/")