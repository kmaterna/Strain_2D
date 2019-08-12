import gps_input_functions as gpsin

# [myVelfield] = gpsin.blacklist("../Example_data/blacklist_stations.txt", myVelfield)


def output_unr(velfield, outdir, outname):
	outfile=open(outdir+outname,'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0 %s \n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i], myVelfield.name[i]));
	outfile.close();
	return;


[myVelfield] = gpsin.read_unr_vel_file("../../Example_data/midas.NA12.txt")
output_unr(myVelfield, "../../Example_data/", "unr_velo.txt")

[myVelfield] = gpsin.clean_velfield(myVelfield, 2, .95, [-125, -121, 37, 42])
[myVelfield] = gpsin.blacklist("../../Example_data/blacklist_stations.txt", myVelfield)
output_unr(myVelfield, "../../Example_data/", "unr_velo_clean.txt")
