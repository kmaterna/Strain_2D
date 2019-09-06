import gps_input_functions as gpsin


def output_ETS(velfield, outname):
	outfile=open(outname,'w');
	for i in range(len(myVelfield.n)):
		outfile.write("%f %f %f %f %f %f 0 %s \n" % (myVelfield.elon[i], myVelfield.nlat[i], myVelfield.e[i], myVelfield.n[i], myVelfield.se[i], myVelfield.sn[i], myVelfield.name[i]));
	outfile.close();
	return;


[myVelfield] = gpsin.read_ETS_vel_file("../../Vel_Data/Other_vels/Bartlow_ETSvels.txt")
output_ETS(myVelfield, "ETS_velo.txt")

# [myVelfield] = gpsin.read_ETS_vel_file("../../Vel_Data/Other_vels/Bartlow_interETSvels.txt")
# output_ETS(myVelfield, "interETS_velo.txt")
