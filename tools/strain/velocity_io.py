
import collections

StationVel = collections.namedtuple('StationVel', ['elon', 'nlat', 'e', 'n', 'u', 'se', 'sn', 'su', 'name']);


def read_stationvels(input_file):
    # Reading a simple velocity format
    # Format: lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)
    print("Reading file %s " % input_file);
    myVelfield = [];
    ifile = open(input_file, 'r');
    for line in ifile:
        if len(line.split()) == 0:
            continue;
        if line.split()[0] == "#":
            continue;
        temp = line.split();
        lon = float(temp[0]);
        lat = float(temp[1]);
        VE = float(temp[2]);
        VN = float(temp[3]);
        VU = float(temp[4]);
        SE = float(temp[5]);
        SN = float(temp[6]);
        SU = float(temp[7]);
        if len(temp) > 8:
            name = temp[8];
        else:
            name = '';
        mystation = StationVel(elon=lon, nlat=lat, e=VE, n=VN, u=VU, se=SE, sn=SN, su=SU, name=name);
        myVelfield.append(mystation);
    ifile.close();
    return myVelfield;


def write_stationvels(myVelfield, output_file):
    # Writing a simple velocity format
    print("writing human-readable velfile in station-vel format, %s" % output_file);
    ofile = open(output_file, 'w');
    ofile.write(
        "# Format: lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)\n");
    for station_vel in myVelfield:
        ofile.write("%f %f %f %f %f %f %f %f %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.u, station_vel.se,
            station_vel.sn, station_vel.su, station_vel.name));
    ofile.close();
    return;


def read_horiz_vels(filename):
    # An even simpler format for 2D data, no error ellipses
    print("reading file %s " % filename);
    myVelfield = [];
    ifile = open(filename, 'r');
    for line in ifile:
        if len(line.split()) == 0:
            continue;
        if line.split()[0] == "#":
            continue;
        temp = line.split();
        lon = float(temp[0]);
        lat = float(temp[1]);
        VE = float(temp[2]);
        VN = float(temp[3]);
        VU = 0;
        se = 0;
        sn = 0;
        su = 0;
        mystation = StationVel(elon=lon, nlat=lat, e=VE, n=VN, u=VU, se=se, sn=sn, su=su, name='');
        myVelfield.append(mystation);
    ifile.close();
    return myVelfield;


def write_simple_gmt_format(myVelfield, outfile):
    # Write function for the even simpler format
    ofile = open(outfile, 'w');
    for item in myVelfield:
        ofile.write("%f %f %f %f %f %f 0.0\n" % (item.elon, item.nlat, item.e, item.n, item.se, item.sn));
    ofile.close();
    return;
