from Tectonic_Utils.read_write import netcdf_read_write
import collections
import os
from . import utilities

StationVel = collections.namedtuple('StationVel', ['elon', 'nlat', 'e', 'n', 'u', 'se', 'sn', 'su', 'name']);


def read_stationvels(input_file):
    """
    Reading basic format for 2D velocity data.
    Format: lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)

    :param input_file: name of velocity file
    :returns: list of stationvels
    """
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
    """
    Writing basic format for 2D velocity data

    :param myVelfield: Velfield object
    :param output_file: name of velocity file
    """
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


def read_gmt_format(filename):
    """
    Reading simplest gmt format for 2D velocity data, usually without error ellipses

    :param filename: name of velocity file
    :returns: list of stationvels
    """
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
        myVelfield.append(StationVel(elon=lon, nlat=lat, e=VE, n=VN, u=0, se=0, sn=0, su=0, name=''));
    ifile.close();
    return myVelfield;


def write_gmt_format(myVelfield, outfile):
    """
    Writing simplest gmt format for 2D velocity data

    :param myVelfield: Velfield object
    :param outfile: name of velocity file
    """
    print("writing velocity output file %s " % outfile);
    ofile = open(outfile, 'w');
    ofile.write("# Format: lon(deg) lat(deg) VE(mm) VN(mm) SE(mm) SN(mm) Corr\n");
    for item in myVelfield:
        ofile.write("%f %f %f %f %f %f 0.0\n" % (item.elon, item.nlat, item.e, item.n, item.se, item.sn));
    ofile.close();
    return;


def write_multisegment_file(polygon_vertices, quantity, filename):
    """
    Write a quantity (e.g., gmt multisegment file) to color each polygon on a map

    :param polygon_vertices: list, a type of data structure
    :param quantity: z-values
    :type quantity: list
    :param filename: output file name
    :type filename: string
    """
    print("Writing output file %s " % filename);
    ofile = open(filename, 'w');
    for i in range(len(quantity)):
        # Write the value associated with the triangle
        ofile.write("> -Z" + str(quantity[i]) + "\n");
        ofile.write(str(polygon_vertices[i, 0, 0]) + " " + str(polygon_vertices[i, 0, 1]) + "\n");
        ofile.write(str(polygon_vertices[i, 1, 0]) + " " + str(polygon_vertices[i, 1, 1]) + "\n");
        ofile.write(str(polygon_vertices[i, 2, 0]) + " " + str(polygon_vertices[i, 2, 1]) + "\n");
    ofile.close();
    return;


# --------- READ FUNCTION ----------- #
def read_multiple_strain_files(MyParams, filename):
    """
    Read strain quantities of `filename` into a dictionary
    Each dictionary key is a strain method
    Each dictionary value is a data structure: [lon, lat, value]
    lon : list of floats
    lat : list of floats
    value : 2D array of floats
    We also guarantee the mutual co-registration of the dictionary elements
    """
    strain_values_dict = {};
    for method in MyParams.strain_dict.keys():
        specific_filename = MyParams.strain_dict[method]+"/"+filename
        assert(os.path.isfile(specific_filename)), FileNotFoundError("Cannot find file " + specific_filename);
        if os.path.isfile(specific_filename):
            [lon, lat, val] = netcdf_read_write.read_any_grd(specific_filename);
            strain_values_dict[method] = [lon, lat, val];
    utilities.check_coregistered_grids(MyParams.range_strain, MyParams.inc, strain_values_dict);
    utilities.check_coregistered_shapes(strain_values_dict);
    return strain_values_dict;
