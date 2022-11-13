import collections
import glob
import os

import xarray as xr

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
        temp = line.split();
        try:
            lon = float(temp[0]);
        except ValueError:
            continue
        lat, VE, VN, VU = float(temp[1]), float(temp[2]), float(temp[3]), float(temp[4]);
        SE, SN, SU = float(temp[5]), float(temp[6]), float(temp[7]);
        if len(temp) > 8:
            name = temp[8];
        else:
            name = '';
        mystation = StationVel(elon=lon, nlat=lat, e=VE, n=VN, u=VU, se=SE, sn=SN, su=SU, name=name);
        myVelfield.append(mystation);
    ifile.close();
    return myVelfield;


def write_stationvels(myVelfield, output_file, header=""):
    """
    Writing basic format for 2D velocity data

    :param myVelfield: Velfield object
    :param output_file: name of velocity file
    :param header: optional string that is written into the header line
    """
    print("writing human-readable velfile in station-vel format, %s" % output_file);
    ofile = open(output_file, 'w');
    ofile.write("# "+header+"\n");
    ofile.write("lon(deg) lat(deg) VE(mm) VN(mm) VU(mm) SE(mm) SN(mm) SU(mm) name(optional)\n");
    for station_vel in myVelfield:
        ofile.write("%f %f %f %f %f %f %f %f %s\n" % (
            station_vel.elon, station_vel.nlat, station_vel.e, station_vel.n, station_vel.u, station_vel.se,
            station_vel.sn, station_vel.su, station_vel.name));
    ofile.close();
    return;


def read_gmt_format(filename):
    """
    Reading simplest gmt format for 2D velocity data, with error ellipses

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
        try:
            _lon = float(temp[0]);
        except ValueError:
            continue
        lon, lat, VE, VN = float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3]);
        se, sn = float(temp[4]), float(temp[5]);
        myVelfield.append(StationVel(elon=lon, nlat=lat, e=VE, n=VN, u=0, se=se, sn=sn, su=0.1, name=''));
    ifile.close();
    return myVelfield;


def write_gmt_format(myVelfield, outfile):
    """
    Writing simplest gmt format for 2D velocity data

    :param myVelfield: Velfield object
    :param outfile: name of velocity file
    """
    print("writing vector output file %s " % outfile);
    ofile = open(outfile, 'w');
    ofile.write("lon(deg) lat(deg) VE(mm) VN(mm) SE(mm) SN(mm) Corr\n");
    for item in myVelfield:
        ofile.write("%f %f %f %f %f %f 0.0\n" % (item.elon, item.nlat, item.e, item.n, item.se, item.sn));
    ofile.close();
    return;


def write_multisegment_file(polygon_vertices, quantity, filename):
    """
    Write a quantity (e.g., gmt multisegment file) to color each polygon on a map

    :param polygon_vertices: list, a type of data structure
    :param quantity: z-values, numpy array
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


# --------- READ FUNCTION FOR MULTIPLE STRAIN NETCDFS ----------- #
def read_multiple_strain_netcdfs(MyParams, plot_type):
    """
    Get all the models (e.g. gpsgridder, geostats, huang, etc.) that have computed plot_type of 
    strain rate and return them as a single xarray Dataset

    Parameters
    ----------
    MyParams: dict - Parameter Dictionary containing strain rate methods/directories in a sub-dict
    plot_type: str - The type of strain rate quantity to return. Can be max_shear, dilatation, etc.
    
    Returns
    -------
    ds_new: xarray Dataset - A dataset containing the plot_type variable from each type of model
    """
    building_dict = {}
    ds = [];
    for key, value in MyParams.strain_dict.items():
        try:
            specific_filename = glob.glob(value + os.sep + '*' + "_strain.nc")[0]
        except:
            breakpoint()
        ds = xr.load_dataset(specific_filename)
        building_dict[key] = ds[plot_type];

    ds_new = xr.Dataset(building_dict, coords=ds.coords)
    return ds_new
