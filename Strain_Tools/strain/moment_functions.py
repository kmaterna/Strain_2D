# Computing moment via Savage and Simpson, 1997

import argparse
import sys
import os
import numpy as np
from . import strain_tensor_toolbox, utilities

help_message = "Compute moment accumulation rate via method of Savage and Simpson (1997) "


def cmd_parser(cmdargs):
    """Simple command line parser for moment accumulation rate calculator. Returns a dictionary of params."""
    p = argparse.ArgumentParser(
          description="\n"+help_message+"\n")
    if len(cmdargs) < 2:
        print("\n"+help_message+"\n")
        print("For full help message, try --help")
        sys.exit(0)
    p.add_argument('--netcdf', type=str, required=True,
                   help='''filename of Netcdf, an output from strain_2D, required''')
    p.add_argument('--outfile', type=str, required=True,
                   help='''filename of desired output text file, required''')
    p.add_argument('--mu', type=float, default=30,
                   help='''shear modulus [GPa], [default 30]''')
    p.add_argument('--depth', type=float, default=10,
                   help='''seismogenic thickness [km] [default 10]''')
    p.add_argument('--landmask', type=str, default='landmask.grd',
                   help='''grdfile showing where is land [default landmask.grd]''')
    p.add_argument('-v', '--verbose', action='count', default=0,
                   help='''controls verbosity [default 0]''')
    p.add_argument('--use_landmask', type=bool, default=1,
                   help='''whether or not to apply landmask [0 or 1; default 1]''')
    config_default = {}
    p.set_defaults(**config_default)
    config = vars(p.parse_args())
    config['outdir'] = os.path.split(config['netcdf'])[0]
    return config


def moment_coordinator(MyParams):
    """
    Coordinates the calculation.

    :param MyParams: a dictionary
    """
    # Input, Compute, Output
    lons, lats, exx, exy, eyy = utilities.read_basic_fields_from_netcdf(MyParams["netcdf"])
    landmask = np.ones(np.shape(exx))
    if MyParams["use_landmask"]:
        landmask = utilities.make_gmt_landmask(np.array(lons), np.array(lats), MyParams["landmask"])
    Mo, moment_map = compute_moments_loop(lons, lats, exx, exy, eyy, landmask, MyParams["mu"], MyParams["depth"])
    write_Mo_outputs(MyParams, Mo, lons, lats, moment_map)
    return Mo


def get_savage_simpson_moment(exx, exy, eyy, mu_GPa, depth_km, area_km2):
    """
    Minimum moment accumulation rate associated with a surface strain rate tensor.
    as per Savage and Simpson 1997, Equation 22
    exx, exy, eyy assumed in units of nanostrain
    With strain in nanostrain and mu in GPa, we don't need to multiply and then divide by 1e9
    """
    [e1, e2, _] = strain_tensor_toolbox.eigenvector_eigenvalue(exx, exy, eyy)
    depth_m = depth_km * 1000
    area_m2 = area_km2 * 1e6
    M0_min = 2 * mu_GPa * depth_m * area_m2 * np.max([np.abs(e1), np.abs(e2), np.abs(e1+e2)])
    return M0_min


def compute_moments_loop(lons, lats, exx, exy, eyy, landmask, mu, depth):
    print("Computing total moment rate of strain rate field from Savage and Simpson (1997).")
    Mo = 0
    xinc_km = (lons[1] - lons[0]) * (111.000*np.cos(np.deg2rad(lats[0])))
    yinc_km = (lats[1] - lats[0]) * 111.000
    area_km2 = xinc_km * yinc_km
    exx = np.multiply(exx, landmask)
    exy = np.multiply(exy, landmask)
    eyy = np.multiply(eyy, landmask)
    moment_map = np.zeros(np.shape(exx))
    for i in range(len(lats)):
        for j in range(len(lons)):
            if landmask[i][j] > 0:
                incremental_moment = get_savage_simpson_moment(exx[i][j], exy[i][j], eyy[i][j], mu, depth, area_km2)
                Mo = Mo + incremental_moment
                moment_map[i][j] = incremental_moment
    return Mo, moment_map


def write_Mo_outputs(MyParams, Mo, lons, lats, moment_map):
    print("Writing file %s " % MyParams["outfile"])
    print("Moment Accumulation Rate: %f e18 N-m / year" % (Mo/1e18))
    ofile = open(MyParams["outfile"], 'w')
    ofile.write("Mu: %f GPa\n" % (MyParams["mu"]))
    ofile.write("Depth: %f km\n" % (MyParams["depth"]))
    ofile.write("Infile: %s\n" % (MyParams["netcdf"]))
    ofile.write("Moment rate accumulation: %f e18 N-m / year\n" % (Mo/1e18))
    ofile.close()

    moment_outfile = os.path.join(MyParams['outdir'], 'moment_rate_map.txt')
    comment = ("# Moment rate accumulation rate, in N-m per year, with depth = " + str(MyParams["depth"]) +
               " km and shear modulus = " + str(MyParams["mu"]) + " GPa")
    write_text_grid_quantity(moment_outfile, lons, lats, moment_map, comment=comment)
    return


def write_text_grid_quantity(outfile, lons_1d, lats_1d, quantity, comment=""):
    """Write a text file representation of a grid"""
    ofile = open(outfile, 'w')
    print("Writing file %s " % outfile)
    ofile.write("# " + comment + "\n")
    for i in range(len(lats_1d)):
        for j in range(len(lons_1d)):
            ofile.write("%f %f %f\n" % (lons_1d[j], lats_1d[i], quantity[i][j]))
    ofile.close()
    return
