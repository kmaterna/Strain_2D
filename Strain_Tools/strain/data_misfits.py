
import numpy as np
from . import velocity_io

# Math regarding the development of 2D misfit and chi2 metrics of data fitting


def compute_misfits(resid_vels, obs_vels, outlier_tolerance=10.0):
    """
    Compute the misfits from a vector of residual velocities, observed velocities

    :param resid_vels: list of station_vels
    :param obs_vels: list of station_vels
    :param outlier_tolerance: float, default 10 mm/yr
    :returns: list of misfits [mm/yr], list of chi2 [unitless]
    """
    misfit_total, chi2_total = [], []
    outlier_count = 0
    for resid, obs in zip(resid_vels, obs_vels):
        if np.abs(resid.e) > outlier_tolerance or np.abs(resid.n) > outlier_tolerance:
            outlier_count += 1
            continue
        misfit_total.append(abs(resid.e))
        misfit_total.append(abs(resid.n))
        chi2_total.append(np.square(resid.e) / np.square(obs.se))
        chi2_total.append(np.square(resid.n) / np.square(obs.sn))
    print("Percentage of residuals outside tolerance: %f" % (100*outlier_count/len(resid_vels)))
    print("Median absolute deviation:", np.median(misfit_total))
    print("Median chi2:", np.round(np.median(chi2_total), 5))
    normalized_chi2 = np.round(0.5 * np.sum(chi2_total)/np.sum(~np.isnan(np.array(chi2_total))),5)
    print("Normalized Chi^2:", normalized_chi2)
    return misfit_total, chi2_total, normalized_chi2


def write_misfits_to_file(misfits, chi2, chi2n, outfile):
    with open(outfile, 'a') as ofile:
        ofile.write("\n")
        ofile.write("Median absolute deviation: %.5f mm/yr\n" % (np.median(misfits)))
        ofile.write("Median chi2: %.5f\n" % (np.median(chi2)))
        ofile.write("Normalized chi2: %.5f\n" % (np.median(chi2n)))
    return


def misfits_coordinator(params):
    """ A driver for the data-fitting misfit computation. Operates on dictionary with file i/o options in its fields"""
    obs_vels = velocity_io.read_stationvels(params["obs_velfile"])
    resid_vels = velocity_io.read_stationvels(params["resid_velfile"])
    misfit_total, chi2_total, normalized_chi2 = compute_misfits(resid_vels, obs_vels)
    write_misfits_to_file(misfit_total, chi2_total, normalized_chi2, params["outfile"])
    return
