
import numpy as np

# Math regarding the development of 2D misfit and chi2 metrics of data fitting
# Depends on lists of StationVels


def compute_misfits(resid_vels, obs_vels, outlier_tolerance=10.0):
    misfit_total, chi2_total = [], [];
    for resid, obs in zip(resid_vels, obs_vels):
        if np.abs(resid.e) > outlier_tolerance or np.abs(resid.n) > outlier_tolerance:
            continue;
        misfit_total.append(abs(resid.e));
        misfit_total.append(abs(resid.n));
        chi2_total.append(np.square(resid.e) / np.square(obs.se));
        chi2_total.append(np.square(resid.n) / np.square(obs.sn));
    print("Median absolute deviation:", np.median(misfit_total))
    print("Median chi2:", np.median(chi2_total));
    return misfit_total, chi2_total;
