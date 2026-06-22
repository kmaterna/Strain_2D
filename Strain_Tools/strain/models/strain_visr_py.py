"""
Pure-Python port of the VISR strain-rate interpolation method.

This reproduces the strain-rate-interpolation path (``function = 3``) of the
original Fortran ``visr`` program by Zheng-Kang Shen, so that VISR can be run
without compiling and linking the Fortran source.  The algorithm follows:

    Shen, Z.-K., M. Wang, Y. Zeng, and F. Wang (2015), Strain determination
    using spatially discrete geodetic data, Bull. Seismol. Soc. Am., 105(4),
    2117-2127, doi:10.1785/0120140247.
    http://scec.ess.ucla.edu/~zshen/visr/visr.html

The implementation mirrors the Fortran ``visr_core``/``cmp_strain``/``llxy``
routines (see visr_test/visr.f) but vectorizes all per-station arithmetic with
NumPy so that the per-grid-point work is a handful of array operations plus a
single 6x6 weighted-least-squares solve.  Results are intended to match the
compiled Fortran to numerical precision.

Differences from the Fortran:
  * Grid points where interpolation fails (too few stations, insufficient
    azimuthal coverage, total weight below threshold, or point on a creep
    fault) are written as NaN rather than skipped.
  * Velocity uncertainties (Se, Sn) propagated from the solution covariance
    are returned, so this method populates the sigma layers.
"""

import numpy as np
from scipy.spatial import Voronoi, ConvexHull

from strain.models.strain_2d import Strain_2d
from .. import utilities

# --- Fixed constants, matching the Fortran source ---
DEG2RAD = np.pi / 180.0
EARTH_A_KM = 6378.137              # NAD83 semi-major axis [km]
FLATTENING = 1.0 / 298.2572        # NAD83 flattening
WT_AZ = 0.25                       # relative weight of mean azimuth (azimuth weighting)
NP_SITE = 6                        # neighbors used for circular-area substitute
CFA1 = 2.0                         # coefficient for circular-area weighting
CFA2 = 2.0                         # voronoi-area cap relative to circular area
CUTOFF_GAUSSIAN = 2.15             # selection cutoff (dr/tau) for gaussian weighting
CUTOFF_QUADRATIC = 10.0            # selection cutoff (dr/tau) for quadratic weighting


class visr_py(Strain_2d):
    """ Pure-Python VISR class for 2d strain rate, with general strain_2d behavior. """

    def __init__(self, params):
        super().__init__(params.inc, params.range_strain, params.range_data,
                         params.xdata, params.ydata, params.outdir)
        self._Name = 'visr_py'
        ms = params.method_specific

        # distance (temporal) weighting scheme
        self._distance_weighting = str(ms.get("distance_weighting", "gaussian")).lower()
        if self._distance_weighting not in ("gaussian", "quadratic"):
            raise ValueError("visr_py 'distance_weighting' must be 'gaussian' or 'quadratic'.")
        # spatial coverage weighting scheme
        self._spatial_weighting = str(ms.get("spatial_weighting", "voronoi")).lower()
        if self._spatial_weighting not in ("azimuth", "voronoi"):
            raise ValueError("visr_py 'spatial_weighting' must be 'azimuth' or 'voronoi'.")
        # smoothing-distance search: "min/max/inc" in km
        if "min_max_inc_smooth" not in ms:
            raise ValueError("visr_py requires 'min_max_inc_smooth' (e.g. '1/100/1').")
        smin, smax, sinc = [float(v) for v in ms["min_max_inc_smooth"].split('/')]
        # candidate taus: min, min+inc, ... strictly less than max (max is the give-up sentinel)
        self._taus = np.arange(int(smin), int(smax), int(sinc), dtype=float)
        if self._taus.size == 0:
            raise ValueError("visr_py 'min_max_inc_smooth' produced no candidate smoothing distances.")
        self._max_tau = float(smax)
        # weighting threshold Wt and velocity-uncertainty floor
        if "weighting_threshold" not in ms:
            raise ValueError("visr_py requires 'weighting_threshold'.")
        self._wt0 = float(ms["weighting_threshold"])
        self._unc_thresh = float(ms.get("uncertainty_threshold", 0.0))

        # creep faults (optional)
        self._num_creep = int(ms.get("num_creeping_faults", 0))
        self._creep_file = ms.get("creep_file", None)

        # selection cutoff depends on the distance weighting scheme
        self._cutoff = CUTOFF_GAUSSIAN if self._distance_weighting == "gaussian" else CUTOFF_QUADRATIC

    def compute(self, myVelfield):
        print("------------------------------\nComputing strain via pure-Python VISR method.")
        Ve, Vn, Se, Sn, rot, exx, exy, eyy = compute_visr_py(
            myVelfield, self._xdata, self._ydata, self._taus, self._cutoff,
            self._wt0, self._unc_thresh, self._distance_weighting, self._spatial_weighting,
            self._num_creep, self._creep_file,
        )
        # Observed/residual velocities within the strain bounding box
        velfield_within_box = utilities.filter_by_bounding_box(myVelfield, self._strain_range)
        model_velfield = utilities.create_model_velfield(self._xdata, self._ydata, Ve, Vn, velfield_within_box)
        residual_velfield = utilities.subtract_two_velfields(velfield_within_box, model_velfield)
        print("Success computing strain via pure-Python VISR method.\n")
        return [Ve, Vn, Se, Sn, rot, exx, exy, eyy, velfield_within_box, residual_velfield]


# ----------------- GEOMETRY -------------------------

def llxy(ref_lat, ref_lon, lat, lon):
    """
    Port of the Fortran ``llxy`` routine: project lat/lon (deg) onto a local
    cartesian frame (km) centered at (ref_lat, ref_lon), using the NAD83 radius
    of curvature at the reference latitude.

    :returns: (x, y) arrays in km, where x is ~east and y is ~north.
    """
    rlatc = ref_lat * DEG2RAD
    rlonc = ref_lon * DEG2RAD
    esq = 2.0 * FLATTENING - FLATTENING ** 2
    q = 1.0 - esq * np.sin(rlatc) ** 2
    r = EARTH_A_KM * np.sqrt(1.0 - esq) / q

    # transformation matrix rows we need (vp1, vp2)
    t11 = np.sin(rlatc) * np.cos(rlonc)
    t12 = np.sin(rlatc) * np.sin(rlonc)
    t13 = -np.cos(rlatc)
    t21 = -np.sin(rlonc)
    t22 = np.cos(rlonc)

    rlat = np.asarray(lat, dtype=float) * DEG2RAD
    rlon = np.asarray(lon, dtype=float) * DEG2RAD
    v1 = np.cos(rlat) * np.cos(rlon)
    v2 = np.cos(rlat) * np.sin(rlon)
    v3 = np.sin(rlat)
    vp1 = t11 * v1 + t12 * v2 + t13 * v3
    vp2 = t21 * v1 + t22 * v2  # t23 == 0
    x = r * vp2
    y = -r * vp1
    return x, y


# ----------------- SPATIAL WEIGHTING -------------------------

def compute_voronoi_areas(pos):
    """
    Voronoi-cell areas for each station, with a circular-area substitute for
    unbounded (hull) cells or cells larger than CFA2x the circular area.
    Mirrors get_voronoi_area_version + cmp_area1 in the Fortran.

    :param pos: (n, 2) array of cartesian station coordinates [km]
    :returns: (n,) array of areas [km^2]
    """
    n = pos.shape[0]
    # circular-area substitute: pi * (CFA1 * mean(NP_SITE nearest distances) / 2)^2
    diff = pos[:, None, :] - pos[None, :, :]
    dmat = np.sqrt(np.sum(diff ** 2, axis=2))
    sorted_d = np.sort(dmat, axis=1)
    nn = sorted_d[:, 1:1 + NP_SITE]            # skip self (column 0)
    carea = np.pi * (CFA1 * np.mean(nn, axis=1) / 2.0) ** 2

    areas = carea.copy()
    if n < 4:
        return areas  # Voronoi/ConvexHull undefined; fall back to circular areas
    vor = Voronoi(pos)
    hull = ConvexHull(pos)
    hull_vertices = set(hull.vertices.tolist())
    for i_station, i_region in enumerate(vor.point_region):
        if i_station in hull_vertices:
            continue
        region = vor.regions[i_region]
        if len(region) == 0 or -1 in region:
            continue
        try:
            actual = ConvexHull(vor.vertices[region]).volume   # 2D "volume" == area
        except Exception:
            continue
        if actual < CFA2 * carea[i_station]:
            areas[i_station] = actual
    return areas


def azimuth_weights(az_deg):
    """
    Azimuthal data-density weighting (Fortran id_wght == 1) for a selected set
    of stations, given their azimuths (deg) from the target point.
    """
    n = az_deg.size
    D = az_deg[None, :] - az_deg[:, None]      # D[i, j] = az[j] - az[i]
    D = np.where(D > 180.0, D - 360.0, D)
    D = np.where(D < -180.0, D + 360.0, D)
    np.fill_diagonal(D, np.nan)
    # smallest positive gap (default 180 if none); largest negative gap (default -180 if none)
    pos = np.where(D > 0.0, D, np.inf)
    neg = np.where(D < 0.0, D, -np.inf)
    daz1 = pos.min(axis=1)
    daz2 = neg.max(axis=1)
    daz1 = np.where(np.isinf(daz1), 180.0, daz1)
    daz2 = np.where(np.isinf(daz2), -180.0, daz2)
    azi_avrg = WT_AZ * 360.0 / n
    azi_tot = (1.0 + WT_AZ) * 360.0
    return (0.5 * (daz1 - daz2) + azi_avrg) * n / azi_tot


# ----------------- CREEP FAULTS -------------------------

def read_creep_faults(creep_file, num_creep, ref_lat, ref_lon):
    """
    Read creep-fault endpoints, project to local cartesian, and compute each
    fault's unit direction (dcs, dsn). Returns endpoint coords and directions.
    """
    alon, alat, blon, blat = [], [], [], []
    with open(creep_file, 'r') as f:
        count = 0
        for line in f:
            vals = line.split()
            if len(vals) < 4:
                continue
            a_lon, a_lat, b_lon, b_lat = (float(vals[0]), float(vals[1]),
                                          float(vals[2]), float(vals[3]))
            if a_lon > 180.0:
                a_lon -= 360.0
            if b_lon > 180.0:
                b_lon -= 360.0
            alon.append(a_lon); alat.append(a_lat); blon.append(b_lon); blat.append(b_lat)
            count += 1
            if count >= num_creep:
                break
    ax, ay = llxy(ref_lat, ref_lon, alat, alon)
    bx, by = llxy(ref_lat, ref_lon, blat, blon)
    dx = bx - ax
    dy = by - ay
    ds = np.sqrt(dx ** 2 + dy ** 2)
    dcs = dx / ds
    dsn = dy / ds
    return ax, ay, bx, by, dcs, dsn


# ----------------- CORE INTERPOLATION -------------------------

def compute_visr_py(myVelfield, xdata, ydata, taus, cutoff, wt0, unc_thresh,
                    distance_weighting, spatial_weighting, num_creep, creep_file):
    """
    Run the VISR strain-rate interpolation over the output grid (xdata, ydata).

    :returns: Ve, Vn, Se, Sn, rot, exx, exy, eyy  (each a 2D array shaped (ny, nx))
    """
    # --- extract station data directly (utilities.getVels has a known se/sn bug) ---
    lon = np.array([s.elon for s in myVelfield], dtype=float)
    lat = np.array([s.nlat for s in myVelfield], dtype=float)
    ux = np.array([s.e for s in myVelfield], dtype=float)
    uy = np.array([s.n for s in myVelfield], dtype=float)
    sx = np.array([s.se for s in myVelfield], dtype=float)
    sy = np.array([s.sn for s in myVelfield], dtype=float)
    # apply velocity-uncertainty floor (Fortran rsga reset)
    sx = np.maximum(sx, unc_thresh) if unc_thresh > 0 else sx
    sy = np.maximum(sy, unc_thresh) if unc_thresh > 0 else sy
    nstn = lon.size

    # --- network mass center and station cartesian coordinates ---
    xmean = lon.mean()
    ymean = lat.mean()
    sxc, syc = llxy(ymean, xmean, lat, lon)
    pos = np.column_stack([sxc, syc])

    # --- voronoi areas if needed ---
    area = compute_voronoi_areas(pos) if spatial_weighting == "voronoi" else None

    # --- creep-fault geometry (per station side is1 is point-independent) ---
    have_creep = num_creep > 0 and creep_file is not None
    if have_creep:
        ax, ay, bx, by, dcs, dsn = read_creep_faults(creep_file, num_creep, ymean, xmean)
        # signed side of each STATION relative to each fault (independent of target point)
        dxa_s = ax[None, :] - sxc[:, None]          # (nstn, nf)
        dya_s = ay[None, :] - syc[:, None]
        y1_stn = -(dya_s * dcs[None, :] - dxa_s * dsn[None, :])
        is1_stn = np.sign(y1_stn)                   # (nstn, nf)

    # --- output grids ---
    nx, ny = xdata.size, ydata.size
    Ve = np.full((ny, nx), np.nan)
    Vn = np.full((ny, nx), np.nan)
    Se = np.full((ny, nx), np.nan)
    Sn = np.full((ny, nx), np.nan)
    rot = np.full((ny, nx), np.nan)
    exx = np.full((ny, nx), np.nan)
    exy = np.full((ny, nx), np.nan)
    eyy = np.full((ny, nx), np.nan)

    # --- grid cartesian coordinates (about the same mass center) ---
    gxx, gyy = np.meshgrid(xdata, ydata)
    gxc, gyc = llxy(ymean, xmean, gyy.ravel(), gxx.ravel())
    gxc = gxc.reshape(ny, nx)
    gyc = gyc.reshape(ny, nx)

    print("Analyzing %d stations to compute %d grid points." % (nstn, nx * ny))

    for j in range(ny):
        for i in range(nx):
            px, py = gxc[j, i], gyc[j, i]
            dx = sxc - px
            dy = syc - py
            dr = np.sqrt(dx ** 2 + dy ** 2)
            ang = np.degrees(np.arctan2(dy, dx))

            # --- precompute per-fault sectors for this target point ---
            block_base = None
            on_fault = False
            if have_creep:
                ang1 = np.degrees(np.arctan2(ay - py, ax - px))
                ang2 = np.degrees(np.arctan2(by - py, bx - px))
                y1p = -((ay - py) * dcs - (ax - px) * dsn)
                same_ang = ang1 == ang2
                if np.any((~same_ang) & (y1p == 0.0)):
                    on_fault = True
                is0 = np.where(same_ang, 0.0, np.sign(y1p))
                strk1 = np.minimum(ang1, ang2)
                strk2 = np.maximum(ang1, ang2)
                dang = strk2 - strk1
                rotmask = dang > 180.0
                new_strk1 = np.where(rotmask, strk2, strk1)
                new_strk2 = np.where(rotmask, strk1 + 360.0, strk2)
                icrp = rotmask
                # station opposite side of fault from target
                opposite = is1_stn != is0[None, :]               # (nstn, nf)
                ang_adj = np.where((icrp[None, :]) & (ang[:, None] < 0.0),
                                   ang[:, None] + 360.0, ang[:, None])
                in_sector = (ang_adj >= new_strk1[None, :]) & (ang_adj <= new_strk2[None, :])
                block_base = np.any(opposite & in_sector, axis=1)   # (nstn,)

            if on_fault:
                continue  # indxx == 4

            # --- adaptive smoothing-distance search ---
            for tau in taus:
                sel = dr <= cutoff * tau
                if block_base is not None:
                    sel = sel & (~block_base)
                nslct = int(np.count_nonzero(sel))
                if nslct < 3:
                    continue                                       # indxx == 1

                az = ang[sel]
                # azimuthal coverage: BOTH the raw span and the [0,360)-wrapped span
                # must exceed 180 degrees (Fortran rejects if either is <= 180).
                az2 = np.where(az < 0.0, az + 360.0, az)
                if (az.max() - az.min()) <= 180.0 or (az2.max() - az2.min()) <= 180.0:
                    continue                                       # indxx == 2

                # spatial coverage weight
                if spatial_weighting == "azimuth":
                    az_w = np.where(az < 0.0, az + 360.0, az)
                    wght = azimuth_weights(az_w)
                else:
                    a_sel = area[sel]
                    wght = a_sel / a_sel.mean()

                drs = dr[sel]
                rt2 = (drs / tau) ** 2
                if distance_weighting == "gaussian":
                    wti = np.exp(-rt2) * wght
                else:
                    wti = wght / (rt2 + 1.0)
                wtt = 2.0 * wti.sum()
                if wtt < wt0:
                    continue                                       # indxx == 3

                # --- build weighted design matrix and solve normal equations ---
                sw = np.sqrt(wti)
                dxs, dys = dx[sel], dy[sel]
                sxs, sys = sx[sel], sy[sel]
                n = nslct
                A = np.zeros((2 * n, 6))
                # unknowns: [vx, vy, exx, exy, eyy, w]
                # ux rows
                A[:n, 0] = sw / sxs
                A[:n, 2] = sw * dxs / sxs
                A[:n, 3] = sw * dys / sxs
                A[:n, 5] = sw * dys / sxs
                # uy rows
                A[n:, 1] = sw / sys
                A[n:, 3] = sw * dxs / sys
                A[n:, 4] = sw * dys / sys
                A[n:, 5] = -sw * dxs / sys
                b = np.empty(2 * n)
                b[:n] = sw * ux[sel] / sxs
                b[n:] = sw * uy[sel] / sys

                N = A.T @ A
                rhs = A.T @ b
                try:
                    cov = np.linalg.inv(N)
                except np.linalg.LinAlgError:
                    continue
                m = cov @ rhs

                # store results (strain/rotation scaled by 1000 -> nanostrain/yr, nano-rad/yr)
                Ve[j, i] = m[0]
                Vn[j, i] = m[1]
                Se[j, i] = np.sqrt(cov[0, 0])
                Sn[j, i] = np.sqrt(cov[1, 1])
                exx[j, i] = m[2] * 1000.0
                exy[j, i] = m[3] * 1000.0
                eyy[j, i] = m[4] * 1000.0
                rot[j, i] = m[5] * 1000.0
                break  # first successful tau wins

    return Ve, Vn, Se, Sn, rot, exx, exy, eyy
