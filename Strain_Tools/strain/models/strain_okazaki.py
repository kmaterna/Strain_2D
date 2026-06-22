"""
Strain rate via basis-function expansion with ABIC (Okazaki et al., 2021).

This is a pure-Python implementation of the global Bayesian method of:

    Okazaki, T., Fukahata, Y. & Nishimura, T. (2021). Consistent estimation of
    strain-rate fields from GNSS velocity data using basis function expansion
    with ABIC. Earth Planets Space 73, 153. doi:10.1186/s40623-021-01474-5

Each velocity component is expanded in 2D tensor-product cubic B-splines on a
regular grid (spacing ``grid_km``). The basis coefficients are estimated by a
smoothness-regularized least squares, where the smoothing hyperparameter is
selected objectively by minimizing Akaike's Bayesian Information Criterion
(ABIC). Velocity gradients (and hence strain rates) follow analytically from
the B-spline derivatives, and uncertainties are propagated from the posterior
model covariance.

This implementation is adapted from the author's reference code (provided in
models/230920_gnss_abic_okazaki/) but restructured for the Strain_2d API:
  * the optimize-then-estimate workflow runs in-memory in a single compute()
    (no pickle/file IO),
  * the observation matrix H is built vectorized over stations,
  * the roughness matrix R is built as a sum of Kronecker products of small
    banded 1D integral matrices (the reference builds it with an O(n_basis^2)
    Python loop),
  * the estimation grid is the Strain_2d output grid (xdata, ydata).
"""

import numpy as np
import scipy.sparse as sparse
from scipy.linalg import cho_factor, cho_solve

from strain.models.strain_2d import Strain_2d
from .. import utilities

# Okazaki's degree->km factor (1000/9 ~= 111.11 km/deg); longitude scaled by cos(lat0).
KM_PER_DEG = 1000.0 / 9.0


class okazaki(Strain_2d):
    """ Okazaki ABIC basis-function-expansion method for 2D strain rate. """

    def __init__(self, params):
        super().__init__(params.inc, params.range_strain, params.range_data,
                         params.xdata, params.ydata, params.outdir)
        self._Name = 'okazaki'
        ms = params.method_specific

        if "grid_km" not in ms:
            raise ValueError("okazaki requires 'grid_km' (B-spline grid spacing in km).")
        self._grid_km = float(ms["grid_km"])
        if self._grid_km <= 0:
            raise ValueError("okazaki 'grid_km' must be positive.")
        # optional: correlation length (km) of the data-error covariance (0 = white noise)
        self._corr_leng = float(ms.get("corr_leng", 0.0))
        # optional: extra rings of basis nodes added beyond the data extent (>=2 for full
        # cubic-spline support at the edges)
        self._node_margin = int(ms.get("node_margin", 2))
        if self._node_margin < 2:
            raise ValueError("okazaki 'node_margin' must be >= 2 for full edge support.")

    def compute(self, myVelfield):
        print("------------------------------\nComputing strain via Okazaki ABIC basis-function method.")
        Ve, Vn, Se, Sn, rot, exx, exy, eyy = compute_okazaki(
            myVelfield, self._xdata, self._ydata, self._grid_km,
            corr_leng=self._corr_leng, node_margin=self._node_margin)
        velfield_within_box = utilities.filter_by_bounding_box(myVelfield, self._strain_range)
        model_velfield = utilities.create_model_velfield(self._xdata, self._ydata, Ve, Vn, velfield_within_box)
        residual_velfield = utilities.subtract_two_velfields(velfield_within_box, model_velfield)
        print("Success computing strain via Okazaki method.\n")
        return [Ve, Vn, Se, Sn, rot, exx, exy, eyy, velfield_within_box, residual_velfield]


# ----------------- COORDINATES -------------------------

def lonlat_to_xy(lon, lat, lon0, lat0):
    """Project lon/lat (deg) to local cartesian (km) about (lon0, lat0), Okazaki-style."""
    x = (np.asarray(lon, dtype=float) - lon0) * np.cos(np.deg2rad(lat0)) * KM_PER_DEG
    y = (np.asarray(lat, dtype=float) - lat0) * KM_PER_DEG
    return x, y


# ----------------- CUBIC B-SPLINE (vectorized) -------------------------

def _spline3(t, g):
    """Cubic B-spline value as a function of t = x - x0 (array), grid spacing g."""
    t = np.asarray(t, dtype=float)
    f = np.zeros_like(t)
    a = np.abs(t)
    # |t| < g : 0.5|t|^3/g^3 - t^2/g^2 + 2/3   (written via the signed-region formulas)
    m1 = a < g
    f[m1] = 0.5 * a[m1] ** 3 - g * t[m1] ** 2 + (2.0 / 3.0) * g ** 3
    # g <= |t| < 2g : -(|t| - 2g)^3 / 6
    m2 = (a >= g) & (a < 2 * g)
    f[m2] = -((a[m2] - 2.0 * g) ** 3) / 6.0
    return f / g ** 3


def _spline3_dx(t, g):
    """First derivative of the cubic B-spline w.r.t. x, as a function of t = x - x0."""
    t = np.asarray(t, dtype=float)
    f = np.zeros_like(t)
    a = np.abs(t)
    s = np.sign(t)
    # |t| < g : derivative = (1.5 t^2 sign - 2 g t) ... use signed form
    m1 = a < g
    f[m1] = 1.5 * s[m1] * t[m1] ** 2 - 2.0 * g * t[m1]
    # g <= |t| < 2g : -0.5 sign (|t| - 2g)^2
    m2 = (a >= g) & (a < 2 * g)
    f[m2] = -0.5 * s[m2] * (a[m2] - 2.0 * g) ** 2
    return f / g ** 3


# ----------------- 1D B-SPLINE INTEGRAL TABLES (for roughness R) -------------------------
# X0 = integral of phi_k phi_l ; X1 = integral of phi'_k phi'_l ; X2 = integral of phi''_k phi''_l
# (with boundary B-splines truncated). Ported verbatim from the reference functions.py.

def _X0(k, l, M):
    if k == 0:
        if l == 0: return 1 / 252
        if l == 1: return 43 / 1680
        if l == 2: return 1 / 84
    elif k == 1:
        if l == 0: return 43 / 1680
        if l == 1: return 151 / 630
        if l == 2: return 531 / 2520
    elif k == 2:
        if l == 0: return 1 / 84
        if l == 1: return 531 / 2520
        if l == 2: return 599 / 1260
    elif k == M - 1:
        if l == M - 1: return 1 / 252
        if l == M - 2: return 43 / 1680
        if l == M - 3: return 1 / 84
    elif k == M - 2:
        if l == M - 1: return 43 / 1680
        if l == M - 2: return 151 / 630
        if l == M - 3: return 531 / 2520
    elif k == M - 3:
        if l == M - 1: return 1 / 84
        if l == M - 2: return 531 / 2520
        if l == M - 3: return 599 / 1260
    if k == l: return 151 / 315
    elif (k - l + 1) * (k - l - 1) == 0: return 397 / 1680
    elif (k - l + 2) * (k - l - 2) == 0: return 1 / 42
    elif (k - l + 3) * (k - l - 3) == 0: return 1 / 5040
    else: return 0.0


def _X1(k, l, M):
    if k == 0:
        if l == 0: return 1 / 20
        if l == 1: return 7 / 120
        if l == 2: return -1 / 10
    elif k == 1:
        if l == 0: return 7 / 120
        if l == 1: return 1 / 3
        if l == 2: return -11 / 60
    elif k == 2:
        if l == 0: return -1 / 10
        if l == 1: return -11 / 60
        if l == 2: return 37 / 60
    elif k == M - 1:
        if l == M - 1: return 1 / 20
        if l == M - 2: return 7 / 120
        if l == M - 3: return -1 / 10
    elif k == M - 2:
        if l == M - 1: return 7 / 120
        if l == M - 2: return 1 / 3
        if l == M - 3: return -11 / 60
    elif k == M - 3:
        if l == M - 1: return -1 / 10
        if l == M - 2: return -11 / 60
        if l == M - 3: return 37 / 60
    if k == l: return 2 / 3
    elif (k - l + 1) * (k - l - 1) == 0: return -1 / 8
    elif (k - l + 2) * (k - l - 2) == 0: return -1 / 5
    elif (k - l + 3) * (k - l - 3) == 0: return -1 / 120
    else: return 0.0


def _X2(k, l, M):
    if k == 0:
        if l == 0: return 1 / 3
        if l == 1: return -1 / 2
        if l == 2: return 0.0
    elif k == 1:
        if l == 0: return -1 / 2
        if l == 1: return 4 / 3
        if l == 2: return -1.0
    elif k == 2:
        if l == 0: return 0.0
        if l == 1: return -1.0
        if l == 2: return 7 / 3
    elif k == M - 1:
        if l == M - 1: return 1 / 3
        if l == M - 2: return -1 / 2
        if l == M - 3: return 0.0
    elif k == M - 2:
        if l == M - 1: return -1 / 2
        if l == M - 2: return 4 / 3
        if l == M - 3: return -1.0
    elif k == M - 3:
        if l == M - 1: return 0.0
        if l == M - 2: return -1.0
        if l == M - 3: return 7 / 3
    if k == l: return 8 / 3
    elif (k - l + 1) * (k - l - 1) == 0: return -3 / 2
    elif (k - l + 3) * (k - l - 3) == 0: return 1 / 6
    else: return 0.0


def _build_1d_integral_matrices(M):
    """Build the three banded 1D integral matrices (M x M) used to assemble R."""
    A0 = np.zeros((M, M))
    A1 = np.zeros((M, M))
    A2 = np.zeros((M, M))
    for k in range(M):
        # the integrals vanish for |k - l| > 3, so only fill the band
        for l in range(max(0, k - 3), min(M, k + 4)):
            A0[k, l] = _X0(k, l, M)
            A1[k, l] = _X1(k, l, M)
            A2[k, l] = _X2(k, l, M)
    return A0, A1, A2


# ----------------- CORE -------------------------

def compute_okazaki(myVelfield, xdata, ydata, grid_km, corr_leng=0.0, node_margin=2):
    """
    Estimate velocity & strain-rate fields via the Okazaki ABIC basis-function method.

    :returns: Ve, Vn, Se, Sn, rot, exx, exy, eyy  (each 2D array shaped (ny, nx))
    """
    lon = np.array([s.elon for s in myVelfield], dtype=float)
    lat = np.array([s.nlat for s in myVelfield], dtype=float)
    ve = np.array([s.e for s in myVelfield], dtype=float)
    vn = np.array([s.n for s in myVelfield], dtype=float)
    n_dat = lon.size

    # local cartesian frame about the data centroid
    lon0, lat0 = lon.mean(), lat.mean()
    sx, sy = lonlat_to_xy(lon, lat, lon0, lat0)

    # output (estimation) grid in cartesian coordinates
    gxx, gyy = np.meshgrid(xdata, ydata)
    ny, nx = gxx.shape
    ex, ey = lonlat_to_xy(gxx.ravel(), gyy.ravel(), lon0, lat0)

    # ---- basis node grid (cover stations AND estimation points, with margin) ----
    xmax = max(np.abs(sx).max(), np.abs(ex).max())
    ymax = max(np.abs(sy).max(), np.abs(ey).max())
    Xran = int(np.ceil(xmax / grid_km)) + node_margin
    Yran = int(np.ceil(ymax / grid_km)) + node_margin
    k_x = grid_km * np.arange(-Xran, Xran + 1)
    l_y = grid_km * np.arange(-Yran, Yran + 1)
    n_X, n_Y = k_x.size, l_y.size
    n_basis = n_X * n_Y
    print("Basis functions: %d x %d = %d   (grid_km=%g)" % (n_X, n_Y, n_basis, grid_km))
    # The 1D B-spline boundary-integral corrections assume the low- and high-edge
    # node regions do not overlap, which requires at least 7 nodes per dimension.
    if n_X < 7 or n_Y < 7:
        raise ValueError("okazaki: grid_km=%g is too coarse for this region "
                         "(basis grid is %dx%d; need >= 7 nodes per dimension). "
                         "Decrease grid_km." % (grid_km, n_X, n_Y))

    # ---- observation matrix H (vectorized; basis index j = l*n_X + k) ----
    H = _basis_matrix(sx, sy, k_x, l_y, grid_km, deriv=None)  # (n_dat, n_basis)

    # ---- roughness matrix R via Kronecker products of banded 1D integral matrices ----
    Ax0, Ax1, Ax2 = _build_1d_integral_matrices(n_X)
    Ay0, Ay1, Ay2 = _build_1d_integral_matrices(n_Y)
    sp = lambda M: sparse.csr_matrix(M)
    R = (sparse.kron(sp(Ay0), sp(Ax2)) + 2.0 * sparse.kron(sp(Ay1), sp(Ax1))
         + sparse.kron(sp(Ay2), sp(Ax0))).toarray()

    P = np.linalg.matrix_rank(R)               # rank (degrees of freedom of the prior)
    sign, lam = np.linalg.slogdet(R)
    if not np.isfinite(lam):
        lam = 0.0                              # only shifts ABIC by a constant

    # ---- data-error covariance E ----
    if corr_leng > 0:
        D = np.sqrt((sx[:, None] - sx[None, :]) ** 2 + (sy[:, None] - sy[None, :]) ** 2)
        E = np.exp(-D / corr_leng)
        Ei = np.linalg.inv(E)
        logdet_E = np.linalg.slogdet(E)[1]
        HtEi = H.T @ Ei
    else:
        Ei = None
        logdet_E = 0.0
        HtEi = H.T

    # normal-equation pieces that don't depend on the hyperparameter
    N0 = HtEi @ H                              # (n_basis, n_basis)
    rhs_x = HtEi @ ve
    rhs_y = HtEi @ vn

    # ---- ABIC hyperparameter search (5-level refinement, as in the reference) ----
    ax_opt, ay_opt, cov_opt, abic_min, loga_opt = _solve_abic(
        N0, R, rhs_x, rhs_y, H, ve, vn, Ei, n_dat, n_basis, P, lam, logdet_E)
    print("ABIC optimum: log10(alpha)=%.3f  ABIC=%.3f" % (loga_opt, abic_min))

    # ---- evaluate fields on the output grid (chunked over estimation points) ----
    Ve = np.full(ex.shape, np.nan)
    Vn = np.full(ex.shape, np.nan)
    Se = np.full(ex.shape, np.nan)
    Sn = np.full(ex.shape, np.nan)
    rot = np.full(ex.shape, np.nan)
    exx = np.full(ex.shape, np.nan)
    exy = np.full(ex.shape, np.nan)
    eyy = np.full(ex.shape, np.nan)

    chunk = 2000
    for c0 in range(0, ex.size, chunk):
        c1 = min(c0 + chunk, ex.size)
        cx, cy = ex[c0:c1], ey[c0:c1]
        Phi = _basis_matrix(cx, cy, k_x, l_y, grid_km, deriv=None)
        Phi_dx = _basis_matrix(cx, cy, k_x, l_y, grid_km, deriv='x')
        Phi_dy = _basis_matrix(cx, cy, k_x, l_y, grid_km, deriv='y')

        covered = np.abs(Phi).sum(axis=1) > 0   # points inside the basis support
        fx = Phi @ ax_opt
        fy = Phi @ ay_opt
        fxx = 1000.0 * (Phi_dx @ ax_opt)
        fxy = 1000.0 * (Phi_dy @ ax_opt)
        fyx = 1000.0 * (Phi_dx @ ay_opt)
        fyy = 1000.0 * (Phi_dy @ ay_opt)

        # posterior velocity variance: diag(Phi cov Phi^T)
        var_v = np.einsum('ij,ij->i', Phi, Phi @ cov_opt)
        var_v = np.clip(var_v, 0.0, None)

        idx = slice(c0, c1)
        Ve[idx] = np.where(covered, fx, np.nan)
        Vn[idx] = np.where(covered, fy, np.nan)
        Se[idx] = np.where(covered, np.sqrt(var_v), np.nan)
        Sn[idx] = Se[idx]
        exx[idx] = np.where(covered, fxx, np.nan)
        eyy[idx] = np.where(covered, fyy, np.nan)
        exy[idx] = np.where(covered, 0.5 * (fxy + fyx), np.nan)
        rot[idx] = np.where(covered, 0.5 * (fyx - fxy), np.nan)

    reshape = lambda a: a.reshape(ny, nx)
    return (reshape(Ve), reshape(Vn), reshape(Se), reshape(Sn),
            reshape(rot), reshape(exx), reshape(exy), reshape(eyy))


def _basis_matrix(px, py, k_x, l_y, grid_km, deriv=None):
    """
    Evaluate the 2D B-spline basis (or a first derivative) at points (px, py).
    Returns a dense (n_pts, n_basis) array with basis index j = l*n_X + k.

    :param deriv: None -> phi; 'x' -> d phi/dx; 'y' -> d phi/dy
    """
    tx = px[:, None] - k_x[None, :]            # (n_pts, n_X)
    ty = py[:, None] - l_y[None, :]            # (n_pts, n_Y)
    if deriv == 'x':
        Sx = _spline3_dx(tx, grid_km)
        Sy = _spline3(ty, grid_km)
    elif deriv == 'y':
        Sx = _spline3(tx, grid_km)
        Sy = _spline3_dx(ty, grid_km)
    else:
        Sx = _spline3(tx, grid_km)
        Sy = _spline3(ty, grid_km)
    # j = l*n_X + k  ->  reshape of (n_pts, n_Y, n_X)
    return (Sy[:, :, None] * Sx[:, None, :]).reshape(px.shape[0], -1)


def _abic_for_alpha(alpha, N0, R, rhs_x, rhs_y, H, ve, vn, Ei, n_dat, n_basis, P, lam, logdet_E):
    """ABIC value (and solution) for a given hyperparameter alpha. Returns (abic, ax, ay, sigma, factor)."""
    A = N0 + alpha * R
    factor = cho_factor(A, lower=True)
    logdet_A = 2.0 * np.sum(np.log(np.abs(np.diag(factor[0]))))
    ax = cho_solve(factor, rhs_x)
    ay = cho_solve(factor, rhs_y)
    rx = ve - H @ ax
    ry = vn - H @ ay
    if Ei is None:
        fit = rx @ rx + ry @ ry
    else:
        fit = rx @ (Ei @ rx) + ry @ (Ei @ ry)
    s = fit + alpha * (ax @ (R @ ax) + ay @ (R @ ay))
    NPM = n_dat + P - n_basis
    sigma = 0.5 * s / NPM
    LL1 = -NPM
    LL2 = -logdet_A + P * np.log(alpha) + lam
    LL3 = -NPM * np.log(np.pi * s / NPM)
    LL4 = -logdet_E
    LL = LL1 + LL2 + LL3 + LL4
    abic = -2.0 * LL + 4.0
    return abic, ax, ay, sigma, factor


def _solve_abic(N0, R, rhs_x, rhs_y, H, ve, vn, Ei, n_dat, n_basis, P, lam, logdet_E):
    """
    5-level refinement search over log10(alpha) (as in the reference solve_ABIC_opt),
    returning the optimal coefficients and posterior model covariance.
    """
    abic_min = np.inf
    loga_opt = 0.0
    best = None
    mode_opt = 0.0
    for n in range(5):
        for i in range(9):
            if n == 0:
                loga = 0.625 * (i + 1) - 3.125
            elif n == 1:
                loga = 0.125 * (i + 1) + mode_opt - 0.625
            elif n == 2:
                loga = 0.025 * (i + 1) + mode_opt - 0.125
            elif n == 3:
                loga = 0.005 * (i + 1) + mode_opt - 0.025
            else:
                loga = 0.001 * (i + 1) + mode_opt - 0.005
            alpha = 10.0 ** loga
            abic, ax, ay, sigma, factor = _abic_for_alpha(
                alpha, N0, R, rhs_x, rhs_y, H, ve, vn, Ei, n_dat, n_basis, P, lam, logdet_E)
            if abic < abic_min:
                abic_min = abic
                loga_opt = loga
                best = (ax, ay, sigma, factor)
        mode_opt = loga_opt

    ax_opt, ay_opt, sigma_opt, factor = best
    # posterior model covariance: sigma^2 * (N0 + alpha R)^-1
    cov_opt = sigma_opt * cho_solve(factor, np.eye(n_basis))
    return ax_opt, ay_opt, cov_opt, abic_min, loga_opt
