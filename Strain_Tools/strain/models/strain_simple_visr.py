"""
Simplified version of visr method.
"""

import numpy as np
import scipy as sp
import scipy.sparse as sparse
from scipy.stats import circmean
from scipy.spatial import Voronoi, ConvexHull
from scipy.spatial.distance import cdist
from scipy.optimize import brentq
from matplotlib.path import Path as MplPath
from typing import Literal
from pyproj import Proj

from strain import utilities
from strain.utilities import getVels
from strain.models.strain_2d import Strain_2d


class simple_visr(Strain_2d):
    """
    Implement a generic 2D strain rate method
    Strain_2d(grid_inc, strain_range, data_range)
    Parameters
    -------------
    grid_inc : list of [float, float]
          xinc, yinc in degrees for the rectilinear grid on which strain is calculated
    strain_range: list of [float, float, float, float]
          west, east, south, north edges of bounding box for grid on which strain is calculated
    xdata: 1d array of lon numbers of output grid points
    ydata: 1d array of lat numbers of output grid points
    data_range: list of [float, float, float, float]
          west, east, south, north edges of box that contains geodetic data
          (potentially larger than strain_range)
          NOT USED by this method

    Methods
    --------------
    compute(): required method that takes a velocity field and will eventually compute strain
    only provided as a template here
    """
    def __init__(self, params):
        # standard initialization
        super().__init__(params.inc, params.range_strain, params.range_data,
                         params.xdata, params.ydata, params.outdir)
        self._Name = "simple_visr"
        # read required method-specific parameters
        # can fail if key is missing or cannot be cast to correct data type
        self.weighting_threshold = float(params.method_specific["weighting_threshold"])
        self.distance_method = str(params.method_specific["distance_method"])
        self.coverage_method = str(params.method_specific["coverage_method"])
        # read optional method-specific parameters
        self.utmzone = params.method_specific.get("utmzone", None)
        self.estimate_within = params.method_specific.get("estimate_within", None)
        # validate parameters
        if self.weighting_threshold <= 0:
            raise ValueError("'weighting_threshold' needs to be positive.")
        if self.distance_method not in ["gaussian", "quadratic"]:
            raise ValueError("'distance_method' needs to be one of ['gaussian', 'quadratic']")
        if self.coverage_method not in ["azimuth", "voronoi"]:
            raise ValueError("'coverage_method' needs to be one of ['azimuth', 'voronoi']")
        if self.utmzone is not None:
            self.utmzone = int(self.utmzone)  # can fail with wrong data type
        if self.estimate_within is not None:
            self.estimate_within = float(self.estimate_within)  # can fail with wrong data type
            if self.estimate_within <= 0:
                raise ValueError("'estimate_within' needs to be None or positive.")

    @staticmethod
    def best_utmzone(longitudes: np.ndarray) -> int:
        """
        Given a list of longitudes, find the UTM zone that is appropriate.

        Parameters
        ----------
        longitudes
            Array of longitudes [°].

        Returns
        -------
            UTM zone at the average input longitude.
        """
        lon_mean = np.rad2deg(circmean(np.deg2rad(longitudes),
                                       low=-np.pi, high=np.pi))
        utmzone = int(np.ceil(((lon_mean + 180) / 6)))
        return utmzone

    @staticmethod
    def _prepare_vel_strain_rot(locations: np.ndarray,
                                velocities: np.ndarray,
                                covariances: np.ndarray | None,
                                utmzone: int,
                                reference: int | list | np.ndarray,
                                ) -> tuple[np.ndarray,
                                           np.ndarray | None,
                                           np.ndarray | None,
                                           np.ndarray]:
        # input checks
        assert (isinstance(locations, np.ndarray) and locations.ndim == 2 and
                locations.shape[1] == 2), \
            "'locations' needs to be a NumPy array with two columns."
        assert (isinstance(velocities, np.ndarray) and velocities.ndim == 2 and
                velocities.shape[1] == 2), \
            "'velocities' needs to be a NumPy array with two columns."
        assert locations.shape[0] == velocities.shape[0], \
            f"Mismatch between locations shape {locations.shape} and velocities " \
            f"shape {velocities.shape}."
        num_stations = locations.shape[0]
        if covariances is not None:
            assert (isinstance(covariances, np.ndarray) and covariances.ndim == 2 and
                    covariances.shape[0] == num_stations and
                    covariances.shape[1] == 2), \
                "Invalid covariance input type or shape."
        # parse reference
        if isinstance(reference, np.ndarray):
            assert reference.ndim == 2 and reference.shape[1] == 2, \
                f"'reference' needs to be a two-column array, got shape {reference.shape}."
        else:
            raise ValueError(f"Invalid input for 'reference': {reference}.")
        # make sure we're not inverting if we don't have enough data points
        if num_stations < 6:
            raise ValueError(f"{num_stations} stations are less stations than "
                             "necessary (6) for a stable computation.")
        # convert to UTM eastings & northings
        transformer = Proj(proj="utm", zone=utmzone)
        EN = np.stack(transformer(locations[:, 0], locations[:, 1]), axis=1)
        ENO = np.stack(transformer(reference[:, 0], reference[:, 1]), axis=1)
        # stack observations
        d = np.concatenate([velocities[:, 0], velocities[:, 1]], axis=0)
        # build covariances
        if covariances is not None:
            W = sparse.diags(1 / np.concatenate([covariances[:, 0], covariances[:, 1]]))
        else:
            W = None
        # done
        return EN, ENO, W, d

    @staticmethod
    def _get_vel_strain_rot_mapping(dE: np.ndarray | None = None,
                                    dN: np.ndarray | None = None,
                                    num_stations: int | None = None,
                                    out: np.ndarray | None = None
                                    ) -> np.ndarray:
        # input checks
        if num_stations is None:
            assert (dE is not None) and (dN is not None)
            num_stations = dE.shape[0]
        # build output array if necessary
        if out is None:
            G = np.zeros((2 * num_stations, 6))
            G[:num_stations, 4] = 1
            G[num_stations:, 5] = 1
        else:
            assert out.shape == (2 * num_stations, 6)
            G = out
        # add data if present
        if dE is not None:
            G[:num_stations, 0] = dE
            G[num_stations:, 2] = dE
        if dN is not None:
            G[:num_stations, 1] = dN
            G[num_stations:, 3] = dN
        # done
        return G

    @staticmethod
    def get_field_vel_strain_rot(locations: np.ndarray,
                                 velocities: np.ndarray,
                                 field: np.ndarray,
                                 weighting_threshold: float,
                                 covariances: np.ndarray | None,
                                 utmzone: int,
                                 distance_method: Literal["gaussian", "quadratic"],
                                 coverage_method: Literal["azimuth", "voronoi"],
                                 estimate_within: float | None
                                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""
        For a set of horizontal velocities on a 2D cartesian grid, estimate the
        interpolated velocity, strain, and rotation at target locations assuming
        a locally homogenous velocity field. See [tape09]_ and [shen15]_ for an introduction.

        This function uses a local approximation to the spherical Earth by
        converting all station and target locations into a suitable UTM zone, and only
        considering the horizontal velocities.

        Parameters
        ----------
        locations
            Array of shape :math:`(\text{num_stations}, 2)` containing the
            longitude and latitude [°] of the observations (stations).
        velocities
            Array of shape :math:`(\text{num_stations}, 2)` containing the
            East and North velocities [m/time] of the observations
        field
            Array of shape :math:`(\text{num_field}, 2)` containing the
            longitude and latitude [°] of the target field coordinates.
        weighting_threshold
            Weighting parameter that determines the smoothness of the output
            (:math:`W_t` in [shen15]_).
        covariances
            Array of shape :math:`(\text{num_stations}, 2)` containing the
            variances in the East and North velocities [m^2/time^2]. Alternatively,
            array of shape :math:`(\text{num_stations}, 3)` additionally
            containing the East-North covariance [m2/time^2].
        utmzone
            The UTM zone to use for the horizontal approximation.
        distance_method
            The method used to calculate the distance weighting starting at a
            target point to all stations (:math:`L_i` in [shen15]_).
        coverage_method
            The method used to calculate the distance weighting starting at a
            target point to all stations (:math:`Z_i` in [shen15]_).
        estimate_within
            If set, only estimate the values at target points that are this distance [m]
            away from the convex hull of all stations.

        Returns
        -------
        v
            Velocity of the field :math:`\mathbf{v}`, shape :math:`(\text{num_field}, 2)`.
        epsilon
            :math:`2 \times 2` strain tensor :math:`\mathbf{\varepsilon}` field,
            shape :math:`(\text{num_field}, 2, 2)`.
        omega
            :math:`2 \times 2` rotation tensor :math:`\mathbf{\omega}` field,
            shape :math:`(\text{num_field}, 2, 2)`.

        References
        ----------

        .. [tape09] Tape, C., Musé, P., Simons, M., Dong, D., & Webb, F. (2009),
           *Multiscale estimation of GPS velocity fields*,
           Geophysical Journal International, 179(2), 945–971,
           doi:`10.1111/j.1365-246X.2009.04337.x <https://doi.org/10.1111/j.1365-246X.2009.04337.x>`_.

        .. [shen15] Shen, Z.-K., Wang, M., Zeng, Y., Wang, F. (2015),
           *Optimal Interpolation of Spatially Discretized Geodetic Data*,
            Bulletin of the Seismological Society of America, 105(4), 2117–2127,
            doi:`10.1785/0120140247 <https://doi.org/10.1785/0120140247>`_.
        """
        # prepare
        assert isinstance(field, np.ndarray) and field.ndim == 2 and field.shape[1] == 2, \
            "'field' has to be a two-column array."
        EN, ENf, Wdata, d = \
            simple_visr._prepare_vel_strain_rot(locations, velocities, covariances, utmzone, field)
        num_stations = EN.shape[0]
        num_field = ENf.shape[0]
        all_vectors = EN[:, :2][None, :, :] - ENf[:, :2][:, None, :]
        if Wdata is None:
            Wdata = np.eye(2 * num_stations)  # C^-1 in Shen paper
        else:  # need dense representation
            Wdata = Wdata.toarray()

        # buld convex hull of network
        chull = ConvexHull(EN)
        chull_path = MplPath(chull.points[chull.vertices])
        # find field points that are inside the station network
        if estimate_within is None:
            estimate_within = np.sort(cdist(ENf[:, :2], ENf[:, :2]), axis=1)[:, 1].min()
        ix_field_inside = np.flatnonzero(chull_path.contains_points(ENf, radius=estimate_within))
        num_field_inside = ix_field_inside.size

        # define distance weighting function for various scale distances (L_i in Shen paper)
        all_dists = cdist(ENf[ix_field_inside, :2], EN[:, :2])
        if distance_method == "gaussian":
            def dist_weight_fun(scale, dist=all_dists): return np.exp(-(dist / scale)**2)
        elif distance_method == "quadratic":
            def dist_weight_fun(scale, dist=all_dists): return 1 / (1 + (dist / scale)**2)

        # define coverage weighting function (Z_i in Shen paper)
        if coverage_method == "azimuth":
            # get direction from all field points to all stations
            all_azimuths = np.arctan2(all_vectors[ix_field_inside, :, 1],
                                      all_vectors[ix_field_inside, :, 0])
            # sort by azimuth
            all_sortings = np.argsort(all_azimuths, axis=1)
            all_azimuths_sorted = np.take_along_axis(all_azimuths, all_sortings, axis=1)
            # calculate the difference in azimuths between each station
            all_azimuths_sorted = np.concatenate(
                [all_azimuths_sorted,
                 (all_azimuths_sorted[:, 0] + 2 * np.pi)[:, None]],
                axis=1)
            all_thetas = np.diff(all_azimuths_sorted, axis=1)
            # sum both sides of azimuthal difference
            all_thetas = np.concatenate(
                [all_thetas[:, -1][:, None],
                 all_thetas],
                axis=1)
            all_thetas = all_thetas[:, :-1] + all_thetas[:, 1:]
            # calculate final weight
            coverage_weight = num_stations * all_thetas / (4 * np.pi)
        elif coverage_method == "voronoi":
            # build Voronoi network of observations
            vor = Voronoi(EN)
            # calculate areas, default is based on closest neighbors
            vor_areas = np.pi * \
                np.mean(np.sort(cdist(EN[:, :2], EN[:, :2]), axis=1)[:, 1:7], axis=1)**2
            # update with actual area if not infinite or larger than twice the default
            for i_station, i_region in enumerate(vor.point_region):
                if i_station not in chull.vertices:
                    indices = vor.regions[i_region]
                    if -1 not in indices:
                        area = ConvexHull(vor.vertices[indices]).volume
                        if area < 2 * vor_areas[i_station]:
                            vor_areas[i_station] = area
            # calculate final weights
            coverage_weight = np.broadcast_to(
                num_stations * vor_areas / np.sum(vor_areas),
                (num_field_inside, num_stations))

        # calculate optimal scale distance
        def total_weight_fun_elem(scale, i):
            return np.sum(dist_weight_fun(scale, all_dists[i, :]) * coverage_weight[i, :]
                        ) - weighting_threshold
        optim_x1 = np.max(all_dists, axis=1)
        optimal_dist_scales = np.array(
            [brentq(total_weight_fun_elem, 1, optim_x1[i], args=(i, ))
            for i in range(num_field_inside)])
        if np.any(np.isclose(optimal_dist_scales, 1)):
            print("Warning: optimal distance scale probably not found!")
        # calculate final optimal weights (L_i in Shen paper)
        distance_weight = dist_weight_fun(optimal_dist_scales[:, None])

        # calculate joint covariance matrix (G_i in Shen paper)
        joint_weight = coverage_weight * distance_weight

        # build dummy G (from Tape paper, not Shen) (matrix A in Shen)
        G = simple_visr._get_vel_strain_rot_mapping(num_stations=num_stations)

        # create empty field output
        v = np.full((num_field, 2), np.nan)
        epsilon = np.full((num_field, 2, 2), np.nan)
        omega = np.full((num_field, 2, 2), np.nan)

        # start loop
        for i_field_sub, i_field in enumerate(ix_field_inside):
            # get local G
            G = simple_visr._get_vel_strain_rot_mapping(
                dE=all_vectors[i_field, :, 0],
                dN=all_vectors[i_field, :, 1],
                out=G)
            # get local W (C^-1 in Shen paper)
            # this is technically different from the paper because of the off-diagonals
            W = Wdata @ np.diag(np.tile(joint_weight[i_field_sub, :], 2))
            # get weighted least squares input
            GtWG = G.T @ W
            GtWd = GtWG @ d
            GtWG = GtWG @ G
            # solve
            m = sp.linalg.lstsq(GtWG, GtWd)[0]
            # extract wanted quantities
            L = m[:4].reshape(2, 2)
            v[i_field, :] = m[4:]  # velocity of point
            epsilon[i_field, :, :] = (L + L.T) / 2  # strain rate
            omega[i_field, :, :] = (L - L.T) / 2  # rotation rate

        # done
        return v, epsilon, omega

    def compute(self, myVelfield):
        print("------------------------------\n"
              "Running the Simple VISR strain computation method")
        # create output meshgrid
        xgrd, ygrd = np.meshgrid(self._xdata, self._ydata)
        field = np.stack([xgrd.ravel(), ygrd.ravel()], axis=1)
        # extract NumPy arrays
        lon, lat, e, n, se, sn = getVels(myVelfield)
        locations = np.stack([lon, lat], axis=1)
        velocities = np.stack([e, n], axis=1) / 1e3  # convert to m
        uncertainties = (np.stack([se, sn], axis=1) / 1e3) ** 2  # convert from s.d. to variance
        print(f"Analyzing {lon.size} observers to compute "
              f"{self._xdata.size * self._ydata.size} field points")
        # determine UTM zone if needed
        if self.utmzone is None:
            self.utmzone = simple_visr.best_utmzone(locations[:, 0])
            print("Calculated best-matching UTM zone")
        print(f"Using UTM zone {self.utmzone}")
        # compute
        print("Running computation...", flush=True, end="")
        v, epsilon, omega = simple_visr.get_field_vel_strain_rot(
            locations, velocities, field, self.weighting_threshold, uncertainties,
            self.utmzone, self.distance_method, self.coverage_method, self.estimate_within)
        print(" done")
        print("Formatting outputs...", flush=True, end="")
        # reformat velocity field
        ve = v[:, 0].reshape(self._ydata.size, self._xdata.size) * 1e3  # [mm/a]
        vn = v[:, 1].reshape(self._ydata.size, self._xdata.size) * 1e3  # [mm/a]
        # reformat epsilon and omega outouts
        rot = omega[:, 1, 0].reshape(self._ydata.size, self._xdata.size) * 1e9  # [rad/Ga]
        exx = epsilon[:, 0, 0].reshape(self._ydata.size, self._xdata.size) * 1e9  # [nanostrain/a]
        exy = epsilon[:, 0, 1].reshape(self._ydata.size, self._xdata.size) * 1e9  # [nanostrain/a]
        eyy = epsilon[:, 1, 1].reshape(self._ydata.size, self._xdata.size) * 1e9  # [nanostrain/a]
        # create model velfield by using the nearest-neighbor approach used elsewhere
        velfield_within_box = utilities.filter_by_bounding_box(myVelfield, self._strain_range)
        model_velfield = utilities.create_model_velfield(
            self._xdata, self._ydata, ve, vn, velfield_within_box)
        residual_velfield = utilities.subtract_two_velfields(velfield_within_box, model_velfield)
        print(" done\n")
        # done
        return [ve, vn, rot, exx, exy, eyy, velfield_within_box, residual_velfield]
