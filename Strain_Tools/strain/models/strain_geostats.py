# Use a geostatistical interpolation scheme.
# Maurer, 2021. (in prep)

from abc import ABC, abstractmethod
import os

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform

from strain.models.strain_2d import Strain_2d
from .. import utilities
from strain.strain_tensor_toolbox import strain_on_regular_grid, calc_strain_uncertainty
from strain.utilities import getVels


class geostats(Strain_2d):
    """ 
    Geostatistical interpolation class for 2d strain rate, with general 
    strain_2d behavior 

    Parameters
    ----------
    params: dict - a dict-like containing at least the key word arguments 
                   'inc' and 'range_strain' for specifying the grid. Can
                   also have variogram parameters specified.
    model: callable of VariogramModel type - any of valid VariogramModel
                   subclasses. If specified, parameters other than inc and 
                   range_strain will not be used. 
    """
    def __init__(self, params=None, model=None):
        Strain_2d.__init__(
                self, 
                params.inc, 
                params.range_strain,
                params.range_data,
                params.xdata,
                params.ydata,
                params.outdir,
            )
        self._Name = 'geostatistical'
        if (model is not None) and callable(model):
            self._model_east = model
            self._model_north = model
        else:
            self._model_type = params.method_specific['model_type']
            self._model_east = self.getVariogram(
                    model_type=self._model_type,
                    sill=np.float64(params.method_specific['sill_east']),
                    rang=np.float64(params.method_specific['range_east']),
                    nugget=np.float64(params.method_specific['nugget_east']),
                    trend=bool(params.method_specific['trend']),
                )
            self._model_north = self.getVariogram(
                    model_type=self._model_type,
                    sill=np.float64(params.method_specific['sill_north']),
                    rang=np.float64(params.method_specific['range_north']),
                    nugget=np.float64(params.method_specific['nugget_north']),
                    trend=bool(params.method_specific['trend']),
                )
        self.setGrid(None)

    def __repr__(self):
        out_str = '''I'm a class for generating strain rates using kriging.
        {}
        '''.format(
            self._model.__repr__()
        )
        return out_str

    def getVariogram(
            self, 
            model_type,
            sill,
            rang,
            nugget=None,
            trend=False
            ):
        """
        Set the parameters of the spatial structure function
        Former docs:
        c0: list                    Single or list of lists containing the
                                    initial parameter guesses for the
                                    variogram model type
    
        Parameters
        ----------
        model_type: str or list     Covariance model type
        sill:
        rang:
        nugget:
        trend: boolean or list      Use a trend or not, can be a list of bool
        """
        try:
            _model = eval(model_type)
            _model = _model()
        except AttributeError:
            raise ValueError('Model "{}" has not been implemented'.format(model_type))
        _model.setParms(
            sill=sill,
            range=rang,
            use_nugget=[True if nugget is not None else False],
            nugget=nugget
        )
        return _model

    def setPoints(self, xy, data):
        """
        Set the data for kriging.

        Parameters
        ----------
        xy: 2D ndarray          XY locations of the input data
        data: 1 or 2D ndarray   Data to krige. If 2D, will do each dim 
                                separately
        """
        if data.ndim != 2:
            raise ValueError(
                'Input data should be an N x 2 matrix of [easting, northing]'
            )
        self._xy = xy
        self._easting = np.array(data)[:, 0]
        self._northing = np.array(data)[:, 1]

    def setGrid(self, XY=None, XY_shape=None):
        """
        Set the query point grid for kriging

        Parameters
        ----------
        XY: 2D ndarray  A 2-D array ordered [x y] of query points.
        XY_shape: tuple of two integers, (n, m), shape of XY
        """
        if XY is None:
            self._lons, self._lats = self._xdata, self._ydata
            self._grid_shape = (len(self._lats), len(self._lons))
            [X, Y] = np.meshgrid(self._lons, self._lats)
            self._XY = np.array([X.flatten(), Y.flatten()]).T
        else:
            self._XY = XY
            self._grid_shape = XY_shape

    def krige_east(self, model, ktype='ok'):
        """
        Interpolate velocities using kriging
        """
        return krige(self._xy, self._XY, self._easting, self._model_east, ktype=ktype)

    def krige_north(self, model, ktype='ok'):
        """
        Interpolate velocities using kriging
        """
        return krige(self._xy, self._XY, self._northing, self._model_north, ktype=ktype)

    def compute(self, myVelfield):
        """Compute the interpolated velocity field"""

        dlon, dlat, e, n, se, sn = getVels(myVelfield)
        xy = np.stack([dlon, dlat], axis=-1)
        data = np.stack([e, n], axis=1)
        self.setPoints(xy=xy, data=data)

        Dest_e, Dsig_e, _ = self.krige_east('ok')
        Dest_n, Dsig_n, _ = self.krige_north('ok')
        
        Ve = Dest_e.reshape(self._grid_shape)
        Vn = Dest_n.reshape(self._grid_shape)
        Se = Dsig_e.reshape(self._grid_shape)
        Sn = Dsig_n.reshape(self._grid_shape)

        # Compute strain rates
        dx, dy = self._grid_inc[0] * 111 * np.cos(np.deg2rad(self._strain_range[2])), self._grid_inc[1] * 111
        exx, eyy, exy, rot = strain_on_regular_grid(dx, dy, Ve, Vn)

        # Right now we aren't calculating uncertainties for any method except geostats.
        # We might want to consider exploring adding uncertainties to other methods if
        # they are amenable; e.g. local average gradient and visr *should* be able to 
        # provide them. Others may not (wavelets, gpsgridder esp). 
        var_dil, var_shear = calc_strain_uncertainty(np.square(Se), np.square(Sn), self._grid_inc[0], self._grid_inc[1], exx, eyy, exy)

        # Report observed and residual velocities within bounding box
        velfield_within_box = utilities.filter_by_bounding_box(myVelfield, self._strain_range);
        model_velfield = utilities.create_model_velfield(self._xdata, self._ydata, Ve, Vn, velfield_within_box);
        residual_velfield = utilities.subtract_two_velfields(velfield_within_box, model_velfield);

        # Return the strain rates etc. in the same units as other methods
        return Ve, Vn, rot*1000, exx*1000, exy*1000, eyy*1000, velfield_within_box, residual_velfield
        

def krige(xy, XY, data, model, ktype='ok'):
    """
    Interpolate velocities using kriging
    """
    # Create the data covariance matrix
    SIG = compute_covariance(model, xy)
    SIG = SIG + np.sqrt(np.finfo(float).eps)*np.eye(SIG.shape[0])
    if not is_pos_def(SIG):
        raise RuntimeError('SIG matrix is not positive definite, probably meaning it is ill-conditioned')

    # create the data/grid covariance and point-wise terms
    sig0 = compute_covariance(model, xy, XY); 
    sig2 = model.getSigma00()

    # Do different things, depending on if SK, OK, or UK is desired
    if ktype == 'sk':
        Dest, Dsig, lam = simple_kriging(SIG, sig0, data, sig2)
    elif ktype == 'ok':
        Dest, Dsig, lam, nu = ordinary_kriging(SIG, sig0, data, sig2)
    elif ktype == 'uk':
        Dest, Dsig, lam = universal_kriging(SIG, sig0, data, sig2, xy, XY)
    else:
        raise ValueError('Method "{}" is not implemented'.format(ktype))

    return Dest, Dsig, lam


def simple_kriging(SIG, sig0, data, sig2):
    """
    Perform simple (i.e. zero-mean) kriging

    Parameters
    ----------
    SIG: N x N ndarray 	- Covariance of all observed data locations
    sig0: N x M ndarray - Covariance of observed vs query locations
    data: N x 1 ndarray - Observed data values
    sig2: 1 x 1 float   - point-wise variance

    Returns
    -------
    Dest: M x 1 ndarray - Expected value of the field at the query locations
    Dsig: M x 1 ndarray - Sqrt of the kriging variance for each query location
    lam:  N x M ndarray - weights relating all of the data locations to all of the query locations
    """
    A = SIG
    B = sig0
    lam = np.linalg.lstsq(A, B, rcond=None)
    Dest = np.dot(lam.T, data)
    Dsig = []
    for k in range(lam.shape[1]):
        Dsig.append(np.sqrt(sig2 - np.dot(lam[:, k], sig0[:, k])))
    return Dest, np.array(Dsig), lam


def ordinary_kriging(SIG, sig0, data, sig2):
    """
    Perform ordinary kriging (stationary non-zero mean)

    Parameters
    ----------
    SIG: N x N ndarray 	- Covariance of all observed data locations
    sig0: N x M ndarray - Covariance of observed vs query locations
    data: N x 1 ndarray - Observed data values
    sig2: 1 x 1 float   - point-wise variance

    Returns
    -------
    Dest: M x 1 ndarray - Expected value of the field at the query locations
    Dsig: M x 1 ndarray - Sqrt of the kriging variance for each query location
    lam:  N x M ndarray - weights relating all of the data locations to all of the query locations
    """
    M, N = sig0.shape
    A = np.block([
        [SIG, np.ones((M, 1))],
        [np.ones((1, M)), 0],
    ])
    B = np.vstack([sig0, np.ones((1, N))])
    lam_nu, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    lam = lam_nu[:-1, :]; nu = -lam_nu[-1, :]
    Dest = np.dot(lam.T, data)
    Dsig = []
    for k in range(lam.shape[1]):
        Dsig.append(np.sqrt(sig2 - np.dot(lam[:, k], sig0[:, k]) + nu[k]))
    return Dest, np.array(Dsig), lam, nu


def universal_kriging(SIG, sig0, data, sig2, xy, XY):
    """
    Perform universal kriging (field can be a linear function of "space". In reality 
    "space" can be any auxilliary variables, such as xy-location, topographic height, etc.)

    Parameters
    ----------
    SIG: N x N ndarray 	- Covariance of all observed data locations
    sig0: N x M ndarray - Covariance of observed vs query locations
    data: N x 1 ndarray - Observed data values
    sig2: 1 x 1 float   - point-wise variance
    xy: N x D ndarray   - Location of the observed data in D-dimensional space
    XY: M x D ndarray   - Location of the query locations in D-dimensional space

    Returns
    -------
    Dest: M x 1 ndarray - Expected value of the field at the query locations
    Dsig: M x 1 ndarray - Sqrt of the kriging variance for each query location
    lam:  N x M ndarray - weights relating all of the data locations to all of the query locations
    """
    raise NotImplementedError


def compute_covariance(model, xy, XY=None):
    """Returns the covariance matrix for a given set of data"""
    if xy.size == 1:
        dist = 0
    elif XY is None:
        dist = squareform(pdist(xy))
    else:
        dist = cdist(xy, XY)

    C = model(dist)
    return C


class VariogramModel(ABC):
    """ Defines the base variogram model class"""
    def __init__(
            self,
            model_type
            ):
        self._model = model_type
        self._params = {}
        self._params['sill'] = None
        self._params['range'] = None
        self._params['nugget'] = None

    def __repr__(self):
      out_str = '''My name is {}
      My sill is {}
      My range is {}
      I'm {} a nugget: {}
      '''.format(
          self._model,
          self._params['sill'],
          self._params['range'],
          ["using" if self._params['nugget'] is not None else "not using"][0],
          self._params['nugget']
        )
      return out_str

    @abstractmethod
    def __call__(self, h):
        pass

    def getParms(self):
        return self._params['sill'], \
               self._params['range'], \
               self._params['nugget']

    def setParms(self, **kwargs):
      for key, value in kwargs.items():
        self._params[key] = value

    def getSigma00(self):
        return self._params['sill']

class Nugget(VariogramModel):
    """Implements a nugget model"""
    def __init__(self, nugget=None):
        VariogramModel.__init__(self, 'Nugget')

    def __call__(self, h):
        if self._params['nugget'] is None:
            raise RuntimeError('You must first specify a nugget')
        return self._params['nugget']*(h == 0)  # params is a scalar for nugget


class Gaussian(VariogramModel):
    """Implements a Gaussian model"""
    def __init__(self,
                 sill=None,
                 range=None,
                 nugget=None
                 ):
        VariogramModel.__init__(self, 'Gaussian')
        self.setParms(sill=sill, range=range, nugget=nugget)

    def __call__(self, h):
        if self._params['sill'] is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if not self._params['nugget'] is not None:
            return self._params['sill']*np.exp(-np.square(h / self._params['range']))
        else:
            return self._params['sill']*np.exp(-np.square(h / self._params['range'])) + self._params['nugget']*(h == 0)


class Exponential(VariogramModel):
    """Implements an Exponential model"""
    def __init__(self,
                 sill=None,
                 range=None,
                 nugget=None
                 ):
        VariogramModel.__init__(self, 'Exponential')
        self.setParms(sill=sill, range=range, nugget=nugget)

    def __call__(self, h):
        if self._params is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if not self._params['nugget'] is not None:
            return self._params['sill']*np.exp(-h/self._params['range'])
        else:
            return self._params['sill']*np.exp(-h/self._params['range']) + self._params['nugget']*(h == 0)  # add nugget


def is_pos_def(A):
    if np.array_equal(A, A.T):
        try:
            np.linalg.cholesky(A)
            return True
        except np.linalg.LinAlgError:
            return False
    else:
        return False