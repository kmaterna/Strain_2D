# Use a geostatistical interpolation scheme.
# Maurer, 2021. (in prep)

from abc import ABC, abstractmethod
import os

import numpy as np
from scipy.spatial.distance import pdist, cdist, squareform

from strain.models.strain_2d import Strain_2d
from strain.strain_tensor_toolbox import strain_on_regular_grid
from strain.utilities import makeGrid, getVels


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
                params.outdir,
            )
        self._Name = 'geostatistical'
        if (model is not None) and callable(model):
            self._model = model
        else:
            self._model_type = params.method_specific['model_type']
            self.setVariogram(
                    model_type = self._model_type,
                    sill = np.float64(params.method_specific['sill']),
                    rang = np.float64(params.method_specific['range']),
                    nugget = np.float64(params.method_specific['nugget']),
                    trend = bool(params.method_specific['trend']),
                )
        self.setGrid(None)

    def __repr__(self):
          out_str = '''I'm a class for generating strain rates using kriging.
          {}
          '''.format(
              self._model.__repr__()
          )
          return out_str

    def setVariogram(
            self, 
            model_type,
            sill,
            rang,
            nugget=None,
            trend=False,
        ):
        ''' 
        Set the parameters of the spatial structure function 
    
        Parameters
        ----------
        model_type: str or list     Covariance model type
        c0: list                    Single or list of lists containing the 
                                    initial parameter guesses for the 
                                    variogram model type
        trend: boolean or list      Use a trend or not, can be a list of bool
        '''
        try:
            _model = eval(model_type)
            self._model = _model()
        except AttributeError:
            raise ValueError('Model "{}" has not been implemented'.format(model_type))
        self._model.setParms(
            sill = sill,
            range = rang,
            use_nugget = [True if nugget is not None else False],
            nugget = nugget
        )

    def setPoints(self, xy, data):
        '''
        Set the data for kriging.

        Parameters
        ----------
        xy: 2D ndarray          XY locations of the input data
        data: 1 or 2D ndarray   Data to krige. If 2D, will do each dim 
                                separately
        '''
        if data.ndim > 2:
            raise ValueError('I can only handle 1- or 2-D data inputs')
        if data.ndim ==2:
            if data.shape[1] > 2:
                raise ValueError(
                        'Input data should be an N x 1 or N x 2 matrix'
                    )
        self._xy = xy
        self._values = np.array(data)
        self._data_dim = self._values.ndim

    def setGrid(self, XY=None, XY_shape=None):
        '''
        Set the query point grid for kriging

        Parameters
        ----------
        XY: 2D ndarray  A 2-D array ordered [x y] of query points. 
        '''
        if XY is None:
            self._lons, self._lats, self._XY = makeGrid(self._grid_inc[0], self._grid_inc[1], self._strain_range)
            self._grid_shape = (len(self._lats), len(self._lons))
        else:
            self._XY = XY
            self._grid_shape = XY_shape

    def krige(self, ktype='ok'):
        ''' 
        Interpolate velocities using kriging
        '''
        # Create the data covariance matrix
        SIG = compute_covariance(self._model, self._xy)
        SIG = SIG + np.sqrt(np.finfo(float).eps)*np.eye(SIG.shape[0])

        # create the data/grid covariance and point-wise terms
        sig0 = compute_covariance(self._model, self._xy, self._XY); 
        sig2 = self._model.getSigma00()

        # Do different things, depending on if SK, OK, or UK is desired
        if ktype == 'sk':
            Dest, Dsig, lam = simple_kriging(SIG, sig0, self._values, sig2)
        elif ktype == 'ok':
            Dest, Dsig, lam, nu = ordinary_kriging(SIG, sig0,self._values, sig2)
        elif ktype == 'uk':
            Dest, Dsig, lam = universal_kriging(SIG, sig0,self._values, sig2, self._xy, self._XY)
        else:
            raise ValueError('Method "{}" is not implemented'.format(self._flag))

        return Dest, Dsig, lam

    def compute(self, myVelfield):
        '''Compute the interpolated velocity field'''

        dlon, dlat, e, n, se, sn = getVels(myVelfield)
        xy = np.stack([dlon, dlat], axis=-1)
        self.setPoints(xy=xy, data = e)
        Dest_e, Dsig_e, _ = self.krige()
        self.setPoints(xy=xy, data = n)
        Dest_n, Dsig_n, _ = self.krige()
        
        # Compute strain rates
        exx, eyy, exy, rot = strain_on_regular_grid(
                self._grid_inc[0]*111, 
                self._grid_inc[1]*111, 
                Dest_e.reshape(self._grid_shape), 
                Dest_n.reshape(self._grid_shape)
            )

        # Return the strain rates etc.
        return self._lons, self._lats, rot, exx, exy, eyy
        

def simple_kriging(SIG, sig0, data, sig2):
    '''
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
    '''
    A = SIG
    B = sig0
    lam = np.linalg.lstsq(A,B, rcond=None)
    Dest = np.dot(lam.T, data)
    Dsig=[]
    for k in range(lam.shape[1]):
        Dsig.append(np.sqrt(sig2 - np.dot(lam[:,k], sig0[:,k])))
    return Dest, np.array(Dsig), lam


def ordinary_kriging(SIG, sig0, data, sig2):
    '''
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
    '''
    M,N = sig0.shape
    A = np.block([
        [SIG, np.ones((M,1))],
        [np.ones((1,M)), 0],
    ])
    B = np.vstack([sig0, np.ones((1,N))])
    lam_nu, _, _, _ = np.linalg.lstsq(A,B, rcond=None)
    lam = lam_nu[:-1,:]; nu = -lam_nu[-1,:]
    Dest = np.dot(lam.T, data)
    Dsig=[]
    for k in range(lam.shape[1]):
        Dsig.append(np.sqrt(sig2 - np.dot(lam[:,k], sig0[:,k]) + nu[k]))
    return Dest, np.array(Dsig), lam, nu


def universal_kriging(SIG, sig0, data, sig2, xy, XY):
    '''
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
    '''
    raise NotImplementedError


def compute_covariance(model, xy, XY=None):
    '''Returns the covariance matrix for a given set of data'''
    if xy.size==1:
        h = 0
    elif XY is None:
        h = squareform(pdist(xy))
    else:
        h = cdist(xy,XY)

    C = model(h)
    return C


class VariogramModel(ABC):
    ''' Defines the base variogram model class'''
    def __init__(
        self,
        model_type,
        use_nugget=False,
      ):
        self._model = model_type
        self._params = {}
        self._params['sill'] = None
        self._params['range'] = None
        self._params['nugget'] = None
        self._use_nugget = use_nugget

    def __repr__(self):
      out_str = '''My name is {}
      My sill is {}
      My range is {}
      I'm {} a nugget
      '''.format(
          self._model,
          self._params['sill'],
          self._params['range'],
          ["using" if self._use_nugget else "not using"][0]
        )
      return out_str

    @abstractmethod
    def __call__(self, h):
        pass

    def getParms(self):
        return self._params['sill'], self._params['range'], self._params['nugget']

    def nugget(self):
      '''This switches the use of nugget'''
      if self._use_nugget:
        self._use_nugget = False
      else:
        self._use_nugget = True

    def setParms(self, **kwargs):
      for key, value in kwargs.items():
        self._params[key] = value

    def getSigma00(self):
        return self._params['sill']

class Nugget(VariogramModel):
    '''Implements a nugget model'''
    def __init__(self, nugget = None):
        VariogramModel.__init__(self, 'Nugget')
    def __call__(self, h):
        if self._params['nugget'] is None:
            raise RuntimeError('You must first specify a nugget')
        return self._params['nugget']*(h != 0) # params is a scalar for nugget


class Gaussian(VariogramModel):
    '''Implements a Gaussian model'''
    def __init__(self,
            sill=None,
            range=None,
            nugget=None,
            use_nugget = False
        ):
        VariogramModel.__init__(self, 'Gaussian', use_nugget = use_nugget)
        self.setParms(sill = sill, range = range, nugget = nugget)

    def __call__(self, h):
        if self._params['sill'] is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if not self._use_nugget == 2:
            return self._params['sill']*(1 - np.exp(-(h**2)/(self._params['range']**2)))
        elif self._use_nugget == 3:
            return self._params['sill']*(1 - np.exp(-(h**2)/(self._params['range']**2))) + self._params['nugget']*(h != 0) # add a nugget


class Exponential(VariogramModel):
    '''Implements an Exponential model'''
    def __init__(self,
            sill=None,
            range=None,
            nugget=None,
            use_nugget = False
        ):
        VariogramModel.__init__(self, 'Exponential', use_nugget = use_nugget)
        self.setParms(sill = sill, range = range, nugget = nugget)

    def __call__(self, h):
        if self._params is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if len(self._params) == 2:
            return self._params['sill']*(1 - np.exp(-h/self._params['range']))
        elif len(self._params) == 3:
            return self._params['sill']*(1 - np.exp(-h/self._params['range'])) + self._params['nugget']*(h != 0) # add a nugget

