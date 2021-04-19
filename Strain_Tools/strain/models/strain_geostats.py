# Use a geostatistical interpolation scheme.
# Maurer, 2021. (in prep)

from abc import ABC, abstractmethod
import os

import numpy as np

from Strain_Tools.strain.models.strain_2d import Strain_2d


class geostats(Strain_2d):
    """ 
    Geostatistical interpolation class for 2d strain rate, with general 
    strain_2d behavior 
    """
    def __init__(self, params):
        strain_2d.Strain_2d.__init__(
                self, params.inc, params.range_strain, params.range_data
            )
        self._Name = 'geostatistical'
        self.setVariogram(
                params.model_type, 
                params.sill, 
                params.range, 
                params.nugget, 
                params.trend
            )
        self.setGrid(params.XY)

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
            self._model = eval(model_type)
        except AttributeError:
            raise ValueError('Model "{}" has not been implemented'.format(model_type)) 
        self._sill = float(sill)
        self._range = float(rang)
        self._nugget = float(nugget)
        self._trend = bool(trend)

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

    def setGrid(self, XY=None):
        '''
        Set the query point grid for kriging

        Parameters
        ----------
        XY: 2D ndarray  A 2-D array ordered [x y] of query points. 
        '''
        if XY is None:
            self._XY = makeGrid(self._grid_inc, self._grid_inc, self._strain_range)
        else:
            self._XY = XY

    def krige(self):
        ''' 
        Interpolate velocities using kriging
        '''
        # Create the data covariance matrix
        SIG = compute_covariance(model, param, xy)
        SIG = SIG + np.sqrt(np.finfo(float).eps)*np.eye(SIG.shape[0])

        # create the data/grid covariance and point-wise terms
        sig0 = compute_covariance(model, param, xy, XY); 
        sig2 = compute_covariance(model, param, 0, 0, 1); 

        # Do different things, depending on if SK, OK, or UK is desired
        if self._flag == 'sk':
            Dest, Dsig, lam = simple_kriging(SIG, sig0,data, sig2)
        elif self._flag == 'ok':
            Dest, Dsig, lam = ordinary_kriging(SIG, sig0,data, sig2)
        elif self._flag == 'uk':
            Dest, Dsig, lam = universal_kriging(SIG, sig0,data, sig2, xy, XY)
        else:
            raise ValueError('Method "{}" is not implemented'.format(self._flag))

        return Dest, Dsig, lam

    def compute(self, myVelfield):
        '''Compute the interpolated velocity field'''
        lon, lat, e, n, se, sn = getVel(myVelfield)
        xy = np.stack([lon, lat], axis=-1)
        self.setPoints(xy=xy, data = e)
        Dest_e, Dsig_e, _ = self.krige()
        self.setPoints(data = n)
        Dest_n, Dsig_n, _ = self.krige()
        
        # Compute strain rates
        dx = self._grid_inc
        dy = self._grid_inc # At the moment I only handle uniform grids
        exx, eyy, exy, rot = strain(dx, dy, Dest_e, Dest_n)

        # Return the strain rates etc.
        return lon, lat, rot, exx, exy, eyy
        

def strain(dx,dy,V1,V2):
    '''Compute strain rate on a regular grid'''
    dV1dx1, dV1dx2 = np.gradient(V1, dx, dy)
    dV2dx1, dV2dx2 = np.gradient(V2, dx, dy)

    # shear strain rate is the symmetric part of the gradient
    e12 = 0.5*(dV1dx2 + dV2dx1)

    # Rotation is the anti-symmetric component of the displacmenet gradient tensor
    rot = 0.5*(dV1dx2 - dV2dx1)

    # all the strain rates in each grid location
    return dV1dx1, dV2dx2, e12, rot


def getVels(velField):
    '''Read velocities from a NamedTuple'''
    for item in myVelfield:
        lon.append(item.elon)
        lat.append(item.elat)
        e.append(item.e)
        n.append(item.n)
        se.append(item.se)
        sn.append(item.sn)
    return np.array(lon), np.array(lat), np.array(e), np.array(n), np.array(se), np.array(sn)


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
    lam = np.linalg.lstsq(A,B)
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
        [SIG, np.ones(M,1)],
        [np.ones(1,M), 0],
    ])
    B = np.vstack([sig0, np.ones(1,N)])
    lam_nu = np.linalg.lstsq(A,B)
    lam = lam_nu[:-1,:]; nu = -lam_nu[-1,:]
    Dest = np.dot(lam.T, data)
    Dsig=[]
    for k in range(lam.shape[1]):
        Dsig.append(np.sqrt(sig2 - np.dot(lam[:,k], sig0[:,k]) + nu))
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


def makeGrid(gridx, gridy, bounds):
    '''Create a regular grid for kriging'''
    x = np.arange(bounds[0], bounds[1], gridx) 
    y = np.arange(bounds[2], bounds[3], gridy) 
    [X, Y] = np.meshgrid(x, y)
    return x, y, np.array([X.flatten(), Y.flatten()]).T


def compute_covariance(model, xy, XY=None):
    '''Returns the covariance matrix for a given set of data'''
    if xy.size==1:
        h = 0
    elif XY is None:
        h = pdist(xy)
    else:
        h = cdist(xy,XY)

    C = model(h)
    return C


class VariogramModel(ABC):
    ''' Defines the base variogram model class'''
    def __init__(self, model_type, useNugget = False):
        self._model = model_type
        self._params = None
        self._nugget = useNugget
    @abstractmethod
    def __call__(self, h):
        pass
    def getParms(self):
        return self._params


class Nugget(VariogramModel):
    '''Implements a nugget model'''
    def __init__(self):
        VariogramModel.__init__(self, 'Nugget')
    def __call__(self, h):
        if self._params is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        return self._params*(h != 0) # params is a scalar for nugget


class Gaussian(VariogramModel):
    '''Implements a Gaussian model'''
    def __init__(self, useNugget = False):
        VariogramModel.__init__(self, 'Gaussian', useNugget)
    def __call__(self, h):
        if self._params is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if len(self._params) == 2:
            return self._params[0]*(1 - np.exp(-(h**2)/(self._params[1]**2)))
        elif len(self._params) == 3:
            return self._params[0]*(1 - np.exp(-(h**2)/(self._params[1]**2))) + self._params[2]*(h != 0) # add a nugget


class Exponential(VariogramModel):
    '''Implements an Exponential model'''
    def __init__(self, useNugget = False):
        VariogramModel.__init__(self, 'Gaussian', useNugget)
    def __call__(self, h):
        if self._params is None:
            raise RuntimeError('You must first specify the parameters of the {} model'.format(self._model))
        if len(self._params) == 2:
            return self._params[0]*(1 - np.exp(-h/self._params[1]))
        elif len(self._params) == 3:
            return self._params[0]*(1 - np.exp(-h/self._params[1])) + self._params[2]*(h != 0) # add a nugget

