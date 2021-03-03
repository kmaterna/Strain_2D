from abc import ABC, abstractmethod

class Strain_2d(ABC):
    """
    Implement a generic 2D strain rate method
    Strain_2d(grid_inc, strain_range, data_range)
    Parameters
    -------------
    grid_inc : list of [float, float]
          xinc, yinc in degrees for the rectilinear grid on which strain is calculated
    strain_range: list of [float, float, float, float]
          west, east, south, north edges of bounding box for grid on which strain is calculated
    data_range: list of [float, float, float, float]
          west, east, south, north edges of box that contains geodetic data (potentially larger than strain_range)

    Methods
    --------------
    compute(): required method that takes a velocity field and will eventually compute strain
    only provided as a template here
    """
    def __init__(self, grid_inc, strain_range, data_range):
        # Initialize general parameters
        self._Name = None
        self._grid_inc = grid_inc
        self._strain_range = strain_range
        self._data_range = data_range

    def Method(self):
        return self._Name

    @abstractmethod
    def compute(self, myVelfield):
        # generic method to be implemented in each method
        pass
