from abc import ABC, abstractmethod

class Strain_2d(ABC):
    """
    Implement a generic 2D strain rate method
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
