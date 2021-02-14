from abc import ABC, abstractmethod

class Strain_2d(ABC):
    """
    Implement a generic 2D strain rate method
    """
    def __init__(self):
        # Initialize general parameters
        self._Name = None
        self._grid_inc = None
        self._strain_range = None
        self._data_range = None

    def Method(self):
        return self._Name

    @abstractmethod
    def compute(self, myVelfield):
        # generic method to be implemented in each method
        pass
