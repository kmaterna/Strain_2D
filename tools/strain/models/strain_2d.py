from abc import ABC, abstractmethod

class Strain_2d(ABC):
    '''
    Implement a generic 2D strain rate method
    '''
    def __init__(self, method='delaunay'):
        # Initialize general parameters
        self._Name = method

    def __str__(self):
        # print method
        raise NotImplementedError

    def Method(self):
        return self._Name

    @abstractmethod
    def compute(self, myVelfield, MyParams):
        # generic method to be implemented in each method
        pass

    def write(self):
        # Write results to files
        raise NotImplementedError
