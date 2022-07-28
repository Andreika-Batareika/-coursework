from abc import ABCMeta, abstractmethod
import numpy as np
class MathModel():
    __metaclass__=ABCMeta
    def __init__(self):
       self.R=6# радиус оболочки ЛАЛВ в м.
       self._dY = np.zeros(0)
    @abstractmethod
    def RightPart(self, t, comeVector):
        pass
    @abstractmethod
    def writePatameters(vx,vy,vz):
        pass
