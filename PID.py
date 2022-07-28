from abc import ABCMeta, abstractmethod
import numpy as np
class PID():
    def __init__(self,hCome):

       self.dt=0.001
       self.h = hCome
       self.error=0
       self.K1=0
       self.K2= 0
       self.K3 =0


    def PIDFunc(self, t, hCome):
        self.error=hCome- self.h
        self.K3 =(self.error-self.K1)/self.dt
        self.K1=self.error
        self.K2= self.K2+self.error*self.dt
        P=-0.001
        I=-0.0000001
        D=-0.01
        return self.K1*P, self.K2*I, self.K3*D
