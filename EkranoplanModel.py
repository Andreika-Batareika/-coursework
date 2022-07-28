import numpy as np
from scipy import interpolate
from MathModel import MathModel
import math


class EkranoplanModel(MathModel):

#Конструктор
    def __init__(self):

        # Вектор-состояния модели
        self._dY = np.zeros([6])

        # Постоянные величины модели
        self.g = 9.81
        self.Jz = 30376.11
        self.B = 7.84
        self.S = 37.66
        self.m = 240
        self.P = 3400
        self.deltaRv = 1e-1
        self.deltaZCP = 1e-2
        self.deltaZEL = 1e-2
        self.W=2145
        self.R=9# радиус оболочки ЛАЛВ в м.



    # правые части

    def RightPart(self, t, comeVector):




         rv=1.29# плотность атмосферного воздуха в кг/м3.
         rg=0.17# плотность гелия в кг/м3.

         b=0.000125# константа, связанная с плотностью воздуха в 1/м
         a=6.5*10**-3# константа, связанная с температурой воздуха в К/м
         c=0.4#коэффициент лобового сопротивления
         mo=240#масса в кг
         V=(4/3)*3.14*self.R**3
         rs=rg+mo/V# суммарная плотность материала ЛАЛВ, массы гелия, и нагрузки
         p1=rv/rs# введенный параметр
         p2=3*c/(8*self.R)# введенный параметр
         T0=300
         g=9.8
         #сила тяжести
        # G = self.m * self.g
         #скоростной напор
         ##q = 0.5 * Ro * self.S * comeVector[4]**2

         # Правые части модели   +q/self.m+
         self._dY[0] = comeVector[3]
         # dTheta / dt, рад / сек
         self._dY[1] = comeVector[4]
         # dX / dt, м / сек
         self._dY[2] = comeVector[5]
         # dH / dt, м / сек
         self._dY[3] = np.random.normal(0,0.001)
         # dOmegaZ / dt, рад / сек^2
         #self._dY[4] = -G/self.m +fA/self.m #+ np.random.normal(0,0.001)
         self._dY[4]=-g+g*p1*math.exp(-b*comeVector[1]*T0/(T0-a*comeVector[1]))-np.sign(comeVector[4])* p1*p2*math.exp(-b*comeVector[1]*T0/(T0-a*comeVector[1]))*comeVector[4]**2        # dV / dt, м / сек^2

         # dPitch / dt, рад / сек
         self._dY[5] = np.random.normal(0,0.001)
         return self._dY
