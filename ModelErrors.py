import numpy as np
from MathModel import MathModel


#Мат. модель экраноплана (продольный канал)
class ModelErrors(MathModel):

#Конструктор
    def __init__(self):

        # Вектор-состояния модели
        self._dY = np.zeros([11])

        # Постоянные величины модели из пач�import numpy as np
from MathModel import MathModel


#Мат. модель экраноплана (продольный канал)
class ModelErrors(MathModel):

#Конструктор
    def __init__(self):

        # Вектор-состояния модели
        self._dY = np.zeros([11])

        # Постоянные величины модели из пачпорта (непонятно откуда)
        self.W_B = 0.0
        self.W_P = 0.0
        self.alpha_CT = 0.005
        self.Sigma_CT = 400.
        self.Sigma_SIN = 100.
        self.w_0 =2*1000#20# 0.00005
        self.alpha_SIN = 0.0001
        self.alpha_TP = 1.0
        self.Sigma_TP = 35.0
        self.K_TP = 0.1
        self.betta = 0.1
        self.Sigma_P = 31.0
        self.Sigma_Relef=100.0
        self.M=0
        self.alpha_DIN = 0.1
        self.Sigma_DIN = 100.0
        self.P=0.2
        self.Sigma_T0 = 0



        self.W_B = 0.0
        self.W_P = 0.0
        self.alpha_CT = 0.005
        self.Sigma_CT = 400.
        self.Sigma_SIN = 100.
        self.w_0 =2*10000.0#20# 0.00005
        self.alpha_SIN = 0.0001
        self.alpha_TP = 1.0
        self.Sigma_TP = 35.0
        self.K_TP = 0.1
        self.betta = 0.1
        self.Sigma_P = 31.0
        self.Sigma_Relef=100.0
        self.M=0
        self.alpha_DIN = 1
        self.Sigma_DIN = 100.0
        self.P=0.2
        self.Sigma_T0 = 0


    # правые части

    def RightPart(self, t, comeVector):


         w=[np.random.normal(0,0.1), np.random.normal(0,0.1), np.random.normal(0,0.1), np.random.normal(0,0.1), np.random.normal(0,0.1), np.random.normal(0,0.1)]
         w=[np.random.normal(0,1), np.random.normal(0,1), np.random.normal(0,1), np.random.normal(0,1), np.random.normal(0,1), np.random.normal(0,1)]
         
         self.M= self.W_P/self.P

         #сила тяжести
        # G = self.m * self.g
         #скоростной напор
         ##q = 0.5 * Ro * self.S * comeVector[4]**2
         #силы лобового сопротивления и ..

         # Правые части модели   +q/self.m+
         
         
         self.W_B=np.abs(self.W_B)
         if(self.W_B>18):
             self.W_B=18
         #if(self.W_B<-18):
         #    self.W_B=-18
         #self.W_B=np.abs(self.W_B)

         self._dY[0] = comeVector[1]
         # dTheta / dt, рад / сек
         self._dY[1] = -(0.67*self.W_B+self.alpha_CT)*comeVector[1]-0.67*self.W_B*self.alpha_CT*comeVector[0]+np.sqrt(2*self.alpha_CT*self.Sigma_CT*self.Sigma_CT)*w[0]
         # dX / dt, м / сек
         self._dY[2] = comeVector[3]+w[1]
         # dH / dt, м / сек
         self._dY[3] = -(self.alpha_SIN+self.w_0)*comeVector[2]-2*self.alpha_SIN*comeVector[3]+(np.sqrt(2*self.alpha_SIN+self.w_0*self.w_0)-2*self.alpha_SIN)*w[1]
         # dOmegaZ / dt, рад / сек^2
         self._dY[4] = comeVector[5] #+ np.random.normal(0,0.001)

         self._dY[5]=- self.alpha_DIN*comeVector[5]+self.K_TP*self.alpha_TP*w[2]      # dV / dt, м / сек^2
         self._dY[5]=- self.alpha_TP*comeVector[5]+self.K_TP*self.alpha_TP*w[2]      # dV / dt, м / сек^2
         
         # dPitch / dt, рад / сек
         self._dY[6] = 0

         self._dY[7] = self.alpha_DIN*w[3]
         self._dY[7] = -self.alpha_DIN*w[3]
         self._dY[8] = -self.betta*comeVector[8]+np.sqrt(2*self.betta*self.Sigma_P*self.Sigma_P)*w[4]
         self._dY[9] = comeVector[10]
         self._dY[10] = -2*self.M*comeVector[10]-self.M*self.M*comeVector[9]+np.sqrt(4*self.M*self.M*self.M*self.Sigma_Relef*self.Sigma_Relef)*w[5]
         #print(self._dY)
         return self._dY



    def writePatameters(vx,vy,vz):
        self.M = np.sqrt(vx*vx+vz*vz)/self.P
        self.W_B = vy
