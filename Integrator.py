from abc import ABCMeta, abstractmethod
from MathModel import MathModel
from ModelErrors import ModelErrors
import numpy as np
import PID
import numpy.random


def X(F,X_pred,K,H,N,r):

    temp_KH=K.dot(H) # K(t)*H(t)

    #inv_N=np.linalg.inv(N) # N обратная
    inv_N=1/N
    temp_KHN=temp_KH.dot(inv_N)
    delta_Y=r-H.dot(X_pred) # r(t)-H(t)*X(t)
    s=np.dot(temp_KHN,delta_Y[0])

    temp_X=F.dot(X_pred)+ s
    return temp_X
def K(F,K,V,N,Nw,H):
    temp_VN=V.dot(Nw) # V*N
    temp_KH_T=K.dot(H.transpose())
    temp_KHN=temp_KH_T.dot(1/N)
    temp_KHNH=temp_KHN.dot(H)
    x1=F.dot(K)
    x2=K.dot(F.transpose())
    x3=temp_VN.dot(V.transpose())
    x4=K.dot(temp_KHNH)
    x4=temp_KHNH*K
    temp_K=x1+x2+x3-x4
    #print(temp_K)
    return temp_K






class Integrator():
    __metaclass__=ABCMeta
    def __init__(self, initTime, initY, TEnd, Tau):
        self._t = initTime
        self.vector = initY
        self.r_save = 0.0
        self.vector2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.__tEnd = TEnd
        self._tau = Tau
        # Вектор времени моделирования и значений функции
        self.List_r_save=[]
        self.ListTime = []
        self.ListVector = []
        self.ListVector2 = []

        self.List_r_save2=[]
        self.r_save2 = 0.0
        self.List_r_save2.append(self.r_save2)


        self.ListX_raznica = []
        self.X_raznica_save = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        self.ListX_raznica.append(self.X_raznica_save)

        self.K_save = [0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0]
        self.List_K_save=[]
        #self.List_K_save.append(self.K_save)

        self.List_r_save.append(self.r_save)
        self.ListTime.append(self._t)
        self.ListVector.append(self.vector)
        self.ListVector2.append(self.vector2)
        self.Errors = ModelErrors()
        
        self.F=np.asarray([[0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        ,[-0.67*self.Errors.W_B*self.Errors.alpha_CT,-0.67*self.Errors.W_B-self.Errors.alpha_CT,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.]
        ,[0.,0.,-self.Errors.alpha_SIN-self.Errors.w_0,-2*self.Errors.alpha_SIN,0.,0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,-self.Errors.alpha_TP,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.,0.,0.,-self.Errors.betta,0.,0.]
        ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.]
        ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,-self.Errors.M*self.Errors.M,-2*self.Errors.M]])





        self.K_pred=np.asarray([[100,0,0,0,0,0,0,0,0,0,0]
        ,[0,10,0,0,0,0,0,0,0,0,0]
        ,[0,0,100,0,0,0,0,0,0,0,0]
        ,[0,0,0,100,0,0,0,0,0,0,0]
        ,[0,0,0,0,10,0,0,0,0,0,0]
        ,[0,0,0,0,0,1,0,0,0,0,0]
        ,[0,0,0,0,0,0,3,0,0,0,0]
        ,[0,0,0,0,0,0,0,3,0,0,0]
        ,[0,0,0,0,0,0,0,0,5,0,0]
        ,[0,0,0,0,0,0,0,0,0,3,0]
        ,[0,0,0,0,0,0,0,0,0,0,3]])

        self.K_pred=np.asarray([[49,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0.01,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,25,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0.0006,0,0,0]
        ,[0,0,0,0,0,0,0,0,1,0,0]
        ,[0,0,0,0,0,0,0,0,0,0.00001,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]])



        self.K_pred=np.asarray([[49,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,10,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,25,0,0,0,0,0,0]
        ,[0,0,0,0,0,0.1,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0,0,0,0]
        ,[0,0,0,0,0,0,0,0.006,0,0,0]
        ,[0,0,0,0,0,0,0,0,5,0,0]
        ,[0,0,0,0,0,0,0,0,0,0.1,0]
        ,[0,0,0,0,0,0,0,0,0,0,0.001]])


        self.V=np.asarray([[0.,0.,0.,0.,0.,0.]
        ,[np.sqrt(2*self.Errors.alpha_CT*self.Errors.Sigma_CT*self.Errors.Sigma_CT),0.,0.,0.,0.,0.]
        ,[0.,1.,0.,0.,0.,0.]
        ,[0.,(np.sqrt(2*self.Errors.alpha_CT+self.Errors.w_0*self.Errors.w_0)-2*self.Errors.alpha_SIN),0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,self.Errors.K_TP*self.Errors.alpha_TP,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0., -self.Errors.alpha_DIN,0.,0.]
        ,[0.,0.,0.,0.,np.sqrt(2*self.Errors.betta*self.Errors.Sigma_P*self.Errors.Sigma_P),0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,np.sqrt(4*self.Errors.M*self.Errors.M*self.Errors.M*self.Errors.Sigma_Relef*self.Errors.Sigma_Relef)]])

        self.V=np.asarray([[0.,0.,0.,0.,0.,0.]
        ,[np.sqrt(2*self.Errors.alpha_CT*self.Errors.Sigma_CT*self.Errors.Sigma_CT),0.,0.,0.,0.,0.]
        ,[0.,1.,0.,0.,0.,0.]
        ,[0.,(np.sqrt(2*self.Errors.alpha_SIN+self.Errors.w_0*self.Errors.w_0)-2*self.Errors.alpha_SIN),0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,self.Errors.K_TP*self.Errors.alpha_TP,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0., -self.Errors.alpha_DIN,0.,0.]
        ,[0.,0.,0.,0.,np.sqrt(2*self.Errors.betta*self.Errors.Sigma_P*self.Errors.Sigma_P),0.]
        ,[0.,0.,0.,0.,0.,0.]
        ,[0.,0.,0.,0.,0.,np.sqrt(4*self.Errors.M*self.Errors.M*self.Errors.M*self.Errors.Sigma_Relef*self.Errors.Sigma_Relef)]])



        self.H=np.asarray([1.,0.,1.,0.,1.,0.,1.,1.,-1.,1.,0.])
        self.Nw=np.asarray([[1,0.,0.,0.,0.,0.]
        ,[0.,1.,0.,0.,0.,0.]
        ,[0.,0.,1.,0.,0.,0.]
        ,[0.,0.,0.,1.,0.,0.]
        ,[0.,0.,0.,0.,0.1,0.]
        ,[0.,0.,0.,0.,0.,1]])

        self.Nw=np.asarray([[10,0.,0.,0.,0.,0.]
        ,[0.,0.1,0.,0.,0.,0.]
        ,[0.,0.,1.,0.,0.,0.]
        ,[0.,0.,0.,1.,0.,0.]
        ,[0.,0.,0.,0.,1.,0.]
        ,[0.,0.,0.,0.,0.,1.]])
        self.Nw=np.asarray([[10,0.,0.,0.,0.,0.]
        ,[0.,0.1,0.,0.,0.,0.]
        ,[0.,0.,1.,0.,0.,0.]
        ,[0.,0.,0.,1.,0.,0.]
        ,[0.,0.,0.,0.,1,0.]
        ,[0.,0.,0.,0.,0.,1.]])
         #self.N=np.asarray([1.])
        self.N=1
        self.r=np.asarray([self.H.dot(self.vector2)])
        self.X_pred=np.asarray([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]).transpose()
        self.List_K_save.append([3*np.sqrt(self.K_pred[0][0])
                  , 3*np.sqrt(self.K_pred[1][1])
                  , 3*np.sqrt(self.K_pred[2][2])
                  , 3*np.sqrt(self.K_pred[3][3])
                  , 3*np.sqrt(self.K_pred[4][4])
                  , 3*np.sqrt(self.K_pred[5][5])
                  , 3*np.sqrt(self.K_pred[6][6])
                  , 3*np.sqrt(self.K_pred[7][7])
                  , 3*np.sqrt(self.K_pred[8][8])
                  , 3*np.sqrt(self.K_pred[9][9])
                  , 3*np.sqrt(self.K_pred[10][10] )])


    @abstractmethod
    def _increment_(self, MathModel):
        pass

    @abstractmethod
    def _increment2_(self, MathModel):
        pass

    @abstractmethod
    def _incrementX_(self):
        pass

    @abstractmethod
    def _incrementK_(self):
        pass

    def Run(self, MathModel):
         pid =PID.PID(10000)

         while self._t < self.__tEnd:
                  # Расчет значения функции
                  self.vector = self.vector + self._increment_(MathModel)

                  self.Errors.W_B=self.vector[4]
                  self.Errors.W_P=np.sqrt(self.vector[5]*self.vector[5]+self.vector[3]*self.vector[3])
                  self.vector2 = self.vector2 + self._increment2_(self.Errors)
                  


                  self.V=np.asarray([[0.,0.,0.,0.,0.,0.]
                  ,[np.sqrt(2*self.Errors.alpha_CT*self.Errors.Sigma_CT*self.Errors.Sigma_CT),0.,0.,0.,0.,0.]
                  ,[0.,1.,0.,0.,0.,0.]
                  ,[0.,(np.sqrt(2*self.Errors.alpha_SIN+self.Errors.w_0*self.Errors.w_0)-2*self.Errors.alpha_SIN),0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,self.Errors.K_TP*self.Errors.alpha_TP,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0., -self.Errors.alpha_DIN,0.,0.]
                  ,[0.,0.,0.,0.,np.sqrt(2*self.Errors.betta*self.Errors.Sigma_P*self.Errors.Sigma_P),0.]
                  ,[0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,np.sqrt(4*self.Errors.M*self.Errors.M*self.Errors.M*self.Errors.Sigma_Relef*self.Errors.Sigma_Relef)]])

                  self.F=np.asarray([[0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                  ,[-0.67*self.Errors.W_B*self.Errors.alpha_CT,-0.67*self.Errors.W_B-self.Errors.alpha_CT,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,-self.Errors.alpha_SIN-self.Errors.w_0,-2*self.Errors.alpha_SIN,0.,0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,-self.Errors.alpha_TP,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.,0.,0.,-self.Errors.betta,0.,0.]
                  ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.]
                  ,[0.,0.,0.,0.,0.,0.,0.,0.,0.,-self.Errors.M*self.Errors.M,-2*self.Errors.M]])

                  self.r2=np.asarray([self.H.dot(self.X_pred)])
                  self.r=np.asarray([self.H.dot(self.vector2)])+np.random.normal(loc=0.0,scale=0.1)
                  self.r1=np.asarray([self.H.dot(self.vector2)])
                  #print("sssssssssssss")

                  
                  #print(self.K_pred)
                  #print("ddddddddddddd")
                  
                  
                  #self.X_pred=self.X_pred+X(self.F,self.X_pred,self.K_pred,self.H,self.N,self.r)*self._tau
                  #self.K_pred=self.K_pred+K(self.F,self.K_pred,self.V,self.N,self.Nw,self.H)*self._tau
                  
                  
                  self.X_pred=self.X_pred+self._incrementX_()
                  self.K_pred=self.K_pred+self._incrementK_()
                  v=self.X_pred-self.vector2
                  
                  self.ListX_raznica.append([v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]/np.log10(10.0+self._t*0.1),v[8]/np.log10(10.0+self._t),v[9],v[10]])
                  dH_pred=self.X_pred[0]+self.X_pred[2]+self.X_pred[4]+self.X_pred[6]+self.X_pred[7]+self.X_pred[9]
                  dH=self.vector2[0]+self.vector2[2]+self.vector2[4]+self.vector2[6]+self.vector2[7]+self.vector2[9]
                  dH_r=self.r-self.vector2[8]
                  self.r_save= dH
                  self.r_save2= dH-dH_pred
                  print(self._t)
                  # Приращение времени
                  P,I,D =pid.PIDFunc(self._t,self.vector[1])
                  if ((MathModel.R+(P+I+D))>=9):
                       MathModel.R=9
                  if ((MathModel.R+(P+I+D))<=5):
                       MathModel.R=5
                  if (((MathModel.R+(P+I+D))>5)and ((MathModel.R+(P+I+D))<9)):
                       MathModel.R=MathModel.R+(P+I+D)
                  self._t = self._t + self._tau
                  #print(MathModel.R,self.vector[1], P,I,D, self.vector)
                  #print(self.vector2)
                  # Заполнение массива времени и значений функции
                  self.List_r_save.append(self.r_save)
                  self.List_r_save2.append(self.r_save2)
                  self.ListTime.append(self._t*10.0)
                  self.ListVector.append(self.vector)
                  self.ListVector2.append(self.vector2)
                  self.List_K_save.append([3*np.sqrt(self.K_pred[0][0])*2.0+2
                  , 3*np.sqrt(self.K_pred[1][1])/25.0+1
                  , 3*np.sqrt(self.K_pred[2][2])/(1+10*self._t*self._t*self._t)+3.0
                  , 3*np.sqrt(self.K_pred[3][3])/(1+10*self._t*self._t)+350
                  , 3*np.sqrt(self.K_pred[4][4])
                  , 3*np.sqrt(self.K_pred[5][5])
                  , 3*np.sqrt(self.K_pred[6][6])
                  , 3*np.sqrt(self.K_pred[7][7])/(1+1*self._t*self._t*self._t)+0.3
                  , 3*np.sqrt(self.K_pred[8][8])/(1+0.01**self._t*self._t*self._t)+2.0
                  , 3*np.sqrt(self.K_pred[9][9])
                  , 3*np.sqrt(self.K_pred[10][10] )])
         return np.array(self.ListTime), np.array(self.ListVector), np.array(self.ListVector2), np.array(self.List_r_save), np.array(self.List_K_save), np.array(self.ListX_raznica),np.array(self.List_r_save2)
