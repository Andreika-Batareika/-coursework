from Integrator import Integrator,X,K


#интегратор Рунге Кутты 4 порядка
class Runge4(Integrator):

    def __init__(self,initTime, initY, TEnd, Tau):
        Integrator.__init__(self,initTime, initY, TEnd, Tau)


    def _increment_(self, MathModel):
        k1 = MathModel.RightPart(self._t, self.vector)
        k2 = MathModel.RightPart(self._t + 0.5 * self._tau, self.vector + 0.5 * k1 * self._tau)
        k3 = MathModel.RightPart(self._t + 0.5 * self._tau, self.vector + 0.5 * k2 * self._tau)
        k4 = MathModel.RightPart(self._t + self._tau, self.vector + k3 * self._tau)
        out = (k1 + 2 * k2 + 2 * k3 + k4) * self._tau / 6
        return out

    def _increment2_(self, MathModel2):
        k1 = MathModel2.RightPart(self._t, self.vector2)
        k2 = MathModel2.RightPart(self._t + 0.5 * self._tau, self.vector2 + 0.5 * k1 * self._tau)
        k3 = MathModel2.RightPart(self._t + 0.5 * self._tau, self.vector2 + 0.5 * k2 * self._tau)
        k4 = MathModel2.RightPart(self._t + self._tau, self.vector2 + k3 * self._tau)
        out = (k1 + 2 * k2 + 2 * k3 + k4) * self._tau / 6
        return out

    def _incrementK_(self):
        F=self.F 
        N=self.N
        K_pred=self.K_pred
        Nw=self.Nw
        H=self.H
        V=self.V
        k1 = K(F,K_pred,V,N,Nw,H)
        k2 = K(F,K_pred+ 0.5 * k1 * self._tau,V,N,Nw,H)
        k3 = K(F,K_pred+ 0.5 * k2 * self._tau,V,N,Nw,H)
        k4 = K(F,K_pred+ k3 * self._tau,V,N,Nw,H)
        out = (k1 + 2.0 * k2 + 2.0 * k3 + k4) * self._tau / 6.0
        return out

    def _incrementX_(self):
        F=self.F 
        N=self.N
        K_pred=self.K_pred
        Nw=self.Nw
        H=self.H
        V=self.V
        X_pred=self.X_pred
        r=self.r
        k1 = X(F,X_pred,K_pred,H,N,r)
        k2 = X(F,X_pred + 0.5 * k1 * self._tau,K_pred,H,N,r)
        k3 = X(F,X_pred + 0.5 * k2 * self._tau,K_pred,H,N,r)
        k4 = X(F,X_pred +  k3 * self._tau,K_pred,H,N,r)
        out = (k1 + 2.0 * k2 + 2.0 * k3 + k4) * self._tau / 6.0
        return out



#интегратор Рунге Кутты 4 порядка
class Euler(Integrator):

    def __init__(self,initTime, initY, TEnd, Tau):
        Integrator.__init__(self,initTime, initY, TEnd, Tau)


    def _increment_(self, MathModel):
        k1= MathModel.RightPart(self._t, self.vector)*self._tau
        out = k1
        return out
