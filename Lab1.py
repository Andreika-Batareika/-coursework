import numpy as np
import matplotlib.pyplot as modulPlot
from runge import Runge4
from EkranoplanModel import EkranoplanModel


E_Mod = EkranoplanModel()
# начальные условия
startPosition = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 2.0])
#интегратор
integrator = Runge4(0, startPosition, 100.0, 0.001)
#моделирование
time, result, resultErrors, dr,K11,x_raznica,resultErrors2 = integrator.Run(E_Mod)

print(x_raznica)

# График X экраноплана
modulPlot.subplot(231)
modulPlot.plot(time, result[:,0], color='r', linewidth=1)
modulPlot.title("X(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X, м')
modulPlot.grid(True)

# График Y экраноплана
modulPlot.subplot(232)
modulPlot.plot(time, result[:,1] , color='r', linewidth=1)
modulPlot.title("Y(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('Y, м')
modulPlot.grid(True)

# График Z экраноплана
modulPlot.subplot(233)
modulPlot.plot(time, result[:,2], color='r', linewidth=1)
modulPlot.title("Z(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X, м')
modulPlot.grid(True)

# График Vx экраноплана
modulPlot.subplot(234)
modulPlot.plot(time, result[:,3], color='r', linewidth=1)
modulPlot.title(" Vx(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('Vx, м/c')
modulPlot.grid(True)

# График Vy экраноплана
modulPlot.subplot(235)
modulPlot.plot(time, result[:,4], color='r', linewidth=1)
modulPlot.title(" Vy(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('Vy, м/c')
modulPlot.grid(True)

# График Vz экраноплана
modulPlot.subplot(236)
modulPlot.plot(time, result[:,5], color='r', linewidth=1)
modulPlot.title(" Vz(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('Vz, м/c')
modulPlot.grid(True)
modulPlot.show()
###########################################################3
modulPlot.subplot(231)
modulPlot.plot(time, resultErrors[:,0], color='r', linewidth=1)
modulPlot.title("X[0](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[0], м')
modulPlot.grid(True)

modulPlot.subplot(232)
modulPlot.plot(time, resultErrors[:,1], color='r', linewidth=1)
modulPlot.title("X[1](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[1], м')
modulPlot.grid(True)

modulPlot.subplot(233)
modulPlot.plot(time, resultErrors[:,2], color='r', linewidth=1)
modulPlot.title("X[2](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[2], м')
modulPlot.grid(True)

modulPlot.subplot(234)
modulPlot.plot(time, resultErrors[:,3], color='r', linewidth=1)
modulPlot.title("X[3](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[3], м')
modulPlot.grid(True)

modulPlot.subplot(235)
modulPlot.plot(time, resultErrors[:,4], color='r', linewidth=1)
modulPlot.title("X[4](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[4], м')
modulPlot.grid(True)

modulPlot.subplot(236)
modulPlot.plot(time, resultErrors[:,5], color='r', linewidth=1)
modulPlot.title("X[5](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[5], м')
modulPlot.grid(True)
modulPlot.show()

modulPlot.subplot(231)
modulPlot.plot(time, resultErrors[:,6], color='r', linewidth=1)
modulPlot.title("X[6](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[6], м')
modulPlot.grid(True)

modulPlot.subplot(232)
modulPlot.plot(time, resultErrors[:,7], color='r', linewidth=1)
modulPlot.title("X[7](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[7], м')
modulPlot.grid(True)


modulPlot.subplot(233)
modulPlot.plot(time, resultErrors[:,8], color='r', linewidth=1)
modulPlot.title("X[8](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[8], м')
modulPlot.grid(True)


modulPlot.subplot(234)
modulPlot.plot(time, resultErrors[:,9], color='r', linewidth=1)
modulPlot.title("X[9](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[9], м')
modulPlot.grid(True)





modulPlot.subplot(2,3,5)
modulPlot.plot(time, resultErrors[:,10], color='r', linewidth=1)
modulPlot.title("X[10](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('X[10], м')
modulPlot.grid(True)
modulPlot.show()


modulPlot.subplot(1,1,1)
modulPlot.plot(time, dr, color='r', linewidth=1)
modulPlot.plot(time, resultErrors2/2.0, color='b', linewidth=1)



modulPlot.title("dr(t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dr, м')
modulPlot.grid(True)
modulPlot.show()



modulPlot.subplot(231)
modulPlot.plot(time, K11[:,0], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,0], color='b', linewidth=1)

modulPlot.plot(time, -1*K11[:,0], color='r', linewidth=1)
modulPlot.title("K[0](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[0], м')
modulPlot.grid(True)

modulPlot.subplot(232)
modulPlot.plot(time, K11[:,1], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,1], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,1], color='r', linewidth=1)
modulPlot.title("K[1](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[1], м')
modulPlot.grid(True)


modulPlot.subplot(233)
modulPlot.plot(time, K11[:,2], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,2], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,2], color='r', linewidth=1)
modulPlot.title("K[2](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[2], м')
modulPlot.grid(True)


modulPlot.subplot(234)
modulPlot.plot(time, K11[:,3], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,3], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,3], color='r', linewidth=1)
modulPlot.title("K[3](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[3], м')
modulPlot.grid(True)





modulPlot.subplot(2,3,5)
modulPlot.plot(time, K11[:,4], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,4], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,4], color='r', linewidth=1)
modulPlot.title("K[4](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('Xd[4], м')
modulPlot.grid(True)
modulPlot.show()







modulPlot.subplot(231)
modulPlot.plot(time, K11[:,5], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,5], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,5], color='r', linewidth=1)
modulPlot.title("K[5](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[5], м')
modulPlot.grid(True)

modulPlot.subplot(232)
modulPlot.plot(time, K11[:,6], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,6], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,6], color='r', linewidth=1)
modulPlot.title("K[6](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[6], м')
modulPlot.grid(True)


modulPlot.subplot(233)
modulPlot.plot(time, K11[:,7], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,7]/4.0, color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,7], color='r', linewidth=1)
modulPlot.title("K[7](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[7], м')
modulPlot.grid(True)


modulPlot.subplot(234)
modulPlot.plot(time, K11[:,8], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,8], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,8], color='r', linewidth=1)
modulPlot.title("K[8](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[8], м')
modulPlot.grid(True)





modulPlot.subplot(2,3,5)
modulPlot.plot(time, K11[:,9], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,9], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,9], color='r', linewidth=1)
modulPlot.title("K[9](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[9], м')
modulPlot.grid(True)


modulPlot.subplot(2,3,6)
modulPlot.plot(time, K11[:,10], color='r', linewidth=1)
modulPlot.plot(time, x_raznica[:,10], color='b', linewidth=1)
modulPlot.plot(time, -1*K11[:,10], color='r', linewidth=1)
modulPlot.title("K[10](t)")
modulPlot.xlabel('t, с')
modulPlot.ylabel('dX[10], м')
modulPlot.grid(True)
modulPlot.show()