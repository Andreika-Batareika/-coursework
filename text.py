from numpy import*
from scipy.integrate import odeint
import matplotlib.pyplot as plt
g=9.81# ускорение свободного падения на земле в м/с2.
rv=1.29# плотность атмосферного воздуха в кг/м3.
rg=0.17# плотность гелия в кг/м3.
R=8# радиус оболочки ЛАЛВ в м.
b=0.000125# константа, связанная с плотностью воздуха в 1/м
a=6.5*10**-3# константа, связанная с температурой воздуха в К/м
c=0.4#коэффициент лобового сопротивления
mo=240#масса в кг
V=(4/3)*pi*R**3
rs=rg+mo/V# суммарная плотность материала ЛАЛВ, массы гелия, и нагрузки
p1=rv/rs# введенный параметр
p2=3*c/(8*R)# введенный параметр
T0=300
def fun(y, t):
         y1, y2= y
         return [y2,-g+g*p1*exp(-b*y1*T0/(T0-a*y1))-p1*p2*exp(-b*y1*T0/(T0-a*y1))*y2**2]
t =arange(0,1100,0.01)
y0 = [0.0,0.0]
[y1,y2]=odeint(fun, y0,t, full_output=False).T
plt.title("Характеристики  подъёма ЛАЛВ  \n Объём: %s м3. Масса : %s кг. \n Подъёмная сила: %s kН. "%(round(V,0),mo,round(0.001*g*rv*V,0)))
plt.plot(t/60,y2,label='Максимальная высота подъёма: %s км. \n Максимальная скорость: % s м/с .\n С учётом температуры воздуха'%(round(max(from numpy import*
from scipy.integrate import odeint
import matplotlib.pyplot as plt
g=9.81# ускорение свободного падения на земле в м/с2.
rv=1.29# плотность атмосферного воздуха в кг/м3.
rg=0.17# плотность гелия в кг/м3.
R=8# радиус оболочки ЛАЛВ в м.
b=0.000125# константа, связанная с плотностью воздуха в 1/м
a=6.5*10**-3# константа, связанная с температурой воздуха в К/м
c=0.4#коэффициент лобового сопротивления
mo=240#масса в кг
V=(4/3)*pi*R**3
rs=rg+mo/V# суммарная плотность материала ЛАЛВ, массы гелия, и нагрузки
p1=rv/rs# введенный параметр
p2=3*c/(8*R)# введенный параметр
T0=300
def fun(y, t):
         y1, y2= y
         return [y2,-g+g*p1*exp(-b*y1*T0/(T0-a*y1))-p1*p2*exp(-b*y1*T0/(T0-a*y1))*y2**2]
t =arange(0,1100,0.01)
y0 = [0.0,0.0]
[y1,y2]=odeint(fun, y0,t, full_output=False).T
plt.title("Характеристики  подъёма ЛАЛВ  \n Объём: %s м3. Масса : %s кг. \n Подъёмная сила: %s kН. "%(round(V,0),mo,round(0.001*g*rv*V,0)))
plt.plot(t/60,y2,label='Максимальная высота подъёма: %s км. \n Максимальная скорость: % s м/с .\n С учётом температуры воздуха'%(round(max(y1)/1000,2),round(max(y2),2)))
def fun(y, t):
         y1, y2= y
         return [y2,-g+g*p1*exp(-b*y1)-p1*p2*exp(-b*y1)*y2**2]
[y1,y2]=odeint(fun, y0,t, full_output=False).T
plt.plot(t/60,y2,label='Максимальная высота подъёма: %s км. \n Максимальная скорость: % s м/с \n Без учёта температуры воздуха'%(round(max(y1)/1000,2),round(max(y2),2)))
plt.ylabel('Высота в м')
plt.xlabel(' Время в минутах')
plt.legend(loc='best')
plt.grid(True)
plt.show()