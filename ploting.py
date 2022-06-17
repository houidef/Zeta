from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

s=20
#plot( lambda x: 2*log(abs(zeta(s+j*zetazero(x)))/abs(zeta(0.5+j*zetazero(x)))),scatter = [0,10])
n=40
x = range(1,n+1)
y = [0] * n
z = [0] * n
for k in x:
    z[k-1]=zetazero(k).imag
for k in x:
    y[k-1]=(1-abs(zeta(s+j*z[k-1]))**2)/(1-abs(zeta(0.51+j*z[k-1]))**2)


print(y[0],exp(y[0]),exp(y[1]),exp(y[19]))
figure, axes = plt.subplots()

d1 = [0] * n
d2 = [0] * n
for k in x:
   d1[k-1] = 1
   d2[k-1]=  (k**(abs(0.51-s)+j*z[k-1])).real

axes.scatter(x, y)

#poly = np.polyfit(x,y,1)

#poly_x = np.linspace(0,20,40)
#poly_y = sum(coeff*poly_x**i for i,coeff in enumerate(reversed(poly)))
#for coeff in enumerate(reversed(poly)):
#    print(coeff)

#axes.plot(poly_x,poly_y,color="blue")
axes.plot(x,d1,color="green")
axes.plot(x,d2,color="y")
figure.savefig("polyfit_example.png")
plt.show()