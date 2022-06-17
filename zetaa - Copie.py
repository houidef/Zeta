from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

s=20
n=100
m0=20
t = np.linspace(0.1,20,n)
yr = [0]*n
yi = [0]*n
y = np.zeros((m0, n))
t0 = [0]*m0
'''
for i,k in enumerate(t):
   r = log((zeta(s+j*k)**2)*gamma(s+j*k))+pi/2*k**0.999#0.995
   yr[i-1] = r.real
   yi[i-1] = 0*r.imag
'''
type = 'real'
for i in range(0,m0-1):
    if(type == 'imag'): t0[i]=i#zetazero(i+1).imag#(9.27)/(40)*i#
    else: t0[i]=0+(1)/(20)*i
for i,k in enumerate(t):
   for m in range(0,m0-1):
      if(type == 'imag'): y[m][i] = log((abs(zeta(k+j*t0[m]))-1/abs(zeta(k+j*t0[m]))).imag**(2**0.5)).real
      else : y[m][i] = (((log((abs(zeta(k+j*t0[m]))-1/abs(zeta(k+j*t0[m]))).real**(2**0.5)))).real+k)
   #print(y[1][i-1],y[2][i-1],y[3][i-1],y[4][i-1],y[5][i-1])

poly = np.polyfit(t,y[1],1)
poly_y = sum(coeff*t**i for i,coeff in enumerate(reversed(poly)))
for coeff in enumerate(reversed(poly)):
    print(coeff)

figure, axes = plt.subplots()
for m in range(1,m0-1):
    if(type == 'imag'): axes.plot(t,(y[m])) #imaginaire
    else: axes.plot(t,y[m]) #real
    print(m,y[m][99])
#axes.plot(t,poly_y,color="blue")
plt.show()