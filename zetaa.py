from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

s=20
n=100
m0=10
t = np.linspace(-1,-2,n)
yr = [0]*n
yi = [0]*n
y = np.zeros((m0, n))
t0 = [0]*m0
p = [0]*m0
def theta(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   return 1/gamma(1+z/2) #(3*gamma(z/2))#1/(z1*z2)
'''  
for i,k in enumerate(t):
   r = log((zeta(s+j*k)**2)*gamma(s+j*k))+pi/2*k**0.999#0.995
   yr[i-1] = r.real
   yi[i-1] = 0*r.imag
'''
type = 'real'
for i in range(0,m0-1):
    p[i]=zetazero(i+1)
    if(type == 'imag'): t0[i]=zetazero(i+1).imag#(9.27)/(40)*i#
    else: t0[i]=0#zetazero(i+1).imag#zetazero(i+1).imag
for i,k in enumerate(t):
   #y[0][i] = abs(zeta((k+j*t0[1])/2-0.5)-zeta((k+j*t0[1])/2))
   for m in range(0,m0-1):
      s = k+j*t0[m]
      l = 0.5772156649015328606
      c = log(2*pi)-1-l/2
      if(type == 'imag'): y[m][i] = log((abs(zeta(k+j*t0[m]))-1/abs(zeta(k+j*tx[m]))).imag**(2**0.5)).real
      else : y[m][i] = abs((gamma(s/2)*zeta(s)*((s-1)*s/2)))# 0*abs(exp(c*s)*theta(s)*(1-s/tx[0])*exp(s/tx[0])*(1-s/tx[1])*exp(s/tx[1]) /(2*(s-1)*(1-s/0.5)**0*(1+s/0.5)**0))    #*(k+j*t0[m]-0.5-j*tx[1])**2*(k+j*t0[m]-0.5-j*tx[2])**2*(k+j*t0[m]-0.5-j*tx[3])**2)/(k+j*t0[m]-0.5)**0))
   #((log(abs(zeta(s)*zeta(1-s)*(s-0.5)**8/((s-0.5-j*tx[0])**2*(s-0.5-j*tx[1])**2*(s-0.5-j*tx[2])**2*(s-0.5-j*tx[3])**2)))))
   #((log(abs(zeta(k+j*t0[m])*zeta(1-k-j*t0[m])))))
   #abs((1-2**(1-(s/2-1/2)))*1/(s/2-1/2)*zeta(s/2-1/2)-(1-2**(1-s/2))*1/(s/2)*zeta(s/2))

poly = np.polyfit(t,y[1],1)
poly_y = sum(coeff*t**i for i,coeff in enumerate(reversed(poly)))
for coeff in enumerate(reversed(poly)):
    print(coeff)

figure, axes = plt.subplots()
for m in range(0,m0-1):
    if(type == 'imag'): axes.plot(t,(y[m])) #imaginaire
    else: plot(t,y[m]) #real
    print(m,y[m][99])
#axes.plot(t,poly_y,color="blue")
plt.show()