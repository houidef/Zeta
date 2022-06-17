from mpmath import *
import matplotlib.pyplot as plt
import numpy as np
p0=zetazero(1)#0.5+j*14#
p_0=p0.conjugate()
t0=zetazero(1).imag
t1=t0-0.5
m0=10
q = [0]*m0
for i in range(0,m0-1):
    q[i]=zetazero(i+1).imag

def xi(s):
   return (zeta(s)*gamma(s/2))	
def f(s,p):
   return 2*(1-((1-s)*s)/(p**2+1/4))**0.5*exp(((1-s)*s)/(p**2+1/4))
def func(s):
   #return abs(2/((1-s/p_0)*(1-s/p0))*(s-1)*s*zeta(s)*zeta(1-s)) #((pi**(-s/2)*gamma(s/2)*zeta(s)*((s-1)*s/2)))) 
   #return abs((2/((1-s/p0)*(1-(1-s)/p0)*(1-s/p_0)*(1-(1-s)/p_0))*(s-1)*s)**0.5*zeta(s)*zeta(1-s))
   #return abs((1/((1-s/p0)*(1-(1-s)/p0)*(1-s/p_0)*(1-(1-s)/p_0))**2*(s-1)**1*s**1)*zeta(s)*zeta(1-s))
   return (2*f(s,q[0])*f(s,q[1])*f(s,q[2])*f(s,q[3])) #abs((1-s/p0)*(1-s/p_0))#(1-(1-s)/p0)*(1-(1-s)/p_0)*((s-1)**2*s**2))
#plot([lambda s: func(s+j*t0),lambda s: zeta(0.5+j*t0,derivative=1)*(s-0.5) ],[0.4,0.6])#14*(2*s-1)/(p_0*p0)+0*func(s+j*t1)
plot([lambda t: 0*func(0.5+j*t),lambda t:xi(0.5+j*t)],[-20,20])
'''
for i in range(1,10):
   x=zetazero(i)
   print(zeta(x,derivative=1),zeta(x,derivative=1)/x.imag)
'''