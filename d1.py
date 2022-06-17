from mpmath import *
import matplotlib.pyplot as plt
import numpy as np
p0=0.5+j*14#zetazero(1)
p_0=p0.conjugate()
t0=13.5#zetazero(1).imag
t1=t0-0.5

def func(s):
   #return abs(2/((1-s/p_0)*(1-s/p0))*(s-1)*s*zeta(s)*zeta(1-s)) #((pi**(-s/2)*gamma(s/2)*zeta(s)*((s-1)*s/2)))) 
   #return abs((2/((1-s/p0)*(1-(1-s)/p0)*(1-s/p_0)*(1-(1-s)/p_0))*(s-1)*s)**0.5*zeta(s)*zeta(1-s))
   #return abs((1/((1-s/p0)*(1-(1-s)/p0)*(1-s/p_0)*(1-(1-s)/p_0))**2*(s-1)**1*s**1)*zeta(s)*zeta(1-s))
   return abs((1-s/p0)*(1-s/p_0))#(1-(1-s)/p0)*(1-(1-s)/p_0)*((s-1)**2*s**2))
plot([lambda s: func(s+j*t0),lambda s: t0*(2*s-1)/(p_0*p0) ],[-5.9,5.9])#14*(2*s-1)/(p_0*p0)+0*func(s+j*t1)
#plot([lambda t: func(0.5+j*t),lambda t: func(0.6+j*t)],[13.99,14.01])
print(func(0.5))