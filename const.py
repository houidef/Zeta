from mpmath import *
import matplotlib.pyplot as plt
import numpy as np
p0=zetazero(1)#0.5+j*14#
p_0=p0.conjugate()
t0=zetazero(1).imag
t1=t0-0.5
m0=1000
#q = [0]*m0
S1=0
for i in range(0,m0-1):
    h=zetazero(i+1)
    print(h.imag)
	#S1=S1+(1/(h*h.conjugate())**2)
	#if(i%10==0) : print(i/10,S1)
print(S1)

