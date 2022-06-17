from mpmath import *
import matplotlib.pyplot as plt
import numpy as np
def myzeta(x):
     return zeta(complex(0.2,x)).imag 
def zt0(xx,s,h):
     return 1/(1-s)*xx**(1-s)*h**(1-s) 
def zt(s):
     return zeta(s)-zeta(s) *(1-2**(1-s))
mp.dps = 25; mp.pretty = True
x = np.linspace(0,50,500)
Yr=[]
Yi=[]
for i in x:
    j = gamma(complex(0.5,i))
    Yr.append(j.real)
    Yi.append(j.imag)
# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


plt.plot(x,Yr, 'g')
plt.plot(x,Yi, 'y')
#ax.scatter(f,0)
for k in range(1,10):
    ax.scatter(zetazero(k).imag,0)
# show the plot
plt.show()