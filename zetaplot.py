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
Yz=[]
Yzj=[]
for i in x:
    j = ((abs(zeta(complex(0.5,i)))-1/abs(zeta(complex(0.5,i)))).real)
    k = 0*log((1-abs(zeta(complex(0.5,i)))).real**(2**0.5)) #0.6397#zt(850,complex(0.5,i),0.01)-zt(0.0001,complex(0.5,i),0.01)
    Yr.append(j.real)
    Yi.append(j.imag)
    #print(j.imag)
    Yz.append(k.real)
    Yzj.append(k.imag)
# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#f = findroot(myzeta,14)
#h0= zt(1,complex(0.5,f),0.1)-zt(1000,complex(0.5,f),0.1)
#h1=zt(1,complex(0.2,f),0.1)-zt(1000,complex(0.2,f),0.1)
#print('zr=',h0)
#print('zr=',h1/h0)
#print(f,zeta(complex(0.5,f)))
#print(f,zeta(complex(0.2,f))/zeta(complex(0.5,f)))
# plot the function
plt.plot(x,Yr, 'r')
plt.plot(x,Yi, 'b')
plt.plot(x,Yz, 'g')
plt.plot(x,Yzj, 'y')
#ax.scatter(f,0)
for k in range(1,10):
    ax.scatter(zetazero(k).imag,0)
# show the plot
plt.show()