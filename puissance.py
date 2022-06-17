from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

def puis(x):
   t=0
   for i in range(1,1000):
     t=t+(x/i)**i
   return t

def puis1(x):
   t=0
   for i in range(1,1000):
     t=t+(x/i)**(i-1)
   return t
  
 
x = np.linspace(1,20,500)
#y=[]
#for c in x:
#   y.append(complex(0,c))
y=x*j
yr=[]
yi=[]
for c in x:
    h = gamma(0.5+j*c)
    yr.append(h.real)
    yi.append(h.imag)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
	
#plt.plot(x,puis(y).imag, 'g')
#plt.plot(x,yr, 'r')
#plt.plot(x,yi, 'b')
for i in range(0,1615):
   print(puis(i))
   #print(i,puis(i*0.01*j).real*i*0.01,puis(i*0.01*j).imag*i*0.01,np.exp(i*0.01*j).real,np.exp(i*0.01*j).imag)

#plt.show()