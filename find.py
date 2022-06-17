from mpmath import *
import matplotlib.pyplot as plt
import numpy as np


s=0.5
def d0(s) :
   return zeta(s)

def d1(s):
   return zeta(s,derivative=1)
def f0(s):
   return d1(s)/d0(s)+(0.5*digamma((s)/2)-  1/(1-s)+   1/(s) -0.5*log(pi))*0
   
def f(t):
  return (zeta(s+j*t,derivative=1)*zeta(s-j*t)-zeta(s+j*t)*zeta(s-j*t,derivative=1)).imag

def findt(i):
   t0= zetazero(i-2).imag+0.0001
   t1= zetazero(i-1).imag-0.0001
   #print('==>ok0:',t0,t1,f0(0.5+j*t0).imag,f0(0.5+j*t1).imag)
   t=(t0+t1)/2
   #print('==>(t):',t,f0(0.5+j*t).imag,abs(t0-t1))
   while((abs(f0(0.5+j*t).imag)>=0.01)and(abs(t0-t1)>=0.00001)):
     t=(t0+t1)/2
     if((f0(0.5+j*t0).imag*f0(0.5+j*t).imag)<0):
      t1=t
      #print('==>left:',t0,t1,f0(0.5+j*t0).imag,f0(0.5+j*t1).imag)
     else:
      t0=t
      #print('==>right:',t0,t1,f0(0.5+j*t0).imag,f0(0.5+j*t1).imag)
   return (t0,t1)
def findz(i):
   if(i==1): t0= 0;t1=1
   else:
      if(i==2): t0= 2;t1=3
      else : t0,t1 = findt(i)
   try:
      r = findroot(lambda x:f0(0.5+j*x).imag, (t0,t1), solver ='illinois', verbose=False)
   except ValueError:
      r=0
   return r
def f2(s) :
   return zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
print(1/2*f2(1/2),1/2*zeta(1/2)*(1/2-1)*1/2*gamma(1/4)*pi**(-1/4))
print(1+euler**2-1/8*pi**2+2*stieltjes(1))
'''   
for i in range(150,250): #99
   print(findz(i))
   #print('---------------------------------------')
   #print('==>',findz(i),zetazero(i-2).imag,zetazero(i-1).imag)
'''
#plot([lambda t: f(t),lambda t: f0(0.5+j*t).imag],[14.3,15.25])