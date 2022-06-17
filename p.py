from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

def f1(s):
   return 1/(s-0.5-j*1)+1/(s-0.5+j*1)
def f2(s):
   Sum=(s-0.5)
   for n in range(1,10):
      S=1
      for i in range(0,2*n+1):
          S=S*(s- i/(2*n))
          #print(n,i,i/(2*n))
      Sum=Sum+0.5*S
   return -Sum
s=0.5
def d0(s) :
   return zeta(s)

def d1(s):
   return zeta(s,derivative=1)
def func0(s):
   return d1(s)/d0(s)+0.5*digamma((s)/2)-  1/(1-s)+   1/(s) -0.5*log(pi)
def func1(s):
   return  d1(s)/d0(s)
   
print(d1(0.5)/d0(0.5),digamma(1/4),+0.5*pi/2+1.5*ln(2)+0.5*euler+0.5*log(pi))
print(d1(0.0000)/d0(0.0000),ln(2*pi))
print(d1(1/4)/d0(1/4),(func0(1/4)*2-(log(2*pi)-1-euler/2-log(pi)/2))/2,func0(1/4)-(0.5*(-pi/2-4*ln(2)-(pi+ln(2**0.5+1)-ln(2**0.5-1))/2**0.5-euler)+4-4/3-0.5*log(pi)))
#plot([lambda t: 0*f1(s+j*t),lambda t: f2(s+j*t)],[-2,2])
plot([lambda t: func0(s+j*t).conjugate(),lambda t: func0(s-j*t)],[15,18])
print(func0(0.5+j*17.8825820769367))