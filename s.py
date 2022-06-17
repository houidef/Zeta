from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

m0=100
p = [0]*m0
z0 = zetazero(1)
t0 = zetazero(1).imag
t1=zetazero(2).imag
t2=zetazero(3).imag
t3=zetazero(4).imag
t4=zetazero(5).imag
t5=zetazero(6).imag
t6=zetazero(7).imag
t7=zetazero(8).imag
def theta(t):
   return arg(gamma(complex(0.25,0.5*t)))-0.5*t*log(pi)

def theta0(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   z3 = 1/pi**(z.conjugate()/2)*gamma(z.conjugate()/2)
   z4 = 1/pi**((1-z).conjugate()/2)*gamma((1-z).conjugate()/2)
   return log(z1)
def deriv0(s) :
   return (exp(j*theta(s.imag))*zeta(s))
def f_(s) :
   return (exp(j*theta(s.imag))*zeta(s)).imag

def deriv1(s):
   return (-exp(j*theta(s.imag))*zeta(s,derivative=1))

def deriv2(s):
   return (exp(j*theta(s.imag))*zeta(s,derivative=2))

def d0(s) :
   return zeta(s)

def d1(s):
   return zeta(s,derivative=1)

def d2(s):
   return zeta(s,derivative=2)
def f(s) :
   return gamma(s/2) #zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
def fx(s) :
   return digamma(s/2) #zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
def dd(s):
    h=0.000001
    return (f(s+h)-f(s))/h
def fs(s,t):
  return 0.5*(1*(3/((s+j*t)-(0.5+j*t0))+3/((s+j*t)-(0.5-j*t0))-1*digamma((s+j*t)/2)-2/(s+j*t)-2/(1-s+j*t)-0.0*(j*3.5-j*t))).imag 
def el(s,t0):
  return 1/(s-(0.5+j*t0))+1/(s-(0.5-j*t0))
def my(s):
   return( zeta(s)*((s-1)*s)*pi**(0.03516*s.imag+0.03516*s.real)*gamma(s/2)*pi**(-s/2))
   #return( zeta(s)*((s-1)*s)*(s.imag**2+0.5**2)**(0/2)*pi**((s)*(1-s)/2)*gamma(s/2)*pi**(-s/2))
def my0(s):
   return zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
def my_(s):
  return log(s)+log(s-1)-3.76810843775211 #(log(abs(t-1))-3.7+log(abs(t))-0.015*abs(t))
s=0.01
a=4#2

s=5.5
t=11
def func(s):
   return ((((d1(s)/d0(s)))-((-0.5*digamma((s)/2)+  1/(1-s)-   1/(s) ) +(el(s,t0) + el(s,t1)+ el(s,t2)+ el(s,t3)+ el(s,t4)+ el(s,t5)+ el(s,t6)+ 0*el(s,t7)-(log(2*pi)-1-euler/2-log(pi)/2)*s))))

def func0(s):
   return d1(s)/d0(s)+0.5*digamma((s)/2)-  1/(1-s)+   1/(s) 

print(pi/2,exp(-3.7-0.015),log(-(log(2*pi)-1-euler/2-log(pi)/2)))
#plot([lambda t: ((log((log(my0(30.6+t*j))))-my_(30.6+t*j))) ,lambda t: 0*log((log(my0(0.5+t*j))))],[-35,100])
t=t0
s=0
plot([lambda t: (log(s+t*j-1)+log(s+t*j))*2-10.6 ,lambda t: log(log((((my0(s+t*j))))/(exp(exp(my_(s+t*j))))))],[0,20])
#plot([lambda t:t-t0,lambda t: (log(s+t*j-1)+log(s+t*j))*2-10.6],[14,16])

print('my -x1:',my0(0.5+t0*j))
print('my -x2:',exp(exp(my_(0.5+t0*j))))
print('my -x2:',exp(exp(my_(0.5+1*j))))




'''
print('t0',t0)
print('t1',t1)
print('t2',t2)
print('t3',t3)
print('t4',t4)
print('t5',t5)
'''
#0.03516*j*t
#2/(t**2+1) ==> 0
#2/(2*t**2+0.5) ==> 0.5
#0==>1
#plot([lambda t: ((d1(s+j*t)/d0(s+j*t)+1/(1-s+j*t))).imag,lambda t: ((d1(s+j*t)/d0(s+j*t)+1/(1-s+j*t))).imag/((-0.5*digamma((s+j*t)/2)-1/(s+j*t)).imag +(1/((s+j*t)-(0.5+j*t0))+1/((s+j*t)-(0.5-j*t0))).imag) ,lambda t: (-0.5*digamma((s+j*t)/2)-1/(s+j*t)).imag, lambda t: (1/((s+j*t)-(0.5+j*t0))+1/((s+j*t)-(0.5-j*t0))).imag],[0,14])
s=2
t=0
print(func(s+j*t).real)
print(euler)
#==>plot([lambda t: euler,lambda t: (func(s+j*t).real)],[0,20])
#==>plot([lambda t: 0.0000*(t-3*pi)*(t+3*pi)*t,lambda t: func(s+j*t).imag],[-40,40])
#plot([lambda t: 0.00*(t-3*pi)**3,lambda t: (((d1(s+j*t)/d0(s+j*t))).imag-((-0.5*digamma((s+j*t)/2)+  1/(1-s-j*t)-   1/(s+j*t) +0.0*j*t ).imag +(el(s+j*t,t0) + el(s+j*t,t1)+ el(s+j*t,t2)+ el(s+j*t,t3)+ el(s+j*t,t4)+el(s+j*t,t5)+el(s+j*t,t6)+el(s+j*t,t7)+0.021*j*t).imag))],[-12,12]) #0.03516
#plot([lambda t: (my(s+j*t))],[0,15])
print('---------------------------------')#-  1/(1-s+j*t)-   1/(s+j*t)
'''
s=0.00001
print(d1(s)/d0(s)-1/(1-s)+0.5*digamma(s/2)+1/s-0.5*log(pi))

print(digamma(s/2)+1/s)
print(log(2*pi)-1/(1-s)-0.5*euler-0.5*log(pi))
s=1.00001
print(d1(s)/d0(s)-1/(1-s)+digamma(s/2)+2/s-0.5*log(pi))
print(d1(s)/d0(s)-1/(1-s)+0.5*digamma(s/2)+1/s-0.5*log(pi))
s=0.0001
'''
#plot([lambda t: (d1(s+j*t)/d0(s+j*t)).imag,lambda t: (1/((s+j*t)-(0.5+j*t0))).imag ],[0,t0-0.01])
#plot([lambda t: (d1(s+j*t)/d0(s+j*t)).imag,lambda t: fs(s,t) ],[0,14])
