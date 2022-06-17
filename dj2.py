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
def my0(s):
   return zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)

cst = 3.76810843775211
c0   = log(2*pi)-1-euler/2-log(pi)/2
c1   = log(2*pi)-1-euler/2
def my_(s):
  return log(s)+log(s-1)- cst#(log(abs(t-1))-3.7+log(abs(t))-0.015*abs(t))
def my_i(s):
  #return -s*(s-1)*c0#s*(s-1)*exp(-cst)+(-s*(s-1))**2*exp(-11)*0
  return -s*(s-1)*c0 #c0*s-s**2*c0

def func(s):
   return ((((d1(s)/d0(s)))-((-0.5*digamma((s)/2)+  1/(1-s)-   1/(s) ) +(el(s,t0) + el(s,t1)+ el(s,t2)+ el(s,t3)+ el(s,t4)+ el(s,t5)+ el(s,t6)+ 0*el(s,t7)-(log(2*pi)-1-euler/2-log(pi)/2)*s))))

def func0(s):
   return d1(s)/d0(s)+0.5*digamma((s)/2)-  1/(1-s)+   1/(s) -0.5*log(pi)

def func1(s):
   h=0.0000001
   return (func0(s+h)-func0(s))/h
def f_(s):
   return -s*(s-0.25)*(s-0.5)*(s-0.75)*(s-1)
   
#print(pi/2,exp(-3.7-0.015),log(-(log(2*pi)-1-euler/2-log(pi)/2)))
#plot([lambda t: ((log((log(my0(30.6+t*j))))-my_(30.6+t*j))) ,lambda t: 0*log((log(my0(0.5+t*j))))],[-35,100])
t=0
s=1.500001
#-----------------------------------------------------------------
#plot([lambda  s: exp(my_i(s+t*j)) ,lambda s:(my0(s+t*j))],[0,15])
#------------------------------------------------------------------
#plot([lambda  t: (log(s+t*j)+log(s+t*j-1))*2-11 ,lambda t:log(log((my0(s+t*j))/exp(my_i(s+t*j))))-((log(s+t*j)+log(s+t*j-1))*2-11) ],[0,40])

#plot([lambda  s: (s+t*j)*(s+t*j-1)*0.00,lambda s: log(-log(my0(s+t*j))+my_i(s+t*j))/log((s+t*j)*(s+t*j-1))],[0,10])
#plot([lambda  t: (ln((s+j*t-0.5-j*t0)*(s+j*t-0.5+j*t0))**0.60+1/(s+j*t-0.5-j*t0)+1/(s+j*t-0.5+j*t0)-2)*0,lambda  t: func0(s+j*t)],[14,15])
h=0.0000001
#preuve1:-------------------------------------------------------------------------------
#plot([lambda s:2*c0 ,lambda t: (func0(s+j*t)-func0(s+j*(t+h)))/(j*h)],[-0.1,0.1]) #for t
#plot([lambda s:2*c0 ,lambda s: (func0(s+j*t)-func0(s+h+j*(t)))/(h)],[-0.1,0.1]) #for s
#---------------------------------------------------------------------------------------
#plot([lambda s: (-2*(s+j*t)*c0+c0)*0 ,lambda s: ((func0(s+j*t)-(-2*(s+j*t)*c0+c0)))/ln(1-s+j*t/(0.5+j*t0))],[0,50])
#plot([lambda s: (s+j*t)*(s+j*t-0.5)*(s+j*t-1)*0.000068 ,lambda s: (-((func0(s+j*t)-(-2*(s+j*t)*c0+c0))))],[-1,3])
plot([lambda s: 0,lambda s: (-((func0(s+j*t)-c0*(-2*(s+j*t)+1)+(s+j*t)*(s+j*t-0.5)*(s+j*t-1)*0.0000742555+2*1.44E-7*f_(s+j*t))))],[-0.25,1.25])
