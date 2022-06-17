from mpmath import *
mp.dps = 15; mp.pretty = True

a = [1/10,0.00007425,1.43E-7,1,1,1,1,1,1,1,1,1,1,1,1]
def f2(s):
   Sum=(s-0.5)
   for n in range(1,35):
      S=1
      for i in range(0,2*n+1):
          S=S*(s- i/(2*n))
      Sum=Sum+1/2*S
   return -Sum
   
myf = lambda s : (s-0.5)+s*(s-0.5)*(s-1)+s*(s-1/4)*(s-0.5)*(s-3/4)*(s-1)+s*(s-1/6)*(s-2/6)*(s-0.5)*(s-4/6)*(s-5/6)*(s-1)
xi0 = lambda s : zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
def d0(s) :
   return zeta(s)

def d1(s):
   return zeta(s,derivative=1)
xi1 = lambda s : d1(s)/d0(s)+0.5*digamma((s)/2)-  1/(1-s)+   1/(s) -0.5*log(pi)
c0   = log(2*pi)-1-euler/2-log(pi)/2
c1   = -3.71013755765504e-5
c2   = 1.4442560033602e-7
#xi = lambda s: (s - 1) * mp.pi ** (-0.5 * s) * mp.gamma(1 + 0.5 * s) * mp.zeta(s)
#tmp = taylor(lambda z: log(xi(z / (z - 1))), 0, 10)
tmp0 = [0]*25
tmp1 = [0]*25
i=1
s = 0.4
t0=zetazero(1).imag
b = 1.00340
print(log(pi)/2,c0,exp(1),zetazero(1).imag)
for k in [0.5]:#[1-0.0000001,0.0000001,0.50000000000001,1-0.50000000000001,5.0000001,1-5.0000001,1.5,1-1.5]:
  print('x=',k,'==========================================================================')
  #tmp0[i] = taylor(lambda z: f2(z), k, 10)
  #print(tmp0[i])
  
  tmp1[i] = taylor(lambda z: (d1(z)/d0(z)-  1/(1-z)), k, 10)
  print('d1/d0:',tmp1[i])
  tmp1[i] = taylor(lambda z: (0.5*digamma(z/2)+1/(z)-0.5*log(pi)), k, 10)
  print('gamma:',tmp1[i])
  tmp1[i] = taylor(lambda z: (d1(z)/d0(z)-  1/(1-z))+(0.5*digamma(z/2)+1/(z)-0.5*log(pi)), k, 10)
  print('xi   :',tmp1[i]); print('')
  
  tmp1[i] = taylor(lambda z: (d1(z)/d0(z)-  1/(1-z)**b), k, 10)
  print('d1/d0:',tmp1[i])
  tmp1[i] = taylor(lambda z: (0.5*digamma(z/2)+1/(z)**b-0.5*log(pi)), k, 10)
  print('gamma:',tmp1[i])
  tmp1[i] = taylor(lambda z: (d1(z)/d0(z)-  1/(1-z)**b)+(0.5*digamma(z/2)+1/(z)**b-0.5*log(pi)), k, 10)
  print('xi   :',tmp1[i]); print('')
  
  #tmp1[i] = taylor(lambda z:   xi1(z), k, 20)
  #print('xi  :',tmp1[i])
  i=i+1
'''
for t in range(0,5):
  print(tmp1[1][2*t+1]/2,tmp0[1][2*t+1],(tmp1[1][2*t+1]/tmp0[1][2*t+1])/2)
'''