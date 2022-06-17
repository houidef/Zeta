from mpmath import *
mp.dps = 15; mp.pretty = True
a0  = 1#-0.00577508738538607
a1  = 1# 0.023104993115419
a2  = 50#-1.85862996426348e-5
a3  = 0# 4.80579771336578e-8
a4  = 0#-1.65757920063248e-10 
a5  = 0# 6.4273283012332e-13  
#myf = lambda s : (s-0.5)+s*(s-0.5)*(s-1)+s*(s-1/4)*(s-0.5)*(s-3/4)*(s-1)+s*(s-1/6)*(s-2/6)*(s-0.5)*(s-4/6)*(s-5/6)*(s-1)
myf_ = lambda s : exp(a0)*(1 +  a1*(s-0.5)**2 +  (1/2*a1**2+a2)*(s-0.5)**4 + (1/6*a1**3+a1*a2)*(s-0.5)**6 + (1/24*a1**4+a2**2/2+a1**2*a2/2)*(s-0.5)**8+
(1/120*a1**5+(1/6*a1**3+1/2*a1))*(s-0.5)**10)
test = lambda s : a0+a1*(s-0.5)**2+a2*(s-0.5)**4+a3*(s-0.5)**6+a4*(s-0.5)**8+a5*(s-0.5)**10
myf0 = lambda s : a0+a1*(s-0.5)**2+a2*(s-0.5)**4+a3*(s-0.5)**6+a4*(s-0.5)**8+a5*(s-0.5)**10
myf1 = lambda s : 2*a1*(s-0.5)+4*a2*(s-0.5)**3+6*a3*(s-0.5)**5+8*a4*(s-0.5)**7+10*a5*(s-0.5)**9
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
tmp0 = [0]*22
tmp1 = [0]*22
i=1
s = 0.4
t0=zetazero(1).imag
print(log(pi)/2,c0,exp(1))
for k in [0.5000000000000]:#[1-0.0000001,0.0000001,0.50000000000001,1-0.50000000000001,5.0000001,1-5.0000001,1.5,1-1.5]:
  print('x=',k,'==========================================================================')
  print('--------xi0------------')
  tmp0[i] = taylor(lambda z: myf_(z), k, 10)
  print(tmp0[i])
  tmp1[i] = taylor(lambda z: exp(test(z)), k, 10)
  print(tmp1[i])
  '''
  tmp0[i] = taylor(lambda z: (xi0(z)), k, 10)
  print(tmp0[i])
  print('--------xi1------------')
  tmp1[i] = taylor(lambda z: myf1(z), k, 10)
  print(tmp1[i])
  tmp1[i] = taylor(lambda z: xi1(z), k, 10)#/((z-0.5)*(z-1/6)*(z-2/6)*(z-4/6)*(z-5/6)*(z-1)*z)
  print(tmp1[i])
  '''
  i=i+1

for t in range(0,6):
  print((tmp1[1][2*t]-tmp0[1][2*t])/exp(1))