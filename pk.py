from mpmath import *
mp.dps = 15; mp.pretty = True

d0 = -2.30954781990365e-7
c0= log(2*pi)-1-euler/2-log(pi)/2
#c0 = -0.0231049931151685
c1 = -3.71725992849473e-5
c2 = 1.44173931400035e-7
c3 = -6.63031680250115e-10
c4 = 3.21366415060963e-12


def mf(s):
   Sum=(s-0.5)
   for n in range(1,2):
      S=1
      for i in range(0,2*n+1):
          S=S*(s- i/(2*n))
      Sum=Sum+1/2*S
   return -Sum
   
b = [0.99424155637662825, 0.022971944315145439, 0.00024690403614063605, 1.6647109627710542e-6, 7.9844531026882735e-9, 2.923205152022192e-11, 8.549080091073689e-14, 2.0619252269236015e-16, 4.1995396162990572e-19, 7.3562821910015406e-22, 1.1245715174645487e-24]
def f0_(s):
   S=0
   for i in range(0,10):
      S = S+b[i]*(s-0.5)**(2*i)
f0  = lambda s: zeta(s,derivative=1)/zeta(s)+0.5*digamma((s)/2)-  1/(1-s)+   1/(s) -0.5*log(pi)
f0_ = lambda s: 0.0462099862308379*(s-0.5)-7.43451985705394e-5*(s-0.5)**3+2.88347862801947e-7*(s-0.5)**5-1.32606336050598e-9*(s-0.5)**7+6.4273283012332e-12*(s-0.5)**8
f0__= lambda s: c0*(-2*s+1)-s*(2*s-1)*(s-1)*3.71278E-05-1.44E-7*(-s*(s-0.25)*(2*s-1)*(s-0.75)*(s-1))+1/2*6.63031680250115e-10*(-s*(s-1/6)*(s-2/6)*(2*s-1)*(s-4/6)*(s-5/6)*(s-1))

a = [0.994242, 0.0229719, 0.000246904, 1.66471e-6, 7.98445e-9, 2.92321e-11, 8.54908e-14, 2.06193e-16, 4.19954e-19, 7.35628e-22, 1.12457e-24, 1.51835e-27, 1.82867e-30, 1.98122e-33, 1.94494e-36, 1.74092e-39, 1.4287e-42, 1.08019e-45, 7.55691e-49, 4.91082e-52, 2.97475e-55]
def f1_(s):
   S=0
   for i in range(0,15):
      S = S+a[i]*(s-0.5)**(2*i)
   return S
f1x  = lambda s: 2*s*zeta(s)*gamma(s)*(2*pi)**(-s)
f2x  = lambda s: -s*zeta(1-s)
f1  = lambda s: zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2)
#f1_ = lambda s: 0.0462099862308379*(s-0.5)-7.43451985705394e-5*(s-0.5)**3+2.88347862801947e-7*(s-0.5)**5-1.32606336050598e-9*(s-0.5)**7+6.4273283012332e-12*(s-0.5)**8
f1__= lambda s: exp(-c0*s*(s-1)-1.85503644842093e-5*(s*(s-1))**2+4.79850774227238e-8*(s*(s-1))**3-1.64956984796612e-10*(s*(s-1))**4+6.38763596494593e-13*(s*(s-1))**5-2.64615575248503e-15*(s*(s-1))**6)
f2__= lambda s: (1-c0*s*(s-1)+0.000248155522911367*(s*(s-1))**2+1.6728051129419e-6/2*(s*(s-1))**3+792091585683e-9/4*(s*(s-1))**4+0*6.71097376669968e-10*(s*(s-1))**5)

t__= lambda s: exp(10+9*s*(s-1)+8*(s*(s-1))**2+7*(s*(s-1))**3+6*(s*(s-1))**4+  5*(s*(s-1))**5+4*(s*(s-1))**6+0*(s*(s-1))**7)

t=zetazero(1).imag
s=0.5
h= 0.0000000001
"""
for k in range(1,10):
  print(f2x(k+h),f1x(k+h),bernoulli(k))
"""
"""
temp = taylor(lambda s: (s)/(exp(s)-1), 0+h, 10)
print([temp[i]*factorial(i) for i in range(0,10)])

temp = taylor(lambda s: (-s)/(exp(-s)-1), 0+h, 10)
print([temp[i]*factorial(i) for i in range(0,10)])
"""
y = lambda s : 20*s*(s-1)+6*s*(s-1)*s*(s-1/2)
def u(t,s):
 return (t)**s/(exp(t**s)-1)
print(ln(f1(1+h))*(1/h))
temp = taylor(lambda s: ln((zeta(s))*(zeta(s).conjugate())), (1+h),20)
#temp = taylor(lambda s: zeta(s)*(s-1)*s*gamma(s/2)*pi**(-s/2), h,20)
#temp = taylor(lambda s: ln(f1(s)), 0.5+0.5**0.5+0.5**0.75,2) #0.707175
#temp = taylor(lambda s: 0*(zeta(1-s)+1/(s))+1*(zeta(0.5+s)+1/(1-(0.5+s))), -0.5**0.5,20)
#temp = taylor(lambda s: 0*(zeta(1-s)*(s))+1*(zeta(s)*(1-s)), 5, 20)
#temp = taylor(lambda s: 0*(zeta(1-s))+1*(zeta(s)), 0.5, 20)
#temp = taylor(lambda s: digamma(s), 6, 20)
print([temp[i] for i in range(0,20)])
for i in range(0,20):
  print(temp[i])

"""
for k in range(1,2):
  temp = taylor(lambda s: (f2x(s))**k, 2, 10)
  #print(abs(temp[3]))	
  #print([temp[i] for i in range(0,10)])
for i in range(0,10):
  print(temp[i])
print('-------------')
"""
"""
for k in range(-5,5):
  temp = taylor(lambda s: (f1x(s))**k, 1+h, 5)
  print([temp[i] for i in range(0,5)])
"""
#temp = taylor(lambda s: ln(t__(s)), 0, 10)
#print([temp[i] for i in range(0,10)])
'''
for k in [1,0.5]:
  print(k,'************************************************************')
  k=k+h
  temp = taylor(lambda s: (f0(s)), k, 12)
  print([temp[i] for i in range(0,12)])
  
  print('==========')
  k0 = temp[0]
  k1 = temp[1]
  k2 = temp[2]-   k1
  k3 = temp[3]- 2*k2
  k4 = temp[4]-(3*k3+k2)
  k5 = temp[5]-(4*k4+3*k3)
  k6 = temp[6]-(5*k5+6*k4+k3)
  k7 = temp[7]-(6*k6+10*k5+4*k4)
  print(k0,k1,k2,k3,k4,k5,k6,k7)
  print('==========')
  k0 =  temp[0]
  k1 = -temp[1]
  k2 =  temp[2]-   k1
  k3 = -temp[3]- 2*k2
  k4 =  temp[4]-(3*k3+k2)
  k5 = -temp[5]-(4*k4+3*k3)
  k6 =  temp[6]-(5*k5+6*k4+k3)
  k7 = -temp[7]-(6*k6+10*k5+4*k4)
  print(k0,k1,k2,k3,k4,k5,k6,k7)
  print('==========')
  
  print('==========')
  k0 =  temp[0] + 1/4 * temp[2] +    1/4**2 * temp[4] +1 * 1/4**3 * temp[6] +1/4**4     * temp[8]  + 1/4**5      * temp[10]+1/4**6 * temp[12]
  k1 =  temp[2] + 2/4 * temp[4] + 3* 1/4**2 * temp[6] +4 * 1/4**3 * temp[8] +0.01953125 * temp[10] + 0.005859375 * temp[12]
  k2 =  temp[4] + 3/4 * temp[6] + 6* 1/4**2 * temp[8] +10* 1/4**3 * temp[10]+0.05859375 * temp[12]
  k3 =  temp[6] + 4/4 * temp[8] + 10*1/4**2 * temp[10]+20* 1/4**3 * temp[12] 
  k4 =  temp[8] + 5/4 * temp[10]+ 15*1/4**2 * temp[12]
  k5 =  temp[10]+ 6/4 * temp[12] 
  k6 =  temp[12]
  k7 = 0
  print(k0,k1,k2,k3,k4,k5,k6,k7)#(9-k1)/temp[12]
  print('==========')
'''
"""
temp = taylor(lambda s: (t__(s)), 1, 10)
print([temp[i]/temp[0] for i in range(0,10)])
print('==========')
for i in range(1,7):
  temp[i] = temp[i]/temp[0]
k0 = ln(temp[0])
k1 = temp[1]
k2 = temp[2]-0.5*(temp[1]+1)**2+0.5
k3 = temp[3]- 2*k2
k4 = temp[4]-(3*k3+k2)
k5 = temp[5]-(4*k4+3*k3)
k6 = temp[6]-(5*k5+6*k4+k3)
k7 = temp[7]-(6*k6+10*k5+4*k4)
print(k0,k1, k2,10/6*temp[1])
"""
#temp = taylor(lambda s: ln(t__(s)), 0.5, 10)
#print([temp[i] for i in range(0,10)])
#temp = taylor(lambda s: ln(f1(s)), 0.5000001, 10)
#print([temp[i]*1/(2*0.5000001-1) for i in range(0,10)])

'''
for s0 in [1,0,0.5]:
  print(s0,':------------------------------------------------')
  s0 = s0+h
  print('*****')
  ac0 = ((ln(f1(s0))+c0*(s0*(s0-1)))*ln(s0*(s0-1))**-1)
  print('f0[0]:',ac0)  
  """
  bc0 = (ln(t__(s0))-5)/(s0*(s0-1))
  print('f0[0]:',bc0)
  print('*****')
  s1=s0+h
  ac0 = (ln(t__(s1))-ln(t__(s0)))/(s1*(s1-1)-s0*(s0-1))
  print('f0[0]:',ac0)
  bc0 = (ln(f1(s1))-ln(f1(s0)))/(s1*(s1-1)-s0*(s0-1))
  print('f0[0]:',bc0)
  print('*****')
  s2=s0+2*h
  ac1 = (ln(t__(s2))-ln(t__(s1)))/(s2*(s2-1)-s1*(s1-1))
  ac2 = (ac1-ac0)/s2#/(s2*(s2-1)-s0*(s0-1))
  print('f0[0]:',ac0,ac1,ac2)
  bc1 = (ln(t__(s2))-ln(t__(s1)))/(s2*(s2-1)-s1*(s1-1))
  bc2 = (bc1-bc0)/(s2*(s2-1)-s0*(s0-1))
  print('f0[0]:',bc2)
  print('*****')
  """
  print('----')
'''

"""
print(1/f1x(0.000001),bernoulli(0))
print(1/f1x(1.000001),bernoulli(1)) 
for i in range(2,20):
  print(f1x(i),bernoulli(i)) 
print('----------')
h=0.000000001
for i in range(2,20):
  print(f1x(-2*i+h),bernoulli(i)) 
"""
s= 0.8
t=zetazero(1).imag-0.1
#plot([lambda t: f1(s+j*t),lambda t: f1__(s+j*t),lambda t: 0*f1_(s+j*t)],[14,19])
#plot([lambda t: f0(s+j*t),lambda t: f0__(s+j*t),lambda t: 0*c0*mf(s+j*t)],[0,15])
#plot([lambda s: f1x(s+j*t)],[0,1])