from mpmath import *

#Re(s)>1 
def zeta0(s):
   t=0
   h=0.001
   for i in range(int(1/h),10000):
     t += int(h*i)*(h*i)**(-s-1)
   return s*t*h
#Re(s)>0   
def zeta2(s):
   t=0
   signe = 1
   for i in range(1,100000):
     x=i**(-s)
     t+=signe*abs(x)
     signe = -signe
   return 1/(1-2**(1-s))*t

#Re(s)>0
def zeta3(s):
   t=0
   h=0.0001
   for i in range(1,10000):
     t += (-ln(i)-ln(h))**(s-1)/(1+i*h)
   return 1/((1-2**(1-s))*gamma(s))*t*h

def zeta33(s):
   t=0
   h=0.0001
   for i in range(1,100000):
     t += i**(s-1)/(1+exp(i*h))
     #print(t)
   return 1/((1-2**(1-s))*gamma(s))*t*h**s
   
#Re(s)>1
def zeta4(s):
   t=0
   h=0.0001
   for i in range(1,10000):
     t += (-ln(i*h))**(s-1)/(1+i*h)
   return 1/((1-2**(1-s))*gamma(s))*t*h
   
def inv(s,m):
   t=0
   h=0.001
   for i in range(1,1000):
     t += (i*h)**m*abs(ln(i*h))**(s-1)
   return 1/gamma(s)*t*h

def inv1(s,m):
   t=0
   h=0.0001
   for i in range(1,10000):
     t += i**m*abs(ln(i)+ln(h))**(s-1)
   return 1/gamma(s)*t*h**(m+1)
borwein_cache = {}

def borwein_coefficients(n):
    if n in borwein_cache:
        return borwein_cache[n]
    ds = [0] * (n+1)
    d = 1
    s = ds[0] = 1
    for i in range(1, n+1):
        d = d * 4 * (n+i-1) * (n-i+1)
        d //= ((2*i) * ((2*i)-1))
        s += d
        ds[i] = s
    borwein_cache[n] = ds
    return ds

def zeta1(s):#, prec, rnd=round_fast, alt=0):
    wp = 100
    n = int(wp/2.54 + 5)
    d = borwein_coefficients(n)
    t = 0
    for k in range(0,n):
        eman = (k+1)**-s
        w = (d[k] - d[n]) * eman
        if k & 1:
            t -= w
        else:
            t += w
    t = t / (-d[n])
    q = 1-2**(1-s)
    return t/q
s= 2
m=5
#print(1/(m+1)**s)
#print(inv1(s,m))
def theta(t):
   return arg(gamma(complex(0.25,0.5*t)))-0.5*t*log(pi)

def theta0(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   z3 = 1/pi**(z.conjugate()/2)*gamma(z.conjugate()/2)
   z4 = 1/pi**((1-z).conjugate()/2)*gamma((1-z).conjugate()/2)
   print(z1)
   print(z2)
   print(z3)
   print(z4)
   return log(z1)


def xi(s):
   return pi**(-s/2)*gamma(s/2)*zeta(z)
   
def xis(s):
   x=pi**(-s/2)*gamma(s/2)*zeta(s)
   return x/abs(x)
s = 2.5#0.5 + 10*j
#k=zeta2(s)
'''
print('zeta2',1/(1-2**(1-s))*k)
print('zeta : k =',k)
print('zeta',zeta(z)/abs(zeta(z)))
print('zeta 0.5 : ',zeta(s))
print('theta(t) = ',exp(-j*theta(t)),theta(t))
print('Z(t) =  ',exp(j*theta(t))*zeta(z))
print('arg : ',epss(z),epss(z.conjugate()))
print('theta0 : ',exp(theta0(z).real),exp(-j*theta0(z).imag),theta0(z))
'''
print('normal',zeta(s))
print('zeta0 :',zeta0(s))
print('zeta1 :',zeta1(s))
print('zeta2 :',zeta2(s))
print('zeta3 :',zeta3(s))
print('zeta3*:',zeta33(s))
print('zeta4 :',zeta4(s))
