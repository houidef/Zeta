from mpmath import *
import matplotlib.pyplot as plt
import numpy as np

j=complex(0,1)
def zeta2(s):
   t=0
   signe = 1
   for i in range(1,100000):
     x=i**(-s)
     t+=signe*abs(x)
     #print(x.real,x.imag,abs(x),arg(x))
     signe = -signe
   return t#1/(1-2**(1-s))*t

m=5
#print(1/(m+1)**s)
#print(inv1(s,m))
def theta(t):
   return arg(gamma(complex(0.25,0.5*t)))-0.5*t*log(pi)

def theta0(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   #z3 = 1/pi**(z.conjugate()/2)*gamma(z.conjugate()/2)
   #z4 = 1/pi**((1-z).conjugate()/2)*gamma((1-z).conjugate()/2)
   return j*log(z1).imag

def theta_(z):
   z1 = 1/pi**(z.conjugate()/2)*gamma(z.conjugate()/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   return -j*log(z1*z2).imag

def Zt(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   return (z1/z2)#
 
def Zt2(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   return -log(abs(z1/z2))  
def theta2(z):
   z1 = 1/pi**(z/2)*gamma(z/2)
   z2 = 1/pi**((1-z)/2)*gamma((1-z)/2)
   #z3 = 1/pi**(z.conjugate()/2)*gamma(z.conjugate()/2)
   #z4 = 1/pi**((1-z).conjugate()/2)*gamma((1-z).conjugate()/2)
   return j*log(z2).imag

def eps(z):
   return 1/pi**(z/2)*gamma(z/2)*zeta(z)
   
def epss(z):
   x=1/pi**(z/2)*gamma(z/2)*zeta(z)
   return x/abs(x)
   
def prod0(z):
    return exp(theta0(z))*exp(theta0(1-z))*zeta(z)*zeta(1-z)  

#here1 :	
def prod(z):
    return exp(theta0(z.conjugate()))*exp(theta0(1-z))*zeta(z.conjugate())*zeta(1-z)  

def prodw(z):
    return (exp(theta0(z.conjugate()))*exp(theta0(1-z)))
	
def prod1(z):
    return zeta(z)*zeta(1-z) 
def prod2(z):
    return zeta(z)*zeta(z.conjugate())
#here2 :
def prod3(z):
    return zeta(z.conjugate())*zeta(1-z) 	
  	
 # setting the axes at the centre

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


x = np.linspace(5,50,300)
y = np.linspace(0.1,200,100)
Yr=[]
Yb=[]
Yy=[]
Yg=[]


'''
for t in x:
    #k0 =  prod(complex(s,t)) #(exp(j*theta(t))*zeta(complex(0.5,t)))**2
    k0 =  exp((theta0(complex(s,t))))#
    Yr.append(k0.real)
    Yi.append(k0.imag)
    k = ((abs(zeta(complex(0.5,t))))**2)#prod(complex(0.5,t)) #(zeta(complex(0.5,t))**2).conjugate()
    k = exp(j*(arg(zeta(complex(s,t)))))
    Yrr.append(k.real)
    Yii.append(k.imag)
'''
def fun(z):
   return (zeta(z)*zeta(1-z))

def fun1(z):
  return (zeta(z,1,1)*zeta(1-z)+zeta(z)*zeta(1-z,1,1))
  
def derive(z,n):
   h = 0.0001
   if(n!=0):
      return (derive(z,n-1)-derive(z+j*h,n-1))/h
   else:
      return fun(z)

def derives(z,n):
   h = 0.0001
   if(n!=0):
      return (derives(z,n-1)-derives(z+h,n-1))/h
   else:
      return fun(z)
	  
#0.05-0.02	
s = 2
def fac(s,n):
   if(n==0) : return 1
   return 1/(s-0.5)**2*(fac(s,n-1)+1)
for t in x:
    k =  abs(zeta(s+j*t))**2*(Zt((s+j*t)))/abs(Zt((s+j*t)))-1#abs(zeta(s+j*t,derivative=1))/abs(Zt((s+j*t)))#(derive(complex(s,t),0))/(Zt(complex(s,t))) #(derive(complex(s,t),0))#prodw(complex(s,t)) #(exp(j*theta(t))*zeta(complex(0.5,t)))**2
    
	#k = derive(complex(s,t),0)
    '''
    d = derive(complex(s,t),1)
    if(abs(d.real)<0.05):
        print('point1:',t,k.real)
    if(abs(d.imag)<0.05):
        print('point2:',t,k.imag)
    '''
    #print(k0)
	#15/22.85*log(k)-1.64*(s)+5.87
    #k = 1/(s+j*t)**(s)
    #k = (k)/(j*t)**(0.03)#2**10*(k0)/(t)**(s-0.5)#(0.030+0.0015*j*t) #0.514  0.030  prod(complex(0.5,t)) #(zeta(complex(0.5,t))**2).conjugate()
    #print(t,k.real-k0.real,cos(-pi/2*0.3)*1/t**0.3,cos(-pi/2*0.3)**(1/0.3))

    #k = 1/t**(-2*s+1)*k/(t+5)**(2*s-1) 
    Yr.append(k.real)
    Yb.append(k.imag)
    k = ((zeta(s+j*t))**2-1)
    #abs(derive(complex(s,t),0)) #derive(complex(s,t),1)
    h=0.001
    #k0 = 0.27*derive(complex(s,t),1) #equation de houidef
    #k = exp(-j*theta(t)*2)#derive(0omplex(s,t),1) #(derive(complex(s1,t),0)-derive(complex(s1+h,t),0))/h
    #k = 16.55*(derive(complex(0.5,t),0)-derive(complex(0.5-h,t),0))/h
    #k=derive(complex(s,t),2)#1/gamma(complex(s,t)/10)#exp(theta_(complex(s,t)))
    Yy.append(k.real)
    Yg.append(k.imag) #green

"""
t=zetazero(1).imag
def findzeroo():
   m=0
   for r in y:
      m1 = m
      m = derive(complex(r,t),0)
      if(m.imag*m1.imag<0):
         break
   print(r,m.imag)
findzeroo()
"""
"""
t=zetazero(1).imag
for s in y:
    k = abs(derive(complex(s,t),0))
    Yr.append(abs(log(k.real)))
    Yb.append(abs(Zt(complex(s,t))))
    #k=derive(complex(s,t),2)#1/gamma(complex(s,t)/10)#exp(theta_(complex(s,t)))
    #Yy.append(-k.real)
    #Yg.append(k.imag) #green
	
"""


#plot(lambda x: -log(abs(1/pi**(complex(1,x)*2/pi)*gamma(complex(1,x)*2/pi)))-0.998*x+0.6-1/(x**0.6+1),[0.1,100])

   
#plot([lambda x: derive(complex(0.6,x),0).imag,lambda x: derive(complex(0.6,x),1).real],[0.1,100])
plt.plot(x,Yr, 'r')
plt.plot(x,Yb, 'b')
plt.plot(x,Yg, 'g')
plt.plot(x,Yy, 'y')

z = 0.5 + 10j
s = z.real
t = z.imag
k=zeta2(z)
print('zeta2',1/(1-2**(1-s.real))*k)
print('zeta2',1/(1-2**(1-s))*k)
print('zeta : k =',k)
print('zeta',zeta(z)/abs(zeta(z)))
print('zeta 0.5 : ',zeta(s))

print('theta(t) = ',exp(-j*theta(t)),theta(t))
print('Z(t) =  ',exp(j*theta(t))*zeta(z))
print('arg : ',epss(z),epss(z.conjugate()))
print('theta0 : ',exp(-theta0(z)),theta0(z))


#print('zero',zetazero(1))

plt.show()