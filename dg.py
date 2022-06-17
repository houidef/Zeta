from mpmath import *
mp.dps = 15; mp.pretty = True


i=2
k=1
change = +1
while (abs(psi(1,i) - 0.0624517456841562)> 0.00000000000001) :
   print(i,psi(1,i))
   oldchange = change
   elif (psi(1,i) < 0.0624517456841562) :  
     change = -1 
     i = i - k
   else: 
     change = +1;
     i = i + k 
   if(change != oldchange) : k=k/10 