from mpmath import *
import matplotlib.pyplot as plt
import numpy as np
m0=100
for i in range(0,m0-1):
    h=zetazero(i+1).imag
    print(digamma(i+1))