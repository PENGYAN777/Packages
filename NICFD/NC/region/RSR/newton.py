#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:40:08 2024

Newton method

@author: yan
"""
import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import pandas as pd
import math
import time
import os
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')

start = time.time()
fluidname = "HEOS::MD4M"

"""
0. function
"""
def f(d1,u1,p1,d2,u2,p2):
    h1 = CP.CoolProp.PropsSI('Hmass','P', p1, 'Dmass', d1,  fluidname)
    h2 = CP.CoolProp.PropsSI('Hmass','P', p2, 'Dmass', d2,  fluidname)
    ht1 = h1 + 0.5*u1*u1
    ht2 = h2 + 0.5*u2*u2
    return  np.array( [ [d1*u1 - d2*u2], 
                        [d1*u1*u1 + p1 - d2*u2*u2 - p2],
                        [ht1-ht2] ] )

def fp(d1,u1,p1,d2,u2,p2):
    dh2dd = CP.CoolProp.PropsSI('d(Hmass)/d(Dmass)|P','P', p2, 'Dmass', d2,  fluidname)
    dh2dp = CP.CoolProp.PropsSI('d(Hmass)/d(P)|Dmass','P', p2, 'Dmass', d2,  fluidname)
    return np.array( [ [-u2, -d2, 0], 
                       [-u2*u2, -2*d2*u2, -1],
                       [-dh2dd, -u2, -dh2dp] ]   )

def module(x1,x2,x3):
    return np.sqrt(x1**2 + x2**2 + x3**2)

"""
1. initial guess 
"""
# tol = 1e-3 # torelance
# nmax = 100 # number of iteration
# x = np.zeros( (2,nmax) )
# x0 = (0,0) 
# x[0,0] = x0[0]
# x[1,0] = x0[1]

# n = 1

"""
2. while iteration
"""
# while n<nmax-1:
#     a = f(x[0,n-1],x[1,n-1])
#     b = fp(x[0,n-1],x[1,n-1])
#     dx = - np.matmul(np.linalg.inv(b), a)
#     x[0,n] = x[0, n-1] + dx[0]
#     x[1,n] = x[1, n-1] + dx[1]
#     diff = abs( (module(x[0,n],x[1,n]) - module(x[0,n-1],x[1,n-1]) ) / module(x[0,n-1],x[1,n-1]) )
#     print("diff: ", diff)
#     if abs( (module(x[0,n],x[1,n]) - module(x[0,n-1],x[1,n-1]) ) / module(x[0,n-1],x[1,n-1]) )<tol:
#         print(" n = ", n)
#         print("x1,x2,f(x)", x[0,n], x[1,n], f(x[0,n],x[1,n]) )
#         break
#     else:
#         n = n + 1
        


end = time.time()
print("computational time(s): ", end - start)
