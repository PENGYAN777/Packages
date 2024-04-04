#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:40:08 2024

Bisection method

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


"""
0. function
"""
def f(x):
    return  x**4-x-2

"""
1. initial guess a and b
"""
a = 1
b = 2
print("f(a)*f(b)<0: ", f(a)*f(b))

tol = 1e-4 # torelance
nmax = 100 # number of iteration

n = 1

"""
2. while iteration
"""
while n<nmax:
    c = (a+b)/2
    if f(c) == 0 or (b-a)/2<tol:
        print("solution found")
        print("c = ", c, "f(c) = ", f(c))
        break
    else:
        n = n + 1
        if f(a)*f(c)>=0:
            a = c
        else:
            b = c
print(" n = ", n)





end = time.time()
print("computational time(s): ", end - start)
