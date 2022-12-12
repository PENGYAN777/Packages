#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

solve non-linear system of equations: M=f(Z,Zt) and Z=g(Zt) for SPVdW
remeber that Z > 1
"""

from scipy.optimize import fsolve
import math
import numpy as np

def equations(p,t):
    
    M, z = p #output
    #t = 1.2 # Zt
    A = 2# A/A^*
    g=1.4# gamma
    
    return ( M**2-(((z-1)/(t-1))**(1-g)-((g-1)*z+1)/((g-1)*t+1) )/((g-1)/2*g*z*z/((g-1)*z+1))  ,
            ((z-1)/(t-1))**((g+1)/2) - A*M*( ((g-1)*z+1)/((g-1)*t+1)+ ((g-1)/2*g*z*z/((g-1)*z+1)) )**((1+g)/(2-2*g))  
            )
    

t_vec = np.linspace(1.1, 2.0, 10) # start, end ,number of points

for t in t_vec:
    
    print("Zt = :", t)
    print("solving equation")
    M, Z = fsolve(equations, (1.1, 1.1), args=(t,))
    print("M,Z:",M,Z)


