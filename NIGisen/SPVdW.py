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
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def equations(p,zt,gamma,A):
    """
    input: Zt, gamma, A
    output: M, Z
    """
    M, z = p #output
    t = zt # total compressibility factor
    g = gamma # spefific heat ratio
    
    return ( M**2-(((z-1)/(t-1))**(1-g)-((g-1)*z+1)/((g-1)*t+1) )/((g-1)/2*g*z*z/((g-1)*z+1))  ,
            ((z-1)/(t-1))**((g+1)/2) - A*M*( ((g-1)*z+1)/((g-1)*t+1)+ ((g-1)/2*g*z*z/((g-1)*z+1)) )**((1+g)/(2-2*g))  
            )
    
# vector of Zt
t_vec = np.linspace(1.1, 2.0, 10) # start, end ,number of points
t_vec = pd.Series(t_vec)
gamma = 1.4
A = 2  # A/A^*
# create array
M_vec = np.zeros(t_vec.size) # Mach number
Z_vec = np.zeros(t_vec.size) # compressibility factor
tr_vec = np.zeros(t_vec.size)  # Tt/T
pr_vec = np.zeros(t_vec.size)  # Pt/T

for i in t_vec.index:
    #print("Zt = :", t_vec[i])
    #print("solving equation")
    M_vec[i], Z_vec[i] = fsolve(equations, (1.1, 1.1), args=(t_vec[i],gamma,A))
    tr_vec[i] = ((Z_vec[i]-1)/(t_vec[i]-1) )**(1-gamma)
    pr_vec[i] = tr_vec[i]**(-gamma/(1-gamma)) 
    #print("M,Z:",M_vec[i],Z_vec[i])
    
# plot (M,Z)
fig1 = plt.figure( dpi=300)
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(t_vec, M_vec, 'ko', lw=2)

axes2=axes.twinx()
axes2.plot(t_vec, Z_vec, 'ro', lw=2)

axes.set_xlabel('$Z_t$',fontsize=14)
#axes.set_yscale("log")
axes.set_ylabel('M',fontsize=12) 
axes2.set_ylabel('Z',fontsize=12,color="red") 
axes2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_title('$M(Z_t)$ and $Z(Z_t)$ at fixed position',fontsize=17)
#axes.set_xlim(0,10)
#axes.legend(loc=2) # 2 means left top
fig1.savefig("files/MZ_Zt.pdf")

# plot (T_t/T,P_t/P)
fig2 = plt.figure( dpi=300)
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(t_vec, 1/tr_vec, 'ko', lw=2)

axes2=axes.twinx()
axes2.plot(t_vec, 1/pr_vec, 'ro', lw=2)

axes.set_xlabel('$Z_t$',fontsize=14)
#axes.set_yscale("log")
axes.set_ylabel('$T/T_t$',fontsize=12) 
axes2.set_ylabel('$P/P_t$',fontsize=12,color="red") 
axes2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.set_title('$T/T_t(Z_t)$ and $P/P_T(Z_t)$ at fixed position',fontsize=17)
#axes.set_xlim(0,10)
#axes.legend(loc=2) # 2 means left top
fig2.savefig("files/PT_Zt.pdf")
