#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 18:04:05 2023

compute downstream Mach number  for nonideal flow

@author: yan
"""
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

"""
0. get data
"""
z9= pd.read_csv("z9.csv", ",", skiprows=0)
z8= pd.read_csv("z8.csv", ",", skiprows=0)
z7= pd.read_csv("z7.csv", ",", skiprows=0)

"""
1. input theta and upstream Mach number, compute downstream data
"""
theta = np.zeros(50) # Gamma
z9_m2 = np.zeros(50) 
z8_m2 = np.zeros(50) 
z7_m2 = np.zeros(50) 

for i in range(50):
    theta[i] = i*math.pi/180 # rad
    M1 = 1.1 # upstream Mach number
    # for z9
    nu1 = z9.iloc[:,6][np.argmin(abs(z9.iloc[:,5]-M1))] 
    nu2 = nu1 + theta[i]
    M2 = z9.iloc[:,5][np.argmin(abs(z9.iloc[:,6]-nu2))] 
    z9_m2[i] = M2
    # for z8
    nu1 = z8.iloc[:,6][np.argmin(abs(z8.iloc[:,5]-M1))] 
    nu2 = nu1 + theta[i]
    M2 = z8.iloc[:,5][np.argmin(abs(z8.iloc[:,6]-nu2))] 
    z8_m2[i] = M2
    # for z7
    nu1 = z7.iloc[:,6][np.argmin(abs(z7.iloc[:,5]-M1))] 
    nu2 = nu1 + theta[i]
    M2 = z7.iloc[:,5][np.argmin(abs(z7.iloc[:,6]-nu2))] 
    z7_m2[i] = M2
    


"""
2. plot
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(theta , z9_m2 , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(theta , z8_m2 , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(theta , z7_m2 , 'b', lw=lwh, label="$Z_t = 0.7$")


axes.set_xlabel('$\\theta$(rad)',fontsize=12)
axes.set_ylabel('M2',fontsize=12) 
axes.set_title('M2 vs $\\theta$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("M2_theta.pdf")
