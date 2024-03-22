#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 13:23:33 2024

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import CoolProp as CP
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')

"""
1. read csv file
"""
inlet = pd.read_csv("inlet.csv", ",", skiprows=0)
outlet = pd.read_csv("outlet.csv", ",", skiprows=0)
sym = pd.read_csv("symmetry.csv", ",", skiprows=0)
wall = pd.read_csv("wall.csv", ",", skiprows=0)






"""
2. plot
"""
# fig 1
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(inlet.iloc[:,0] , inlet.iloc[:,1], 'r', lw=lwh, label="inlet")
axes.plot(outlet.iloc[:,0] , outlet.iloc[:,1], 'g', lw=lwh, label="outlet")
axes.plot(sym.iloc[:,0] , sym.iloc[:,1], 'b', lw=lwh, label="symmetry")
axes.plot(wall.iloc[:,0] , wall.iloc[:,1], 'k', lw=lwh, label="wall")

axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# sub_axes = plt.axes([0.26, 0.26, 0.25, 0.25]) 

# # sub_axes.plot(m6new.iloc[:,6] , m6new.iloc[:,-3], 'k', lw=lwh, label="$s$")
# sub_axes.plot(m4new.iloc[:,11] , m4new.iloc[:,-3], 'k', lw=lwh, label="$s$")
# sub_axes.set_ylabel('$s[J/K/mol]$',fontsize=10) 
# sub_axes.set_ylim([754.5, 760])

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$Y[mm]$',fontsize=12) 
axes.set_title('Geometry of the wedge',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("shock_geo.pdf")