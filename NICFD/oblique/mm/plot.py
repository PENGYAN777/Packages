#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:59:53 2023

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import math
"""
read csv file
"""

z9= pd.read_csv("z9.csv", ",", skiprows=0)
z8= pd.read_csv("z8.csv", ",", skiprows=0)
z7= pd.read_csv("z7.csv", ",", skiprows=0)


"""
plot 
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z9.iloc[:,2]*180/math.pi , z9.iloc[:,3]*180/math.pi  , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(z8.iloc[:,2]*180/math.pi , z8.iloc[:,3]*180/math.pi  , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(z7.iloc[:,2]*180/math.pi , z7.iloc[:,3]*180/math.pi  , 'b', lw=lwh, label="$Z_t = 0.7$")


axes.set_xlabel('$\\theta(degree)$',fontsize=12)
axes.set_ylabel('$\\beta$(degree)',fontsize=12) 
axes.set_title('$\\theta$ vs $\\beta$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("theta_beta.pdf")

