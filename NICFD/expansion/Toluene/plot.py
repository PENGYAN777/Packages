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
"""
read csv file
"""

z34= pd.read_csv("z34.csv", ",", skiprows=0)

"""
plot 
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z34.iloc[:,3] , z34.iloc[:,5] , 'k', lw=lwh, label="$Z_t = 0.34$")

axes.set_xlabel('$\\rho/\\rho_c$',fontsize=12)
axes.set_ylabel('Mach',fontsize=12) 
axes.set_title('Mach number vs $\\rho/\\rho_c$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("Mach_rho.pdf")

fig2 = plt.figure( dpi=300)
lwh = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z34.iloc[:,5] , z34.iloc[:,6] , 'k', lw=lwh, label="$Z_t = 0.34$")

axes.set_xlabel('Mach',fontsize=12)
axes.set_ylabel('$\\nu$',fontsize=12) 
axes.set_title('$\\nu$ vs Mach number',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("Mach_nu.pdf")