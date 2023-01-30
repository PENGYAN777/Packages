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
md = pd.read_csv("paper_mrho.csv", ",", skiprows=0)
num = pd.read_csv("paper_num.csv", ",", skiprows=0)
"""
plot 
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z34.iloc[:,3] , z34.iloc[:,5] , 'k', lw=lwh, label="$Z_t = 0.34$")
axes.plot(md.iloc[:,0] , md.iloc[:,1] , 'ko', lw=lwh/2, label="Cramer et.al 1992")


axes.set_xlabel('$\\rho/\\rho_c$',fontsize=12)
axes.set_ylabel('Mach',fontsize=12) 
axes.set_title('Mach number vs $\\rho/\\rho_c$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("vv_mrho.pdf")

fig2 = plt.figure( dpi=300)
lwh = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z34.iloc[:,5] , z34.iloc[:,6] , 'k', lw=lwh, label="$Z_t = 0.34$")
axes.plot(num.iloc[:,0] , num.iloc[:,1] , 'ko', lw=lwh/2, label="Cramer et.al 1992")

axes.set_xlabel('Mach',fontsize=12)
axes.set_ylabel('$\\nu$',fontsize=12) 
plt.xlim([1, 1.75])
axes.set_title('$\\nu$ vs Mach number',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("vv_num.pdf")