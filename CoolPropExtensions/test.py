#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 11:36:47 2022

@author: yan

learn plot in python
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP

# get fluid name
fluid_name = "nitrogen"

# number of points for (P,T) diagram
n = 50

# get P.T range from triple point to max point. unit[Pa] and [K]
P_triple = CP.PropsSI("ptriple",fluid_name)
T_triple = CP.PropsSI("Ttriple",fluid_name)
P_max = CP.PropsSI("pmax",fluid_name)
T_max = CP.PropsSI("Tmax",fluid_name)
P_cri = CP.PropsSI("pcrit",fluid_name)
T_cri = CP.PropsSI("Tcrit",fluid_name)

#P = np.arange(P_cri , P_max, (P_max-P_triple)/1000)  # arguments: start, stop, step
T = np.arange(T_triple , T_max, (T_max-T_triple)/1000)  
P = CP.PropsSI('P', 'T', T, 'Q', 0, fluid_name)

fig = plt.figure( dpi=300)

#axes = fig.add_axes([0.1, 0.1, 0.4, 0.4]) # left, bottom, width, height (range 0 to 1)
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8]) #, frameon=False)
#axes.plot(T, P, 'r', lw=1, ls='dashed', marker='o', markersize=2, label="curve1") # label = r"$y = \alpha^2$" for latex format
axes.plot(T, P, 'k', lw=2,  label="curve1")
axes.plot(T_triple, P_triple, 'go', lw=2,  label="triple point")
axes.plot(T_cri, P_cri, 'ro', lw=2,  label="critial point")
axes.plot(T_max, P_max, 'bo', lw=2,  label="max point")
axes.set_xlabel('T[K]',fontsize=18)
#axes.set_yscale("log")
axes.set_ylabel('P[Pa]',fontsize=18) 
axes.set_title('P-T diagram',fontsize=18)
#axes.set_xlim(0,10)
axes.legend(loc=2) # 2 means left top
#axes.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
#axes.set_xticks([T_triple,T_max])
#axes.set_xticklabels([r'$T_{triple}$', r'$T_{max}$'], fontsize=18)

# scientif noation
from matplotlib import ticker
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 
axes.yaxis.set_major_formatter(formatter) 


plt.show()
fig.savefig("P-T.pdf")