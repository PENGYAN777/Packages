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
z9= pd.read_csv("z9.csv", ",", skiprows=0)
z8= pd.read_csv("z8.csv", ",", skiprows=0)
z7= pd.read_csv("z7.csv", ",", skiprows=0)
z6= pd.read_csv("z6.csv", ",", skiprows=0)

"""
plot 
"""

fig0 = plt.figure( dpi=300)
lwh = 2
axes = fig0.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure

axes.plot(z9.iloc[:,2] , z9.iloc[:,6]-1-z9.iloc[:,3]/z9.iloc[99,3] , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(z8.iloc[:,2] , z9.iloc[:,6]-1-z8.iloc[:,3]/z8.iloc[99,3]  , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(z7.iloc[:,2] , z9.iloc[:,6]-1-z7.iloc[:,3]/z7.iloc[99,3]  , 'b', lw=lwh, label="$Z_t = 0.7$")
axes.plot(z6.iloc[:,2] , z9.iloc[:,6]-1-z6.iloc[:,3]/z6.iloc[99,3] , 'k--', lw=lwh, label="$Z_t = 0.6$")

axes.set_xlabel('$P/P_t$',fontsize=12)
axes.set_ylabel('$M-1-\sqrt{{\eta}}$',fontsize=12) 
# axes.set_title('$M-1-\sqrt{{\eta}}$ vs $X$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
# axes.set_ylim([-6, 0])




# fig1 = plt.figure( dpi=300)
# lwh = 2
# axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,3] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,3] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,3] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,3] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('$\\rho/\\rho_t$',fontsize=12) 
# axes.set_title('$\\rho/\\rho_t$ vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig1.savefig("isen_rho.pdf")

# fig2 = plt.figure( dpi=300)
# lwh = 2
# axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,4] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,4] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,4] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,4] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('$T/T_t$',fontsize=12) 
# axes.set_title('$T/T_t$ vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig2.savefig("isen_t.pdf")

# fig3 = plt.figure( dpi=300)
# lwh = 2
# axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,5] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,5] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,5] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,5] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('$c/c_t$',fontsize=12) 
# axes.set_title('$c/c_t$ vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig3.savefig("isen_c.pdf")

# fig4 = plt.figure( dpi=300)
# lwh = 2
# axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,6] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,6] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,6] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,6] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('Mach',fontsize=12) 
# axes.set_title('Mach vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig4.savefig("isen_m.pdf")

# fig5 = plt.figure( dpi=300)
# lwh = 2
# axes = fig5.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,7] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,7] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,7] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,7] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('Z',fontsize=12) 
# axes.set_title('Z vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig5.savefig("isen_z.pdf")

# fig6 = plt.figure( dpi=300)
# lwh = 2
# axes = fig6.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,2] , z9.iloc[:,8] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,2] , z8.iloc[:,8] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,2] , z7.iloc[:,8] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,2] , z6.iloc[:,8] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('$\\Gamma$',fontsize=12) 
# axes.set_title('$\\Gamma$ vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig6.savefig("isen_g.pdf")

# fig7 = plt.figure( dpi=300)
# lwh = 2
# axes = fig7.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# # delete last raw since v = 0 at total contidion
# axes.plot(z9.iloc[0:98,2] , z9.iloc[0:98,9] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[0:98,2] , z8.iloc[0:98,9] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[0:98,2] , z7.iloc[0:98,9] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[0:98,2] , z6.iloc[0:98,9] , 'k--', lw=lwh, label="$Z_t = 0.6$")
# axes.set_xlabel('$P/P_t$',fontsize=12)
# axes.set_ylabel('$A/A^*$',fontsize=12) 
# axes.set_title('$A/A^*$ vs $P/P_t$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig7.savefig("isen_A.pdf")