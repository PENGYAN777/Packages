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
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z9.iloc[:,-1] , z9.iloc[:,6] , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(z8.iloc[:,-1] , z8.iloc[:,6] , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(z7.iloc[:,-1] , z7.iloc[:,6] , 'b', lw=lwh, label="$Z_t = 0.7$")
axes.plot(z6.iloc[:,-1] , z6.iloc[:,6] , 'k--', lw=lwh, label="$Z_t = 0.6$")

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$M$',fontsize=12) 
axes.set_title('Mach number vs $X$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig1.savefig("M.pdf")

# fig2 = plt.figure( dpi=300)
# lwh = 2
# axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,-1] , z9.iloc[:,2] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,-1] , z8.iloc[:,2] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,-1] , z7.iloc[:,2] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,-1] , z6.iloc[:,2] , 'k--', lw=lwh, label="$Z_t = 0.6$")

# axes.set_xlabel('$X[mm]$',fontsize=12)
# axes.set_ylabel('$P/P_t$',fontsize=12) 
# axes.set_title('$P/P_t$ vs $X$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig2.savefig("P.pdf")

# fig3 = plt.figure( dpi=300)
# lwh = 2
# axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,-1] , z9.iloc[:,3] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,-1] , z8.iloc[:,3] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,-1] , z7.iloc[:,3] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,-1] , z6.iloc[:,3] , 'k--', lw=lwh, label="$Z_t = 0.6$")

# axes.set_xlabel('$X[mm]$',fontsize=12)
# axes.set_ylabel('$\\rho/\\rho_t$',fontsize=12) 
# axes.set_title('$\\rho/\\rho_t$ vs $X$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig3.savefig("rho.pdf")

# fig4 = plt.figure( dpi=300)
# lwh = 2
# axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,-1] , z9.iloc[:,4] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,-1] , z8.iloc[:,4] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,-1] , z7.iloc[:,4] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,-1] , z6.iloc[:,4] , 'k--', lw=lwh, label="$Z_t = 0.6$")

# axes.set_xlabel('$X[mm]$',fontsize=12)
# axes.set_ylabel('$T/T_t$',fontsize=12) 
# axes.set_title('$T/T_t$ vs $X$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig4.savefig("T.pdf")

fig5 = plt.figure( dpi=300)
lwh = 2
axes = fig5.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z9.iloc[:,-1] , z9.iloc[:,7] , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(z8.iloc[:,-1] , z8.iloc[:,7] , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(z7.iloc[:,-1] , z7.iloc[:,7] , 'b', lw=lwh, label="$Z_t = 0.7$")
axes.plot(z6.iloc[:,-1] , z6.iloc[:,7] , 'k--', lw=lwh, label="$Z_t = 0.6$")

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$Z$',fontsize=12) 
axes.set_title('$Z$ vs $X$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
fig5.savefig("Z.pdf")

# fig6 = plt.figure( dpi=300)
# lwh = 2
# axes = fig6.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# axes.plot(z9.iloc[:,-1] , z9.iloc[:,8] , 'k', lw=lwh, label="$Z_t = 0.9$")
# axes.plot(z8.iloc[:,-1] , z8.iloc[:,8] , 'r', lw=lwh, label="$Z_t = 0.8$")
# axes.plot(z7.iloc[:,-1] , z7.iloc[:,8] , 'b', lw=lwh, label="$Z_t = 0.7$")
# axes.plot(z6.iloc[:,-1] , z6.iloc[:,8] , 'k--', lw=lwh, label="$Z_t = 0.6$")

# axes.set_xlabel('$X[mm]$',fontsize=12)
# axes.set_ylabel('$\\Gamma$',fontsize=12) 
# axes.set_title('$\\Gamma$ vs $X$',fontsize=14)
# axes.legend(loc=0 , prop={'size': 10}) # 
# fig6.savefig("G.pdf")

fig7 = plt.figure( dpi=300)
lwh = 2
axes = fig7.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(z9.iloc[:,-1] , z9.iloc[:,6]-1-z9.iloc[:,3]/z9.iloc[247,3] , 'k', lw=lwh, label="$Z_t = 0.9$")
axes.plot(z8.iloc[:,-1] , z9.iloc[:,6]-1-z8.iloc[:,3]/z8.iloc[247,3]  , 'r', lw=lwh, label="$Z_t = 0.8$")
axes.plot(z7.iloc[:,-1] , z9.iloc[:,6]-1-z7.iloc[:,3]/z7.iloc[247,3]  , 'b', lw=lwh, label="$Z_t = 0.7$")
axes.plot(z6.iloc[:,-1] , z9.iloc[:,6]-1-z6.iloc[:,3]/z6.iloc[247,3] , 'k--', lw=lwh, label="$Z_t = 0.6$")

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$M-1-\sqrt{{\eta}}$',fontsize=12) 
axes.set_title('$M-1-\sqrt{{\eta}}$ vs $X$',fontsize=14)
axes.legend(loc=0 , prop={'size': 10}) # 
axes.set_ylim([-6, 0])
# fig7.savefig("G.pdf")
