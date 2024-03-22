#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 19:01:08 2022

@author: yan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

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


fluidname = "HEOS::MD4M"
# fluidname = "HEOS::D6"
Pc = CP.CoolProp.PropsSI('Pcrit',fluidname)
Tc = CP.CoolProp.PropsSI('Tcrit',fluidname)
dc = CP.CoolProp.PropsSI('rhocrit',fluidname)



mesh3 = pd.read_csv("m3new.csv", ",", skiprows=0)
mesh4 = pd.read_csv("m4new.csv", ",", skiprows=0)
mesh5 = pd.read_csv("m5new.csv", ",", skiprows=0)
mesh6 = pd.read_csv("m6new.csv", ",", skiprows=0)



# fig 1
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(mesh3.iloc[:,-6] , mesh3.iloc[:,7]/Pc, 'k', lw=lwh, label="22k")
axes.plot(mesh4.iloc[:,-6] , mesh4.iloc[:,7]/Pc, 'r', lw=lwh, label="40k")
axes.plot(mesh5.iloc[:,-6] , mesh5.iloc[:,7]/Pc, 'g', lw=lwh, label="74k")
axes.plot(mesh6.iloc[:,-6] , mesh6.iloc[:,7]/Pc, 'b', lw=lwh, label="140k")
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# sub_axes = plt.axes([0.26, 0.26, 0.25, 0.25]) 

# # sub_axes.plot(m6new.iloc[:,6] , m6new.iloc[:,-3], 'k', lw=lwh, label="$s$")
# sub_axes.plot(m4new.iloc[:,11] , m4new.iloc[:,-3], 'k', lw=lwh, label="$s$")
# sub_axes.set_ylabel('$s[J/K/mol]$',fontsize=10) 
# sub_axes.set_ylim([754.5, 760])

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$P/P_c$',fontsize=12) 
axes.set_title('$P/P_c$ along y$=0.3 mm$',fontsize=14)
axes.legend(loc=1) # 

fig1.savefig("shock_gv_p.pdf")

# fig 2
fig2 = plt.figure( dpi=300)
lwh = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(mesh3.iloc[:,-6] , mesh3.iloc[:,3], 'k', lw=lwh, label="22k")
axes.plot(mesh4.iloc[:,-6] , mesh4.iloc[:,3], 'r', lw=lwh, label="40k")
axes.plot(mesh5.iloc[:,-6] , mesh5.iloc[:,3], 'g', lw=lwh, label="74k")
axes.plot(mesh6.iloc[:,-6] , mesh6.iloc[:,3], 'b', lw=lwh, label="140k")
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('Mach',fontsize=12) 
axes.set_title('Mach number along y$=0.3 mm$',fontsize=14)
axes.legend(loc=1) # 

fig2.savefig("shock_gv_m.pdf")

# fig 3
fig3 = plt.figure( dpi=300)
lwh = 2
axes = fig3.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(mesh3.iloc[:,-6] , mesh3.iloc[:,1]/dc, 'k', lw=lwh, label="22k")
axes.plot(mesh4.iloc[:,-6] , mesh4.iloc[:,1]/dc, 'r', lw=lwh, label="40k")
axes.plot(mesh5.iloc[:,-6] , mesh5.iloc[:,1]/dc, 'g', lw=lwh, label="74k")
axes.plot(mesh6.iloc[:,-6] , mesh6.iloc[:,1]/dc, 'b', lw=lwh, label="140k")
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$\\rho$/$\\rho_c$',fontsize=12) 
axes.set_title('$\\rho$ / $\\rho_c$ along y$=0.3 mm$',fontsize=14)
axes.legend(loc=1) # 

fig3.savefig("shock_gv_rho.pdf")

# fig 4
fig4 = plt.figure( dpi=300)
lwh = 2
axes = fig4.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(mesh3.iloc[:,-6] , mesh3.iloc[:,-3], 'k', lw=lwh, label="22k")
axes.plot(mesh4.iloc[:,-6] , mesh4.iloc[:,-3], 'r', lw=lwh, label="40k")
axes.plot(mesh5.iloc[:,-6] , mesh5.iloc[:,-3], 'g', lw=lwh, label="74k")
axes.plot(mesh6.iloc[:,-6] , mesh6.iloc[:,-3], 'b', lw=lwh, label="140k")
axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# sub_axes = plt.axes([0.55, 0.25, 0.25, 0.25]) 
# x = []
# plot the zoomed portion
# for i in mesh4.iloc[:,-3].index:
#     if (mesh4.iloc[i,-3])<0:
#         x.append(i)
# sub_axes.plot(mesh4.iloc[x[0]:x[-1],-6], mesh4.iloc[x[0]:x[-1],-3] , 'b') 
# sub_axes.set_ylim([-0.007, -0.002])
# sub_axes.set_xlim([0, 0.4])
# sub_axes.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

axes.set_xlabel('$X[mm]$',fontsize=12)
axes.set_ylabel('$\\Gamma$',fontsize=12) 
axes.set_title('$\\Gamma$ along y$=0.3 mm$',fontsize=14)
axes.legend(loc=0) # 

fig4.savefig("shock_gv_Gamma.pdf")

# # fig 5
# fig5 = plt.figure( dpi=300)
# lwh = 2
# axes = fig5.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
# # axes.plot(mesh0.iloc[:,-6] , mesh0.iloc[:,-1], 'k', lw=lwh, label="5k")
# # axes.plot(mesh2.iloc[:,-6] , mesh2.iloc[:,-1], 'r', lw=lwh, label="12k")
# # axes.plot(mesh3.iloc[:,-6] , mesh3.iloc[:,-1], 'g', lw=lwh, label="23k")
# axes.plot(mesh4.iloc[:,-6] , mesh4.iloc[:,-1], 'b', lw=lwh, label="41k")


# axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# axes.set_xlabel('$X[mm]$',fontsize=12)
# axes.set_ylabel('s',fontsize=12) 
# axes.set_title('entropy along y$=0.25 mm$',fontsize=14)
# axes.legend(loc=1) # 

# fig5.savefig("shock_gv_s.pdf")