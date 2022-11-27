#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:05:08 2022

@author: yan

plot (P,T) diagram. show contour of compressibility factor Z and 
fundamental derivative of gasdynamics Gamma
"""

import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate


# give fluid name
fluidname = "MM"

# update fluid
fluid = CP.AbstractState("HEOS", fluidname)

pc = fluid.keyed_output(CP.iP_critical)
Tc = fluid.keyed_output(CP.iT_critical)
Tmin =  CP.CoolProp.PropsSI("Ttriple",fluidname)
Tmax =  CP.CoolProp.PropsSI("Tmax",fluidname)
pmax = fluid.keyed_output(CP.iP_max)
pmin = fluid.keyed_output(CP.iP_min)
fillcolor = 'g'

fig = plt.figure(figsize = (6,6))
ax = fig.add_subplot(111)
lw = 3


# ----------------
# Saturation curve
# ----------------
Ts = np.linspace(Tmin, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,fluidname)

# ----------------
# Contour of Z and Gamma
# ----------------
n = 100 # number of points
x = np.linspace(Tmin-50, 1000,n)
y = np.linspace(1e5, pmax,n)
X,Y = np.meshgrid(x,y)
Z =  CP.CoolProp.PropsSI('Z','T',X,'P',Y,fluidname)
Gamma =  CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',X,'P',Y,fluidname)
#print(Z.shape)
#cnt = ax.contour(X, Y, Z, levels=[0.5,0.6,0.7,1,1.1,1.2, 1.5], cmap=matplotlib.cm.RdGy,vmin=abs(Z).min(), vmax=abs(Z).max(), extent=[0, 1, 0, 1])
plt.contour(X, Y, Z, [0.5,0.7,0.9,1, 1.1, 1.3, 1.5], cmap='rainbow')
plt.contourf(X, Y, Gamma, [0.5,0.7,0.9,1, 1.1, 1.3, 1.5], cmap='rainbow')
plt.colorbar()
# ------
# Labels
# ------

plt.plot(Ts,ps,'k',lw = lw, solid_capstyle = 'round', label = "LVS")
ax.legend(loc=4) # 2 means left top
# Critical lines
plt.axvline(Tc, dashes = [2, 2])
plt.axhline(pc, dashes = [2, 2])


plt.ylim(1e5,pmax)
plt.gca().set_yscale('log')
plt.gca().set_xlim(Tmin-50, 800)
plt.ylabel('Pressure [Pa]')
plt.xlabel('Temperature [K]')
plt.title('P-T diagram for siloxane MM')
plt.tight_layout()
fig.savefig("Contour.pdf")
print("plotcontour.py called")