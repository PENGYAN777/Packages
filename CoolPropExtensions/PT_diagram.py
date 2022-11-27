#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:05:08 2022

@author: yan

plot (P,T) diagram. need to adjust axis range and position of text for different fluids.
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

# ------
# Labels
# ------

plt.plot(Ts,ps,'orange',lw = lw, solid_capstyle = 'round')

# Critical lines
plt.axvline(Tc, dashes = [2, 2])
plt.axhline(pc, dashes = [2, 2])

# Labels
plt.text(700, 1e7, 'supercritical',ha= 'center')
plt.text(700, 1e3, 'supercritical_gas', rotation = 90,ha= 'center')
plt.text(450, 1e7, 'supercritical_liquid', rotation = 0, ha = 'right')
plt.text(300, 3e5, 'liquid', rotation = 45)
plt.text(400, 5e3, 'gas', rotation = 45)


plt.ylim(pmin,pmax)
plt.gca().set_yscale('log')
plt.gca().set_xlim(Tmin-50, 800)
plt.ylabel('Pressure [Pa]')
plt.xlabel('Temperature [K]')
plt.title('P-T diagram for siloxane MM')
plt.tight_layout()
fig.savefig("P-T.pdf")
print("PT-diagram.py called")