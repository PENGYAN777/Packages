#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 22:33:59 2024

plot RH, Rayleigh line

@author: yan
"""

import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import pandas as pd

# give fluid name
fluidname = "HEOS::MD4M"
Pc = CP.CoolProp.PropsSI('Pcrit',fluidname)
Tc = CP.CoolProp.PropsSI('Tcrit',fluidname)
dc = CP.CoolProp.PropsSI('rhocrit',fluidname)
vc = 1/dc
w = CP.CoolProp.PropsSI('ACENTRIC',fluidname)
print("fluid name:", fluidname)
R = CP.CoolProp.PropsSI('GAS_CONSTANT',fluidname)
print("Universal gas constant:", R)
MW = CP.CoolProp.PropsSI('M',fluidname)
print("molar mass:", MW)
Rs = R/MW
print("specific gas constant:", Rs)

# ----------------
# Saturation curve
# ---------------
# liquid phase
psl = np.linspace(Pc*0.5, Pc*1.0, 500)
dsl = CP.CoolProp.PropsSI('Dmass','P',psl,'Q',0,fluidname)
vsl = 1/dsl
# vapor phase
psv = np.linspace(Pc*0.5, Pc*1.0, 500)
dsv = CP.CoolProp.PropsSI('Dmass','P',psv,'Q',1,fluidname)
vsv = 1/dsv

# ----------------
# Gamma=0 curve
# ---------------
# consider single capor phase only
pg = np.linspace(Pc*0.7, Pc*1.0,100) #
dg = CP.CoolProp.PropsSI('Dmass','P',pg,'Q',1,fluidname) 
vg = 1/dg
vg = pd.Series(vg)
pg = pd.Series(pg)
pg0 = [] # p w.r.t Gamma = 0
vg0 = []
for i in vg.index:
    pgl = np.linspace(pg[i]*1.005, Pc*1.05,100) # list all possible p
    pgl = pd.Series(pgl)
    G = np.zeros(pgl.size) 
    for j in pgl.index:
        G[j] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',pgl[j],'Dmass', 1/vg[i],fluidname)
        if abs(G[j])<0.001:
            # print("G: ", G[j])
            pg0.append(pgl[j])
            vg0.append(vg[i])
            
pg0 = np.asarray(pg0) # convert list to array
vg0 = np.asarray(vg0)

# ---------------
# initial state1, P,v, M
P1 = Pc*0.95
d1 = CP.CoolProp.PropsSI('Dmass','P',P1,'Q',1,fluidname)/1.07
v1 = 1/d1
# v1 = 1/d1
c1 = CP.CoolProp.PropsSI('A','P', P1, 'Dmass', d1,  fluidname)
h1 =  CP.CoolProp.PropsSI('Hmass','P', P1, 'Dmass', d1,  fluidname)

# ----------------
# RH line, use Eq. (2.2)
n1 = 30
P2 = np.linspace(Pc*0.8, Pc*1.02 ,n1) # post-shock Mach
P2 = pd.Series(P2)
v2 = np.zeros(P2.size) 
for j in P2.index:
    P = P2[j]
    # initial guess pf v2
    n2 = 200
    v = np.linspace(vc*1.1,vc*3 ,n2) 
    v = pd.Series(v)
    diff = np.zeros(v.size) 
    for i in v.index:
        diff[i] = CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', 1/v[i], fluidname) - h1 - 0.5*(P-P1)*(v1+v[i])
    print("min index:", np.argmin(abs(diff)), "min diff: ", diff[np.argmin(abs(diff))])
    v2[j] = v[np.argmin(abs(diff))]
# ----------------
# Rayleigh line, use Eq. (2.4)
M1 = 1.1
u1 = M1*c1
v2r = np.zeros(P2.size) # v2 for Rayleigh line
for i in P2.index:
      P = P2[i]
      v2r[i] = (1-(P-P1)*v1/u1/u1)*v1



# ----------------
# Plot
# ---------------
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2        

# LVS
axes.plot(vsl/vc,psl/Pc,'k',lw = lw, label = "LVS")
axes.plot(vsv/vc,psv/Pc,'k',lw = lw)

# Gamma = 0 curve
axes.plot(vg0/vc,pg0/Pc,'r',lw = lw, label = "$\Gamma=0$")

# # initial state
axes.plot(1/d1/vc,P1/Pc,'bo',lw = lw/2, label = "IS")
# # axes.plot(1/d3/vc,P3/Pc,'bo',lw = lw/2)

# # RHL
axes.plot(v2/vc,P2/Pc,'b',lw = lw, label = "RHL")

# Rayleigh line
axes.plot(v2r/vc,P2/Pc,'b--',lw = lw, label = "RL")

axes.set_ylim([0.7, 1.1])
axes.set_xlim([1.0, 3.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('RHL and Rayleigh line')
axes.legend(loc=0) # 2 means left top
fig1.savefig("RHL_NC.pdf") 