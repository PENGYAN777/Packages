#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:30:25 2024

plot (P,v) diagram. 

@author: yan
"""
# G = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P,'T', T,fluidname)
# P = CP.CoolProp.PropsSI('P','T', T, 'Q', 1,  fluidname)

import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import pandas as pd
import math

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


# ----------------
# DSL
# ---------------
# PAd = [] # pressure of Double sonic A point
# vAd = []
# PBd = [] # pressure of Double sonic B point
# vBd = []

M1 = 1.0
s = 515
n1 = 10
v1 = np.linspace(vc*1.1, vc*2.5,n1) 
v1 = pd.Series(v1)
P1 = np.zeros(v1.size) 
n2 = 50
M2 = np.zeros((n1,n2))
for j in v1.index:
    P1[j] = CP.CoolProp.PropsSI('P','Smass', s, 'Dmass', 1/v1[j],  fluidname)
    c1 = CP.CoolProp.PropsSI('A','P', P1[j], 'Dmass', 1/v1[j],  fluidname)
    h1 = CP.CoolProp.PropsSI('Hmass','P', P1[j], 'Dmass', 1/v1[j],  fluidname)
    u1 = M1*c1
    ht1 = h1 + 0.5 * u1*u1
    # compute Rayleigh line
    P2 = np.linspace(Pc*0.7, min(Pc*0.99,P1[j]) ,n2) 
    P2 = pd.Series(P2)
    v2r = np.zeros(P2.size) # v2 for Rayleigh line

    for i in P2.index:
        P = P2[i]
        d2s =  CP.CoolProp.PropsSI('Dmass','P',P,'Q',1,fluidname) # for saturation
        v2s = 1/d2s
        v2r[i] = (1-(P-P1[j])*v1[j]/u1/u1)*v1[j]
        if v2r[i]<=v2s:
            # print("two pahse, v2s/vc and vtr/vc: ", v2s/vc, v2r[i]/vc)
            continue
        c2 = CP.CoolProp.PropsSI('A','P', P, 'Dmass', 1/v2r[i],  fluidname)
        u22 =(2* ht1 - CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', 1/v2r[i],  fluidname))
        u2 = math.sqrt(u22)
        M2[j,i] = u2/c2
    # break     
        


# PAd.append(P1[i])
# vAd.append(v1[i])
# PBd.append(P2[j])
# vBd.append(v2r[j])
        
               
# Pd = PAd +  PBd
# vd = vAd +  vBd           
# Pd = np.asarray(Pd) # convert list to array
# vd = np.asarray(vd)
# ----------------
# Plot
# ---------------


# fillcolor = 'g'

# fig = plt.figure(figsize = (6,6))
# ax = fig.add_subplot(111)
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2
# LVS
axes.plot(vsl/vc,psl/Pc,'k',lw = lw, label = "LVS")
axes.plot(vsv/vc,psv/Pc,'k',lw = lw)


# criticla point
axes.plot(vc/vc,Pc/Pc,'ko',lw = lw , label = "Critical point")

# Gamma  curve
axes.plot(vg0/vc,pg0/Pc,'b',lw = lw, label = "$\Gamma=0$")

# isentrope
# axes.plot(v1/vc,P1/Pc,'k--',lw = lw/2, label = "s")

# Rayleigh line
axes.plot(v2r/vc,P2/Pc,'b--',lw = lw, label = "RL")

# DSL
# axes.plot(vd/vc,Pd/Pc,'ro',lw = lw/2, label = "DSL")

axes.set_ylim([0.5, 1.1])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("RSR.pdf") 