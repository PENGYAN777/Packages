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
# isentropes
# ---------------
# s1
v1 = vc*2.5
p1 = Pc*1.05
s1 = CP.CoolProp.PropsSI('Smass','P', p1, 'Dmass', 1/v1,  fluidname)

v1 = np.linspace(vc*0.5, vc*4.5,500) # P<Pc
v1 = pd.Series(v1)
p1 = np.zeros(v1.size) 
for i in v1.index:
    p1[i] = CP.CoolProp.PropsSI('P','Smass', s1-190, 'Dmass', 1/v1[i],  fluidname)
    

# s2
v2 = np.linspace(vc*0.5, vc*4.5,500) # P<Pc
v2 = pd.Series(v2)
p2 = np.zeros(v2.size) 
for i in v1.index:
    p2[i] = CP.CoolProp.PropsSI('P','Smass', s1-100, 'Dmass', 1/v2[i],  fluidname)
    
# s3
v3 = np.linspace(vc*0.5, vc*4.5,500) # P<Pc
v3 = pd.Series(v3)
p3 = np.zeros(v3.size) 
for i in v3.index:
    p3[i] = CP.CoolProp.PropsSI('P','Smass', s1-175, 'Dmass', 1/v3[i],  fluidname)



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

# isentrope

axes.plot(v1/vc,p1/Pc,'k--',lw = lw/2)
axes.plot(v2/vc,p2/Pc,'k--',lw = lw/2)
axes.plot(v3/vc,p3/Pc,'k--',lw = lw/2,label="isentropes")
# Labels
plt.text(2.8, 0.8, '$s=s_\\tau$',ha= 'center', size = 10)


# criticla point
axes.plot(vc/vc,Pc/Pc,'ko',lw = lw , label = "Critical point")

# Gamma  curve
axes.plot(vg0/vc,pg0/Pc,'b',lw = lw, label = "$\Gamma=0$")


axes.set_ylim([0.5, 1.2])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("LSV.pdf") 