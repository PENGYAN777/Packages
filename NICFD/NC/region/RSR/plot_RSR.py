#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:30:25 2024

plot (P,v) diagram. of rarefaction shock region 

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
import time

start = time.time()

"""
0. Fluid info and shock angle beta, delfction angle theta
"""
# fluidname = "HEOS::MD4M"
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



"""
1. Basic curve
"""
"""
1.1. Saturation curve
"""
LSV = pd.read_csv("LSV.csv", ",", skiprows=0)



"""
1.2 Gamma = 0 curve
"""
GAMMA = pd.read_csv("GAMMA.csv", ",", skiprows=0)


"""
2. Double sonic locus
"""

PAD = []
VAD = []
PBD = []
VBD = []

M1 = 1.0

s1 = np.linspace(507,520,5)
s1 = pd.Series(s1)
for k in s1.index:
    n1 = 50
    v1 = np.linspace(vc*1.2, vc*1.5 ,n1)
    v1 = pd.Series(v1)
    P1 = np.zeros(v1.size)
    v2 = np.zeros(v1.size)
    P2 = np.zeros(v1.size)
    for j in v1.index:
        d1 = 1/v1[j]
        P1[j] = CP.CoolProp.PropsSI('P','Dmass',d1,'Smass',s1[k],fluidname) 
        c1 = CP.CoolProp.PropsSI('A','P', P1[j], 'Dmass', d1,  fluidname)
        h1 =  CP.CoolProp.PropsSI('Hmass','P', P1[j], 'Dmass', d1,  fluidname)
        u1 = c1*M1
        ht1 = h1 + 0.5*u1*u1
        # Rayleigh line
        n2 = 100
        P2r = np.linspace(Pc*0.99, Pc*0.7 ,n2) # post-shock Mach
        P2r = pd.Series(P2r)
        v2r = np.zeros(P2r.size) 
        diff = np.zeros(P2r.size) 
        for i in P2r.index:
             P = P2r[i]
             v2r[i] = (1-(P-P1[j])*v1[j]/u1/u1)*v1[j]
             d2 = 1/v2r[i]
             ds = CP.CoolProp.PropsSI('Dmass','P', P, 'Q', 1,  fluidname)
             if v2r[i]<1/ds:
                 # print("two phase")
                 continue
             c2 = CP.CoolProp.PropsSI('A','P', P, 'Dmass', d2,  fluidname)
             h2 =  CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', d2,  fluidname)
             u2 = c2
             diff[i] = abs(h2 + 0.5*u2*u2 - ht1)/ht1
        i = np.argmin(diff)
        if i>0:
            print("k, j, i, diff: ", k, j, i , diff[i])
            PAD.append(P1[j])
            VAD.append(v1[j])
            PBD.append(P2r[i])
            VBD.append(v2r[i])
            break
    



"""
3. Post sonic shock originating from pre-shock state along saturation curve
"""
PL = []
VL = []

vmin = vc*1.2
vmax = vc*2.26
n3 = 5
v1 = np.linspace(vmin, vmax ,n3)
v1 = pd.Series(v1)
P1 = np.zeros(v1.size)
for j in v1.index:
    print("j: ", j)
    M1 = 1.001
    d1 = 1/v1[j]
    P1[j] = CP.CoolProp.PropsSI('P','Dmass', d1, 'Q', 1,  fluidname)
    c1 = CP.CoolProp.PropsSI('A','P', P1[j], 'Q', 1,  fluidname)
    h1 =  CP.CoolProp.PropsSI('Hmass','P', P1[j], 'Dmass', d1,  fluidname)
    u1 = c1*M1
    ht1 = h1 + 0.5*u1*u1
    # Rayleigh line
    n4 = 100
    P2r = np.linspace(Pc*0.8, Pc*0.75 ,n4) # post-shock Mach
    P2r = pd.Series(P2r)
    v2r = np.zeros(P2r.size) 
    diff = np.zeros(P2r.size) 
    for i in P2r.index:
        diff[i] = 100
        P = P2r[i]
        v2r[i] = (1-(P-P1[j])*v1[j]/u1/u1)*v1[j]
        d2 = 1/v2r[i]
        ds = CP.CoolProp.PropsSI('Dmass','P', P, 'Q', 1,  fluidname)
        if v2r[i]<1/ds:
            # print("two phase")
            continue
        T2 = CP.CoolProp.PropsSI('T','P', P, 'Dmass', d2,  fluidname)
        c2 = CP.CoolProp.PropsSI('A','P|gas', P, 'T', T2,  fluidname)
        h2 =  CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', d2,  fluidname)
        u2 = c2
        diff[i] = abs(h2 + 0.5*u2*u2 - ht1)/ht1
    i = np.argmin(diff)
    if i>0:
        print("LVS i, diff: ",   i , diff[i])
        PL.append(P2r[i])
        VL.append(v2r[i])
    



"""
X. plot
"""

fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2
"""
X.1 Saturation curve
"""
axes.plot(LSV.iloc[:,2],LSV.iloc[:,3],'b',lw = lw, label = "LVS")


"""
X.2 Gamma = 0
"""
axes.plot(GAMMA.iloc[:,2],GAMMA.iloc[:,3],'k--',lw = lw, label = "$\Gamma=0$")

"""
X.3 Isentropy
"""
# axes.plot(v1/vc,P1/Pc,'k--',lw = lw/2, label = "s")

"""
X.4 Double sonic locus
"""
PAD = np.array(PAD)
VAD = np.array(VAD)
PBD = np.array(PBD)
VBD = np.array(VBD)

axes.plot(VAD/vc,PAD/Pc,'k',lw = lw, label = "DSL")
axes.plot(VBD/vc,PBD/Pc,'k',lw = lw)

"""
X.5 Last portion
"""
PL = np.array(PL)
VL = np.array(VL)
axes.plot(VL/vc,PL/Pc,'k+',lw = lw)


# DSL
# axes.plot(vd/vc,Pd/Pc,'ro',lw = lw/2, label = "DSL")

axes.set_ylim([0.5, 1.1])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("RSR.pdf") 

end = time.time()
print("computational time(s): ", end - start)