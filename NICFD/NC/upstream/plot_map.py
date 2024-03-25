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
1.3 Double sonic locus
"""
DSL = pd.read_csv("DSL.csv", ",", skiprows=0)

"""
1.4 Isentropy
"""
v = vc*2.5
p = Pc*1.05
s = CP.CoolProp.PropsSI('Smass','P', p, 'Dmass', 1/v,  fluidname)
# S_tau
vt = np.linspace(vc*0.5, vc*4.0,500) # P<Pc
vt = pd.Series(vt)
pt = np.zeros(vt.size) 
for i in vt.index:
    pt[i] = CP.CoolProp.PropsSI('P','Smass', s-176, 'Dmass', 1/vt[i],  fluidname)
# S_vle    
vv = np.linspace(vc*0.5, vc*4.0,500) # P<Pc
vv = pd.Series(vv)
pv = np.zeros(vv.size) 
for i in vv.index:
    pv[i] = CP.CoolProp.PropsSI('P','Smass', s-193, 'Dmass', 1/vv[i],  fluidname)

"""
2. PSmax
"""
vsmax = []
psmax = []
smax = np.linspace(s-193, s-180, 1)
smax = pd.Series(smax)
for k in smax.index:
    """
    2.1 isentropy P,v
    """
    v1 = np.linspace(vc*4.0, vc*1.5,10) # P<Pc
    v1 = pd.Series(v1)
    p1 = np.zeros(v1.size) 
    for i in v1.index:
        p1[i] = CP.CoolProp.PropsSI('P','Smass', smax[k], 'Dmass', 1/v1[i],  fluidname)
    
    """
    2.2 find PSmax at which first post-compression sonic shock is observed when increasing P with max intensity or exhibit Gamma=0
    """
    for j in v1.index:
        d1 = 1/v1[j]
        P1 = p1[j]
        T1 = CP.CoolProp.PropsSI('T','P', P1, 'Dmass', d1,  fluidname)
        c1 = CP.CoolProp.PropsSI('A','P|gas', P1, 'T', T1,  fluidname)
        h1 =  CP.CoolProp.PropsSI('Hmass','P', P1, 'Dmass', d1,  fluidname)
        M1 = 1.1
        u1 = c1*M1
        ht1 = h1 + 0.5*u1*u1
        n2 = 200
        P2 = np.linspace(P1, Pc*0.99 ,n2) # post-shock 
        P2 = pd.Series(P2)
        v2r = np.zeros(P2.size) # v2 for Rayleigh line
        M2 = np.zeros(P2.size) 
        for i in P2.index:
             v2r[i] = (1-(P2[i]-P1)*v1[j]/u1/u1)*v1[j]
             d2 = 1/v2r[i]
             ds = CP.CoolProp.PropsSI('Dmass','P', P2[i], 'Q', 1,  fluidname)
             if v2r[i]<1/ds:
                 # print("two phase")
                 continue
             h2 = CP.CoolProp.PropsSI('Hmass','P', P2[i], 'Dmass', d2,  fluidname)
             c2 = CP.CoolProp.PropsSI('A','P', P2[i], 'Dmass', d2,  fluidname)
             u2 = math.sqrt(2*abs(ht1-h2))
             M2[i] = u2/c2
        i = np.argmin(abs(M2-1.0)) 
        T = CP.CoolProp.PropsSI('T','P', P2[i], 'Dmass', 1/v2r[i],  fluidname)
        G = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P|gas', P2[i], 'T', T,  fluidname)
        print("k, j, i, M2, G2 ", k, j, i, M2[i], G)
        print("P1/Pc, v1/vc, P2/Pc, v1/vc, ", p1[j]/Pc, v1[j]/vc, P2[i]/Pc, v2r[i]/vc)
    break
        # if abs(G)<0.01:
        #     print("k, j, i, M2, G2 ", k, j, i, M2[i], G)
        #     vsmax.append(v1[j])
        #     psmax.append(p1[j])
        #     break




"""
X. plot
"""

fig1 = plt.figure( dpi=300)
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2
"""
X.1 Saturation curve
"""
axes.plot(LSV.iloc[:,2],LSV.iloc[:,3],'r',lw = lw, label = "LVS")


"""
X.2 Gamma = 0
"""
axes.plot(GAMMA.iloc[:,2],GAMMA.iloc[:,3],'k',lw = lw, label = "$\Gamma=0$")

"""
X.3 Isentropy
"""
axes.plot(vt/vc,pt/Pc,'k--',lw = lw/2, label = "$s_\\tau$")
axes.plot(vv/vc,pv/Pc,'k-.',lw = lw/2, label = "$s_{vle}$")
# axes.plot(v1/vc,p1/Pc,'b--',lw = lw/2)

"""
X.4 Double sonic locus
"""

axes.plot(DSL.iloc[:,2], DSL.iloc[:,3],'b',lw = lw, label = "DSL")
# axes.plot(DSL.iloc[:,4], DSL.iloc[:,5],'k',lw = lw)

"""
X.5 PSmax
"""
vsmax = np.array(vsmax)
psmax = np.array(psmax)
# axes.plot(vsmax/vc, psmax/Pc,'b+',lw = lw)



axes.set_ylim([0.5, 1.1])
axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Upstream state map for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("map.pdf") 

end = time.time()
print("computational time(s): ", end - start)