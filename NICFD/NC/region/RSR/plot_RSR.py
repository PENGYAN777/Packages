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
"""
2.1 Find M2 = 1 from Rayleigh line
"""
M1 = 1.0
s1 = 515
v1 = 1.35*vc
d1 = 1/v1
P1 = CP.CoolProp.PropsSI('P','Dmass',d1,'Smass',s1,fluidname) 
c1 = CP.CoolProp.PropsSI('A','P', P1, 'Dmass', d1,  fluidname)
h1 =  CP.CoolProp.PropsSI('Hmass','P', P1, 'Dmass', d1,  fluidname)
u1 = c1*M1
ht1 = h1 + 0.5*u1*u1
# Rayleigh line
n1 = 100
P2r = np.linspace(P1*0.99, Pc*0.85 ,n1) # post-shock Mach
P2r = pd.Series(P2r)
v2r = np.zeros(P2r.size) 
M2 = np.zeros(P2r.size) 
for i in P2r.index:
     P = P2r[i]
     v2r[i] = (1-(P-P1)*v1/u1/u1)*v1
     d2 = 1/v2r[i]
     c2 = CP.CoolProp.PropsSI('A','P', P, 'Dmass', d2,  fluidname)
     h2 =  CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', d2,  fluidname)
     u2 = (ht1-h2)*2
     u2 = math.sqrt(u2)
     M2[i] = u2/c2 
i = np.argmin(abs(M2-1.0)) 
v2 = v2r[i]
P2 = P2r[i]
print("M1 > M2: ", M1,M2[i])



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
axes.plot(LSV.iloc[:,2],LSV.iloc[:,3],'k',lw = lw, label = "LVS")


"""
X.2 Gamma = 0
"""
axes.plot(GAMMA.iloc[:,2],GAMMA.iloc[:,3],'b',lw = lw, label = "$\Gamma=0$")

"""
X.3 Isentropy
"""
axes.plot(v1/vc,P1/Pc,'k--',lw = lw/2, label = "s")

"""
X.4 Rayleigh line
"""
axes.plot(v1/vc,P1/Pc,'bo',lw = lw, label = "Pre")
axes.plot(v2/vc,P2/Pc,'b+',lw = lw, label = "Post")

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