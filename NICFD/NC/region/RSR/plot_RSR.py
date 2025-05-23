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
import os
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')

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
2. Last portion of RSR
"""
LSP = pd.read_csv("LSP.csv", ",", skiprows=0)

"""
X. plot
"""

fig1 = plt.figure( dpi=300)
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2.5
"""
X.1 Saturation curve
"""
axes.plot(LSV.iloc[:,2],LSV.iloc[:,3],'r',lw = lw, label = "LVS")


"""
X.2 Gamma = 0
"""
axes.plot(GAMMA.iloc[:,2],GAMMA.iloc[:,3],'b',lw = lw, label = "$\Gamma=0$")

"""
X.3 Isentropy
"""
# axes.plot(v1/vc,P1/Pc,'k--',lw = lw/2, label = "s")

"""
X.4 Double sonic locus
"""
bv = np.concatenate( (  np.array([1.19]), DSL.iloc[:,2] ) )
bp = np.concatenate( ( np.array([0.99]), DSL.iloc[:,3] ) )
axes.plot(bv, bp,'k',lw = lw, label = "DSL")
axes.plot(DSL.iloc[:,4], DSL.iloc[:,5],'k',lw = lw)

"""
X.5 Last portion
"""

axes.plot(LSP.iloc[:,2], LSP.iloc[:,3], 'k',lw = lw)

"""
X.6 pre and post shock condition
"""
P1 = 0.99*Pc
d1 = 0.80*dc
v1 = 1/d1
axes.plot(v1/vc, P1/Pc, 'k+',lw = lw, label = "Pre")

P2 = 0.887*Pc
d2 = 0.522*dc
v2 = 1/d2
axes.plot(v2/vc, P2/Pc, 'k*',lw = lw, label = "Post")


axes.set_ylim([0.5, 1.1])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("RSR.pdf") 

end = time.time()
print("computational time(s): ", end - start)