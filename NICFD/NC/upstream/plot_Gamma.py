#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:04:05 2024

plot Gamma 

@author: yan
"""

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
0. Fluid info 
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
1. Gamma = 0
"""

# consider single vapor phase only
pg = np.linspace(Pc*0.7, Pc*1.0,100) #
dg = CP.CoolProp.PropsSI('Dmass','P',pg,'Q',1,fluidname) 
vg = 1/dg
vg = pd.Series(vg)
pg = pd.Series(pg)
pg0 = [] # p w.r.t Gamma = 0
vg0 = []
for i in vg.index:
    pgl = np.linspace(pg[i]*1.001, Pc*1.1,1000) # list all possible p
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
axes.plot(vg0/vc,pg0/Pc,'bo',lw = lw, label = "$\Gamma=0$")

axes.set_ylim([0.5, 1.1])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
# fig1.savefig("RSR.pdf") 

"""
X.1 write into csv file
"""
pd.DataFrame(vg0/vc).to_csv('GAMMA.csv', index_label = "Index", header  = ['vg0/vc']) 
data = pd.read_csv("GAMMA.csv", ",")
D =pd.DataFrame({'pg0/Pc': pg0/Pc, })
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("GAMMA.csv")

end = time.time()
print("computational time(s): ", end - start)
