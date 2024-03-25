#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 09:51:36 2024

Write LVS results.

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
1. Saturation curve
"""
psl = np.linspace(Pc*0.5, Pc*1.0, 500)
dsl = CP.CoolProp.PropsSI('Dmass','P',psl,'Q',0,fluidname)
vsl = 1/dsl
# vapor phase
psv = np.linspace(Pc*1.0, Pc*0.5, 500)
dsv = CP.CoolProp.PropsSI('Dmass','P',psv,'Q',1,fluidname)
vsv = 1/dsv

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
axes.plot(vsl/vc,psl/Pc,'k',lw = lw, label = "LVS")
axes.plot(vsv/vc,psv/Pc,'k',lw = lw)

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
vs = np.concatenate((vsl/vc, vsv/vc ))
ps = np.concatenate((psl/Pc, psv/Pc ))
pd.DataFrame(vs).to_csv('LSV.csv', index_label = "Index", header  = ['vs/vc']) 
data = pd.read_csv("LSV.csv", ",")
D =pd.DataFrame({'ps/Pc': ps, })
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("LSV.csv")

end = time.time()
print("computational time(s): ", end - start)


