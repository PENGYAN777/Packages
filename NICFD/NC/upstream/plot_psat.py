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

# smin
vm = np.linspace(vc*0.5, vc*4.0,500) # P<Pc
vm = pd.Series(vm)
pm = np.zeros(vm.size) 
for i in vm.index:
    pm[i] = CP.CoolProp.PropsSI('P','Smass', s-218, 'Dmass', 1/vm[i],  fluidname)

"""
2. PSat
"""

        
    

    





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
# axes.plot(vt/vc,pt/Pc,'k--',lw = lw/2, label = "$s_\\tau$")
# axes.plot(vv/vc,pv/Pc,'k-.',lw = lw/2, label = "$s_{vle}$")
axes.plot(vt/vc,pt/Pc,'k--',lw = lw/2)
axes.plot(vv/vc,pv/Pc,'k--',lw = lw/2)
axes.plot(vm/vc,pm/Pc,'k-.',lw = lw/2)


"""
X.4 Double sonic locus
"""

axes.plot(DSL.iloc[:,2], DSL.iloc[:,3],'b',lw = lw, label = "DSL")
# axes.plot(DSL.iloc[:,4], DSL.iloc[:,5],'k',lw = lw)

"""
X.5 PSmax
"""
# vsmax = np.array(vsmax)
# psmax = np.array(psmax)
# axes.plot(v1/vc, p1/Pc,'bo',lw = lw)
# axes.plot(v2/vc, P2/Pc,'b+',lw = lw)
# axes.plot(vsmax/vc, psmax/Pc,'bo',lw = lw, label = "PSmax")



axes.set_ylim([0.6, 1.1])
axes.set_xlim([0.7, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Upstream state map for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
# fig1.savefig("psmax.pdf") 

"""
X.1 write into csv file
"""
# pd.DataFrame(vsmax/vc).to_csv('PSmax.csv', index_label = "Index", header  = ['vsmax/vc']) 
# data = pd.read_csv("PSat.csv", ",")
# D =pd.DataFrame({'psat/Pc': psmax/Pc, })
# newData = pd.concat([data, D], join = 'outer', axis = 1)
# newData.to_csv("PSat.csv")

end = time.time()
print("computational time(s): ", end - start)