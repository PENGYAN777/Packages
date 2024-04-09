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

"""
2. PSmax
"""
vsmax = []
psmax = []
smax = np.linspace(s-193, s-177, 6)
# smax = [s-193,s-190,  ]
# smax = [s-187, s-181 ]
# smax = [s-178, s-177 ]

smax = pd.Series(smax)
for k in smax.index:
    """
    2.1 P,v along isentropy 
    """
    if k<1:
        v1 = np.linspace(vc*4.0, vc*3.0,10)
    elif k == 1:
        continue
    elif k<4:
        v1 = np.linspace(vc*3.5, vc*2.5,10)
    else:
        v1 = np.linspace(vc*2.5, vc*2.0,10) 
    v1 = pd.Series(v1)
    p1 = np.zeros(v1.size) 
    diffs = np.zeros(v1.size) 
    for j in v1.index:
        p1[j] = CP.CoolProp.PropsSI('P','Smass', smax[k], 'Dmass', 1/v1[j],  fluidname)
        c1 =  CP.CoolProp.PropsSI('A','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        h1 =  CP.CoolProp.PropsSI('Hmass','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        s1 = smax[k]
        """
        2.2 find PSmax 
        """
        """
        2.2.1 Rangking Hugoniot curve
        """
        n1 = 100
        P2 = np.linspace(p1[j], Pc*0.99 ,n1) # post-shock 
        P2 = pd.Series(P2)
        v2 = np.zeros(P2.size) 
        G2 = np.zeros(P2.size) 
        s2 = np.zeros(P2.size) 
        for l in P2.index:
            P = P2[l]
            # initial guess pf v2
            n2 = 100
            v = np.linspace(vc*1.5,vc*4.0 ,n2) 
            v = pd.Series(v)
            diff = np.zeros(v.size) 
            for i in v.index:
                diff[i] = CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', 1/v[i], fluidname) - h1 - 0.5*(P-p1[j])*(v1[j]+v[i])
            # print("min index:", np.argmin(abs(diff)), "min diff: ", diff[np.argmin(abs(diff))])
            v2[l] = v[np.argmin(abs(diff))]
            G2[l] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P2[l],'Dmass', 1/v2[l],fluidname)
            s2[l] = CP.CoolProp.PropsSI('Smass','P', P2[l], 'Dmass', 1/v2[l],  fluidname)
            if G2[l]<0:
                splus1 = s2[np.argmax(abs(s2))]        
                splus2 = s2[np.argmin((G2))]    
                diffs[j] = abs(splus1-splus2)/splus1 
                break
    i = np.argmin(diffs)
    vsmax.append(v1[i])
    psmax.append(p1[i])
    print("k,i, diff ", k , i, diffs[i] )
        

        
    

    





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
axes.plot(vv/vc,pv/Pc,'k-.',lw = lw/2)


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
# axes.plot(v1/vc, p1/Pc,'bo',lw = lw)
# axes.plot(v2/vc, P2/Pc,'b+',lw = lw)
axes.plot(vsmax/vc, psmax/Pc,'bo',lw = lw, label = "PSmax")



axes.set_ylim([0.6, 1.1])
axes.set_xlim([0.7, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Upstream state map for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("psmax.pdf") 

"""
X.1 write into csv file
"""
pd.DataFrame(vsmax/vc).to_csv('PSmax.csv', index_label = "Index", header  = ['vsmax/vc']) 
data = pd.read_csv("PSmax.csv", ",")
D =pd.DataFrame({'psmax/Pc': psmax/Pc, })
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("PSmax.csv")

end = time.time()
print("computational time(s): ", end - start)