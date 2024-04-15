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
from scipy.interpolate import griddata
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
DSL = DSL[::-1]

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
1.5 PSmax
"""
smax = pd.read_csv("PSmax.csv", ",", skiprows=0)

"""
1.6 PSat
"""
sat = pd.read_csv("PSat.csv", ",", skiprows=0)

    
"""
2. Mach transion
"""
l1 = 4
l2 = 2
V1 =np.zeros((l1,l2)) 
P1 =np.zeros((l1,l2))
M1 =np.zeros((l1,l2))  
s1 = np.array([s-200, s-193, s-185, s-178,])
s1 = pd.Series(s1)
for k in s1.index:
    count = 0 # record the type of shock, classical 0 or non-classical 1 
    """
    2.1 P,v along isentropy and within upstream map
    """
    #[3.0, 2.5] for s-200
    #[3.7, 1.2] for s-193
    #[2.0, 1.4] for s-180
    if k == 0:
        v1 = np.linspace(vc*3.0, vc*2.5, l2) 
    if k == 1:
        v1 = np.linspace(vc*3.7, vc*1.2, l2) 
    if k == 2:
        v1 = np.linspace(vc*2.4, vc*1.26, l2) 
    if k == 3:
        v1 = np.linspace(vc*2.2, vc*1.4, l2) 
    v1 = pd.Series(v1)
    p1 = np.zeros(v1.size) 
    m1 = np.zeros(v1.size) 
    c1 = np.zeros(v1.size) 
    ps = np.zeros(v1.size) 
    vs = np.zeros(v1.size) 
    for j in v1.index:
        p1[j] = CP.CoolProp.PropsSI('P','Smass', s1[k], 'Dmass', 1/v1[j],  fluidname)
        h1 =  CP.CoolProp.PropsSI('Hmass','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        phase1 = CP.CoolProp.PropsSI('Phase','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        if phase1 == 6:
            continue
        c1[j] =  CP.CoolProp.PropsSI('A','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        G1 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',p1[j],'Dmass', 1/v1[j],fluidname)
        if G1<0:
            count = 1
        """
        2.2 Rangking Hugoniot curve
        """
        n1 = 200
        if count == 0:
            v2 = np.linspace(v1[j], vc*1.5 ,n1) # classical shock
        if count == 1:
            v2 = np.linspace(v1[j], vc*3.7 ,n1) # non-classical shock
        v2 = pd.Series(v2)
        P2 = np.zeros(v2.size) 
        s2 = np.zeros(v2.size) 
        G2 = np.zeros(v2.size) 
        for l in v2.index:
            # initial guess pf v2
            n2 = 200
            if count == 0:
                P = np.linspace(p1[j],Pc*1.0 ,n2) 
            if count == 1:
                P = np.linspace(p1[j],Pc*0.6 ,n2) 
            P = pd.Series(P)
            diff = np.zeros(P.size) 
            for i in P.index:
                diff[i] = CP.CoolProp.PropsSI('Hmass','P', P[i], 'Dmass', 1/v2[l], fluidname) - h1 - 0.5*(P[i]-p1[j])*(v1[j]+v2[l])
            # print("min index:", np.argmin(abs(diff)), "min diff: ", diff[np.argmin(abs(diff))])
            P2[l] = P[np.argmin(abs(diff))]
            phase2 = CP.CoolProp.PropsSI('Phase','P', P2[l], 'Dmass', 1/v2[l],  fluidname)
            if phase2 == 6: # two phase
                G2[l] = 0
            else:
                G2[l] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P2[l],'Dmass', 1/v2[l],fluidname)
            if count == 0:
                if G2[l]<0:
                    s2[l] = CP.CoolProp.PropsSI('Smass','P', P2[l], 'Dmass', 1/v2[l],  fluidname)
            if count == 1:
                if G2[l]>0:
                    s2[l] = CP.CoolProp.PropsSI('Smass','P', P2[l], 'Dmass', 1/v2[l],  fluidname)
        """
        2.3 find Smax
        """
        i = np.argmax(abs(s2))
        ps[j] = P2[i]
        vs[j] = v2[i]
        hs =  CP.CoolProp.PropsSI('Hmass','P', ps[j], 'Dmass', 1/vs[j],  fluidname)
        cs = CP.CoolProp.PropsSI('A','P', ps[j], 'Dmass', 1/vs[j],  fluidname)
        ms = 1.0
        us = cs*ms
        hts = hs + 0.5*us*us
        u1 = np.sqrt( 2*(hts-h1) )
        m1[j] = u1/c1[j]
        # m1[j] = v1[j]/c1[j]*np.sqrt( (ps-p1[j]) / (v1[j]-vs) )
        print("k, j, m1, ", k, j , m1[j])
        V1[k,j] = v1[j]
        P1[k,j] = p1[j]
        M1[k,j] = m1[j]


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
axes.plot(vm/vc,pm/Pc,'k--',lw = lw/2)


"""
X.4 DSL + PSmax
"""
bv = np.concatenate( ( np.array([GAMMA.iloc[0,2]]), sat.iloc[:,2], 
                      smax.iloc[:,2], DSL.iloc[:,2], np.array([1.19]),  ) )
bp = np.concatenate( ( np.array([GAMMA.iloc[0,3]]), sat.iloc[:,3],
                      smax.iloc[:,3], DSL.iloc[:,3], np.array([0.99]),  ) )

axes.plot(bv, bp, 'b',lw = lw, label = "boundary")

"""
X.5 Mach number
"""

axes.plot(v1/vc,p1/Pc,'ko',lw = lw/2)
axes.plot(vs/vc,ps/Pc,'k+',lw = lw/2)

"""
X.6 Contour of Mach number
"""
levels = [1.0, 1.05, 1.10, 1.15, 1.2]
plt.contourf(V1/vc, P1/Pc, M1, levels,  cmap='rainbow', linewidths=lw/10)
plt.colorbar()

axes.set_ylim([0.6, 1.1])
axes.set_xlim([0.7, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Upstream state map for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("map.pdf") 

"""
X.1 write into csv file
"""
# pd.DataFrame(v1/vc).to_csv('contour.csv', index_label = "Index", header  = ['v1/vc']) 
# data = pd.read_csv("contour.csv", ",")
# D =pd.DataFrame({'p1/Pc': p1/Pc, 'Mtr': m1, })
# newData = pd.concat([data, D], join = 'outer', axis = 1)
# newData.to_csv("contour.csv")

end = time.time()
print("computational time(s): ", end - start)