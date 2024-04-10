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
    
s1 = np.linspace(s-193, s-200, 1)
s1 = pd.Series(s1)
for k in s1.index:
    """
    2.1 P,v along isentropy and within upstream map
    """
    v1 = np.linspace(vc*3.7, vc*1.2,1) #[3.7, 1.2]
    v1 = pd.Series(v1)
    p1 = np.zeros(v1.size) 
    m1 = np.zeros(v1.size) 
    c1 = np.zeros(v1.size) 
    ps = np.zeros(v1.size) 
    vs = np.zeros(v1.size) 
    for j in v1.index:
        p1[j] = CP.CoolProp.PropsSI('P','Smass', s1[k], 'Dmass', 1/v1[j],  fluidname)
        h1 =  CP.CoolProp.PropsSI('Hmass','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        c1[j] =  CP.CoolProp.PropsSI('A','P', p1[j], 'Dmass', 1/v1[j],  fluidname)
        """
        2.2 Rangking Hugoniot curve
        """
        n1 = 100
        v2 = np.linspace(v1[j], vc*1.5 ,n1) # post-shock 
        v2 = pd.Series(v2)
        P2 = np.zeros(v2.size) 
        s2 = np.zeros(v2.size) 
        for l in v2.index:
            # initial guess pf v2
            n2 = 200
            P = np.linspace(p1[j],Pc*1.0 ,n2) 
            P = pd.Series(P)
            diff = np.zeros(P.size) 
            for i in P.index:
                diff[i] = CP.CoolProp.PropsSI('Hmass','P', P[i], 'Dmass', 1/v2[l], fluidname) - h1 - 0.5*(P[i]-p1[j])*(v1[j]+v2[l])
            # print("min index:", np.argmin(abs(diff)), "min diff: ", diff[np.argmin(abs(diff))])
            P2[l] = P[np.argmin(abs(diff))]
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
        print("j, m1, ", j , m1[j])



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

end = time.time()
print("computational time(s): ", end - start)