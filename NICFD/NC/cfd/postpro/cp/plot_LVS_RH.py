#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:30:25 2024

plot (P,v) diagram. 

@author: yan
"""
# G = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P,'T', T,fluidname)
# P = CP.CoolProp.PropsSI('P','T', T, 'Q', 1,  fluidname)

import matplotlib
import math
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import pandas as pd
import time

start = time.time()

"""
0. Fluid info and shock angle beta, delfction angle theta
"""
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

beta_cfd = math.atan((0.72-0.2)/0.7)
# beta_cfd = beta_cfd*180/math.pi # degree
theta = math.atan(0.2/0.7)

"""
1. Saturation curve
"""
# liquid phase
psl = np.linspace(Pc*0.5, Pc*1.0, 500)
dsl = CP.CoolProp.PropsSI('Dmass','P',psl,'Q',0,fluidname)
vsl = 1/dsl
# vapor phase
psv = np.linspace(Pc*0.5, Pc*1.0, 500)
dsv = CP.CoolProp.PropsSI('Dmass','P',psv,'Q',1,fluidname)
vsv = 1/dsv


"""
2. Gamma=0 curve
"""
# consider single capor phase only
pg = np.linspace(Pc*0.7, Pc*1.0,100) #
dg = CP.CoolProp.PropsSI('Dmass','P',pg,'Q',1,fluidname) 
vg = 1/dg
vg = pd.Series(vg)
pg = pd.Series(pg)
pg0 = [] # p w.r.t Gamma = 0
vg0 = []
for i in vg.index:
    pgl = np.linspace(pg[i]*1.005, Pc*1.05,100) # list all possible p
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
3. Ranking Hugoniot curve
"""
"""
3.1 pre-shock state
"""

P1 = Pc*0.99
d1 = dc*0.80
v1 = 1/d1
# v1 = 1/d1
c1 = CP.CoolProp.PropsSI('A','P', P1, 'Dmass', d1,  fluidname)
h1 =  CP.CoolProp.PropsSI('Hmass','P', P1, 'Dmass', d1,  fluidname)
M1 = 1.7
u1 = M1*c1
ht1 = h1 + 0.5*u1*u1 
"""
3.2 Ranking Hugoniot curve h(p,v)
"""
n1 = 30
P2 = np.linspace(Pc*0.7, Pc*0.99 ,n1) # post-shock Mach
P2 = pd.Series(P2)
v2 = np.zeros(P2.size) 
M2 = np.zeros(P2.size) 
c2 = np.zeros(P2.size) 
u2 = np.zeros(P2.size) 
h2 = np.zeros(P2.size) 
for j in P2.index:
    P = P2[j]
    # initial guess pf v2
    n2 = 300
    v = np.linspace(vc*1.2,vc*4.0 ,n2) 
    v = pd.Series(v)
    diff = np.zeros(v.size) 
    for i in v.index:
        diff[i] = CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', 1/v[i], fluidname) - h1 - 0.5*(P-P1)*(v1+v[i])
    # print("min index:", np.argmin(abs(diff)), "min diff: ", diff[np.argmin(abs(diff))])
    v2[j] = v[np.argmin(abs(diff))]
    d2 = 1/v2[j]
    c2[j] = CP.CoolProp.PropsSI('A','P', P, 'Dmass', d2,  fluidname)
    h2[j] = CP.CoolProp.PropsSI('Hmass','P', P, 'Dmass', d2,  fluidname)
    u2[j] = math.sqrt(2*(ht1 - h2[j]))
    M2[j] = u2[j]/c2[j]
"""
3.3 Rayleigh line
"""

u1n = u1*math.sin(beta_cfd)
u1t = u1*math.cos(beta_cfd)
v2r = np.zeros(P2.size) # v2 for Rayleigh line
for i in P2.index:
      P = P2[i]
      v2r[i] = v1-(P-P1)*v1*v1/u1n/u1n
"""
3.4 Post-shock state
"""

PP = Pc*0.887
dP = dc*0.522
cP = CP.CoolProp.PropsSI('A','P', PP, 'Dmass', dP,  fluidname)
hP = CP.CoolProp.PropsSI('Hmass','P', PP, 'Dmass', dP,  fluidname)
uP = math.sqrt(2*(ht1 -hP))
MP = uP/cP
print("M2: ", MP)

"""
3.5 isentrope
"""

s1 = CP.CoolProp.PropsSI('Smass','P', P1, 'Dmass', d1,  fluidname)
s2 = CP.CoolProp.PropsSI('Smass','P', PP, 'Dmass', dP,  fluidname)
ve1 = np.linspace(vc*1.5, vc*3,500) # P<Pc
ve1 = pd.Series(ve1)
pe1 = np.zeros(ve1.size) 
pe2 = np.zeros(ve1.size) 
for i in ve1.index:
    pe1[i] = CP.CoolProp.PropsSI('P','Smass', s1+50, 'Dmass', 1/ve1[i],  fluidname)
    pe2[i] = CP.CoolProp.PropsSI('P','Smass', s2, 'Dmass', 1/ve1[i],  fluidname)

"""
x. plot
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2
"""
x.1 Liquid vapor saturation curve
"""
axes.plot(vsl/vc,psl/Pc,'k',lw = lw, label = "LVS")
axes.plot(vsv/vc,psv/Pc,'k',lw = lw)

"""
x.2 critical point
"""
# axes.plot(vc/vc,Pc/Pc,'ko',lw = lw , label = "Critical point")

"""
x.3 Gamma = 0 curve
"""
axes.plot(vg0/vc,pg0/Pc,'g',lw = lw, label = "$\Gamma=0$")

"""
x.4 initial point of Ranking Hugonoit curve 
"""
axes.plot(1/d1/vc,P1/Pc,'bo',lw = lw/2, label = "Pre-shock")

"""
x.5 Ranking Hugonoit curve 
"""
axes.plot(v2/vc,P2/Pc,'b',lw = lw, label = "Ranking Hugoniot")

"""
x.6 Rayleigh line
"""
axes.plot(v2r/vc,P2/Pc,'b--',lw = lw, label = "Rayleigh line")

"""
x.7 post-shock state
"""
axes.plot(1/dP/vc,PP/Pc,'b+',lw = lw/2, label = "Post-shock")

"""
x.8 isentrope
"""
axes.plot(ve1/vc,pe1/Pc,'k--',lw = lw/2,label="isentropes")
# axes.plot(ve1/vc,pe2/Pc,'k--',lw = lw/2)

# axes.set_ylim([0.5, 1.2])
axes.set_xlim([0.5, 5.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("LSV_RHL.pdf") 


end = time.time()
print("computational time(s): ", end - start)