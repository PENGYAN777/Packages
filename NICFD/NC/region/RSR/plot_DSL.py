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
from newton import f, fp, module
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
2. Double sonic locus
"""

PAD = []
VAD = []
PBD = []
VBD = []

M1 = 1.0

s1 = np.linspace(515,516,1)
s1 = pd.Series(s1)
for k in s1.index:
    n1 = 1
    v1 = np.linspace(vc*1.2, vc*1.5 ,n1)
    v1 = pd.Series(v1)
    P1 = np.zeros(v1.size)
    for j in v1.index:
        d1 = 1/v1[j]
        P1[j] = CP.CoolProp.PropsSI('P','Dmass',d1,'Smass',s1[k],fluidname) 
        c1 = CP.CoolProp.PropsSI('A','P', P1[j], 'Dmass', d1,  fluidname)
        u1 = c1*M1
        """
        2.1 initial guess for newton method
        """
        tol = 1e-3 # torelance
        nmax = 50 # number of iteration
        x = np.zeros( (3,nmax) )
        x0 = (1/vc/2.5, u1/2, Pc*0.7) 
        x[0,0] = x0[0]
        x[1,0] = x0[1]
        x[2,0] = x0[2]
        n = 1
        """
        2.2 while iteration for newton method
        """
        while n<nmax-1:
            a = f(1/v1[j], u1, P1[j], x[0,n-1], x[1,n-1], x[2,n-1])
            b = fp(1/v1[j], u1, P1[j], x[0,n-1], x[1,n-1], x[2,n-1])
            dx = - np.matmul(np.linalg.inv(b), a)
            x[0,n] = x[0, n-1] + dx[0]
            x[1,n] = x[1, n-1] + dx[1]
            x[2,n] = x[2, n-1] + dx[2]
            diff = abs( (module(x[0,n],x[1,n],x[2,n]) - module(x[0,n-1],x[1,n-1],x[2,n-1]) ) 
                       / module(x[0,n-1],x[1,n-1],x[2,n-1]) )
            if diff<tol:
                # print("x1,x2,x3,f(x)", x[0,n], x[1,n],x[2,n], f(1/v1[j], u1, P1[j], x[0,n], x[1,n], x[2,n]) )
                break
            else:
                n = n + 1
        """
        2.3 check M2=1?
        """
        d2 = x[0,n-1]
        u2 = x[1,n-1]
        P2 = x[2,n-1]
        c2 = CP.CoolProp.PropsSI('A','P', P2, 'Dmass', d2,  fluidname)
        M2 = u2/c2
        # if abs(M2-1.0)<1e-3:
        print(" n = ", n, "diff = ", diff, "M2 = ", M2)
        PAD.append(P1[j])
        VAD.append(v1[j])
        PBD.append(P2)
        VBD.append(1/d2)
        break
        
        

    

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
axes.plot(LSV.iloc[:,2],LSV.iloc[:,3],'b',lw = lw, label = "LVS")


"""
X.2 Gamma = 0
"""
axes.plot(GAMMA.iloc[:,2],GAMMA.iloc[:,3],'k--',lw = lw, label = "$\Gamma=0$")

"""
X.3 Isentropy
"""
# axes.plot(v1/vc,P1/Pc,'k--',lw = lw/2, label = "s")

"""
X.4 Double sonic locus
"""
PAD = np.array(PAD)
VAD = np.array(VAD)
PBD = np.array(PBD)
VBD = np.array(VBD)

axes.plot(VAD/vc,PAD/Pc,'ko',lw = lw, label = "pre")
axes.plot(VBD/vc,PBD/Pc,'k+',lw = lw,label = "post")


# DSL
# axes.plot(vd/vc,Pd/Pc,'ro',lw = lw/2, label = "DSL")

axes.set_ylim([0.5, 1.1])
# axes.set_xlim([0.5, 4.0])
axes.set_xlabel('$v/v_c$')
axes.set_ylabel('$P/P_c$')
plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("DSL.pdf") 

"""
XX write into csv file
"""
pd.DataFrame(VAD/vc).to_csv('DSL.csv', index_label = "Index", header  = ['VAD/vc']) 
data = pd.read_csv("DSL.csv", ",")
D =pd.DataFrame({'PAD/Pc': PAD/Pc, 'VBD/vc': VBD/vc, 'PBD/Pc': PBD/Pc, })
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("DSL.csv")

end = time.time()
print("computational time(s): ", end - start)
