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
from bisection import bisection_method
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
x = np.zeros((2,100))
y = np.zeros((2,100))

PAD = []
VAD = []
PBD = []
VBD = []

M1 = 1.0
count = 0
s1 = np.linspace(510,520,1)
s1 = pd.Series(s1)
for k in s1.index:
    n1 = 2
    v1 = np.linspace(vc*1.2, vc*1.3 ,n1)
    v1 = pd.Series(v1)
    P1 = np.zeros(v1.size)
    for j in v1.index:
        d1 = 1/v1[j]
        P1[j] = CP.CoolProp.PropsSI('P','Dmass',d1,'Smass',s1[k],fluidname) 
        G1 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P1[j],'Dmass', d1,fluidname)
        c1 = CP.CoolProp.PropsSI('A','P', P1[j], 'Dmass', d1,  fluidname)
        u1 = c1*M1
        """
        2.1 define function since f is alwasy positive, so can not use bisection method
        """
        def ns(d2,d1,u1,p1):
            h1 = CP.CoolProp.PropsSI('Hmass','P', p1, 'Dmass', d1,  fluidname)
            ht1 = h1 + 0.5*u1*u1
            u2 = d1*u1/d2
            p2 = d1*u1*u1 + p1 - d2*u2*u2
            h2 = CP.CoolProp.PropsSI('Hmass','P', p2, 'Dmass', d2,  fluidname)
            ht2 = h2 + 0.5*u2*u2 
            return (ht1 - ht2)/ht1
        # a = d1/1.1
        # b = d1/5
        # root, iterations = bisection_method(ns, d1,u1,P1[j], a, b)
        # print("Approximate root:", root)
        # print("Number of iterations:", iterations)
        D2 = np.linspace(d1/0.5, d1/3 ,100)
        f = ns(D2, d1, u1, P1[j])
        i = np.argmin(abs(f)) 
        x[j,:] =  D2
        y[j,:] =  f
        """
        2.2 check M2=1?
        """
        d2 = D2[i]
        u2 = d1*u1/d2
        P2 = d1*u1*u1 + P1[j] - d2*u2*u2
        c2 = CP.CoolProp.PropsSI('A','P', P2, 'Dmass', d2,  fluidname)
        M2 = u2/c2
        G2 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics', 'P',P2,'Dmass', d2,fluidname)
        if abs(M2-1.0)<1e-3 and G2>0:
            count = count + 1
            PAD.append(P1[j])
            VAD.append(v1[j])
            PBD.append(P2)
            VBD.append(1/d2)
            print("k, ", k, "diff = ", f[i], "M2 = ", M2)
            break
    
print("count = ", count)     
        
        

    

"""
X. plot
"""

fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
lw = 2


axes.plot(1/x[0,:]/vc, y[0,:],'b',lw = lw, label = "v1/vc=1.2, s=510")
axes.plot(v1[0]/vc,0,'bo',lw = lw,)

axes.plot(1/x[1,:]/vc, y[1,:],'k',lw = lw, label = "v1/vc=1.3, s=510")
axes.plot(v1[1]/vc,0,'ko',lw = lw,)

# axes.set_ylim([-0.0, 0.1])
# axes.set_xlim([1.0, 6.0])
axes.set_xlabel('$v/vc$')
axes.set_ylabel('$f(v)$')
# plt.title('Pv diagram for siloxane MD$_4$M')
axes.legend(loc=0) # 2 means left top
fig1.savefig("f.pdf") 

"""
XX write into csv file
"""
# pd.DataFrame(VAD/vc).to_csv('DSL.csv', index_label = "Index", header  = ['VAD/vc']) 
# data = pd.read_csv("DSL.csv", ",")
# D =pd.DataFrame({'PAD/Pc': PAD/Pc, 'VBD/vc': VBD/vc, 'PBD/Pc': PBD/Pc, })
# newData = pd.concat([data, D], join = 'outer', axis = 1)
# newData.to_csv("DSL.csv")

end = time.time()
print("computational time(s): ", end - start)
