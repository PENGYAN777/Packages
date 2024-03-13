#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:21:53 2024

Compute rarefaction oblique shock angle. 
For compression shock, the algorithm is slightly different, plot the geometry and see un,nt,beta relation

@author: yan
"""
import math

beta_cfd = math.atan(0.52/0.7) 
beta_cfd = beta_cfd*180/math.pi
print("Shock angle (degree) from CFD: ", beta_cfd)

from scipy.optimize import fsolve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import CoolProp as CP

"""
0. FLuid property
"""
# print("--------------------Fluid info ------------------")
fluidname = "MD4M"
print("Fluid name:", fluidname)
R = CP.CoolProp.PropsSI("gas_constant",fluidname)
W = CP.CoolProp.PropsSI("molar_mass",fluidname)
Rs = R/W
print("spefific ags constant: J/Kg/K", Rs)
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
print("critical temperature[K]:", Tc)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
print("critical pressure[Pa]:", Pc)
dc = CP.CoolProp.PropsSI('Dmass','P',Pc,'T',Tc,fluidname) 


"""
1. pre-shock conditions
"""

P1 = 820273
d1 = 209.5
T1 = CP.CoolProp.PropsSI('T','P',P1,'Dmass',d1,fluidname) 
c1 = CP.CoolProp.PropsSI('A','P',P1,'T',T1,fluidname) 
h1 = CP.CoolProp.PropsSI('Hmass','P',P1,'T',T1,fluidname) 
M1 = 1.7
u1 = M1*c1 
ht1 = h1 + 0.5*u1*u1
print("pre-shock ocnditions: P1,T1,D1,M1 ",P1,T1,d1,M1 )
theta = math.atan(0.2/0.7)
print("deflection angle ", theta*180/math.pi)
"""
2. post-shock states
"""

beta =np.linspace(1,80,1000)*math.pi/180 
beta = pd.Series(beta)
p = np.zeros(beta.size) 
d = np.zeros(beta.size) 
u = np.zeros(beta.size) 
h = np.zeros(beta.size) 
diff = np.zeros(beta.size) 
for i in beta.index:
    diff[i] = 100
    u[i] = u1*math.cos(beta[i])/math.cos(beta[i] + theta)
    if u[i]>0:
        d[i] = d1*u1*math.sin(beta[i])/u[i]/math.sin(beta[i] + theta)
    if d[i]>0:
        p[i] = P1 + d1*u1*u1*math.sin(beta[i])*math.sin(beta[i]) - d[i]*u[i]*u[i]*math.sin(beta[i]+theta)*math.sin(beta[i]+theta) 
    if p[i]>0:     
        diff[i] = ht1 - CP.CoolProp.PropsSI('Hmass','P',p[i],'Dmass',d[i],fluidname) -0.5*u[i]*u[i]
i = np.argmin(abs(diff))
print("min index:", i)
P2 = p[i]
d2 = d[i]
T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',d2,fluidname) 
c2 = CP.CoolProp.PropsSI('A','P',P2,'Dmass',d2,fluidname) 
u2 = u[i]
M2 = u2/c2
beta_theory = beta[i]
print("post-shock ocnditions: P2,T2,D2,M2 ",P2,T2,d2,M2 )
print("Shock angle (degree) from theory: ", beta_theory*180/math.pi)



