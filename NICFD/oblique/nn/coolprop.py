#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute oblique shock relatiomnships for CoolProp EOS
input: upstream P,T,M.

"""

import math
import numpy as np
import pandas as pd
import CoolProp as CP
from sympy import sin, cos, tan


"""
0. fluid property
"""
fluidname = "nitrogen"
print("Fluid name:", fluidname)
R = CP.CoolProp.PropsSI("gas_constant",fluidname)
print("universal gas constant:  J/mol/K", R)
W = CP.CoolProp.PropsSI("molar_mass",fluidname)
print("molar mass: kg/mol", W)
Rs = R/W
print("spefific ags constant: J/Kg/K", Rs)
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
print("critical temperature[K]:", Tc)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
print("critical pressure[Pa]:", Pc)
dc = CP.CoolProp.PropsSI('Dmass','P',Pc,'T',Tc,fluidname) 
"""
1. input pre-shock conditions
"""


P1 = Pc*0.4# pre-shock pressure
T1 = 300 # pre-shcok temperature
d1 = CP.CoolProp.PropsSI('Dmass','P',P1,'T',T1,fluidname) # pre-shock density
Z1 = CP.CoolProp.PropsSI('Z','P',P1,'T',T1,fluidname) 
G1 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',P1,'T',T1,fluidname) 
print("Z,Gamma:", Z1,G1)
s1 = CP.CoolProp.PropsSI('Smass','P',P1,'T',T1,fluidname) 
cv = CP.CoolProp.PropsSI('Cvmass','P',P1,'T',T1,fluidname) 
cp = CP.CoolProp.PropsSI('Cpmass','P',P1,'T',T1,fluidname) 
print("heat ratio:", cp/cv)
M1 = 2
print("min shock angle (degree): ", np.arcsin(1/M1)*180/math.pi)
c1 = CP.CoolProp.PropsSI('A','P',P1,'T',T1,fluidname)
U1 = M1*c1
beta = 39.4*math.pi/180 # shock angle
h1 = CP.CoolProp.PropsSI('Hmass','P',P1,'T',T1,fluidname) 
ht1 = h1 + 0.5*U1*U1
print("ht:", ht1)
# compute total pressure/temperature
Pt1 = CP.CoolProp.PropsSI('P','Smass',s1,'Hmass',ht1,fluidname) 
Tt1 = CP.CoolProp.PropsSI('T','Smass',s1,'Hmass',ht1,fluidname) 
htotal1 = CP.CoolProp.PropsSI('Hmass','P',Pt1,'T',Tt1,fluidname) 
print("ht1:", htotal1)

"""
2. find P2
"""

p2 = np.linspace(P1*1.01,P1*3,1000) # post-shock pressure 
p2 = pd.Series(p2)
u2 = np.zeros(p2.size) 
d2 = np.zeros(p2.size) 
diff1 = np.zeros(p2.size)
for i in p2.index:
    if abs(p2[i]-Pc)<0.01*Pc:
        p2[i] = 0.99*Pc
    u2[i] = (P1+d1*U1*U1*sin(beta)*sin(beta)-p2[i])/d1/U1/sin(beta)   
    if u2[i]>0:
        d2[i] =  d1*U1*sin(beta)/u2[i]
        diff1[i] = ht1 - 0.5*u2[i]*u2[i] -0.5*U1*U1*cos(beta)*cos(beta) - CP.CoolProp.PropsSI('Hmass','P',p2[i],'Dmass',d2[i],fluidname) 
diff3 = abs(diff1)
minval = np.min(diff3[np.nonzero(diff3)]) 
diff1 = abs(diff1) - minval   
U2 = u2[np.argmin(abs(diff1))]**2 + U1*U1*cos(beta)*cos(beta) 
U2 = math.sqrt(U2)  
P2 = p2[np.argmin(abs(diff1))]
D2 = d2[np.argmin(abs(diff1))]    
T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',D2,fluidname) 
c2 = CP.CoolProp.PropsSI('A','P',P2,'T',T2,fluidname)  
M2 = U2/c2
h2 = CP.CoolProp.PropsSI('Hmass','P',P2,'Dmass',D2,fluidname) 
htotal2 = h2 + 0.5*U2*U2
print("ht2:", htotal2)
print("(ht2-ht1)/ht1:", (htotal2-htotal1)/htotal1)
s2 = CP.CoolProp.PropsSI('Smass','P',P2,'T',T2,fluidname) 
Pt2 = CP.CoolProp.PropsSI('P','Smass',s2,'Hmass',ht1,fluidname) 
Tt2 = CP.CoolProp.PropsSI('T','Smass',s2,'Hmass',ht1,fluidname) 
print("shock angle (degree):", beta*180/math.pi)
print("M2:", M2)
print("P2/P1:", P2/P1)
print("Pt2/Pt1:", Pt2/Pt1)
# print("T2/T1:", T2/T1)
# print("D2/D1:", D2/d1)

"""
3. find defelction angle theta
"""
n = 1000
theta1 = np.zeros(n) 
diff_theta = np.zeros(n) 
for i in range(n):
    theta1[i] = math.pi/2*i/n
    diff_theta[i] = U1*sin(beta)*tan(beta-theta1[i]) - u2[np.argmin(abs(diff1))]*tan(beta)
theta = theta1[np.argmin(abs(diff_theta))]
print("deflection angle (degree):", theta*180/math.pi)