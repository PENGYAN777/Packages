#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 19:10:17 2022

@author: P.Yan

read csv file with P,T, compute Z, Gamma and wirte in thee csv file
"""

import numpy as np
import CoolProp as CP
import pandas as pd


################# compute Gamma for van der wall model
fluidname = "HEOS::MD4M"
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
dc =  CP.CoolProp.PropsSI("rhocrit",fluidname)
R = CP.CoolProp.PropsSI('GAS_CONSTANT',fluidname)
print("Universal gas constant:", R)
MW = CP.CoolProp.PropsSI('M',fluidname)
# print("molar mass:", MW)
Rs = R/MW

p = 0.922*Pc
t = 0.9836*Tc
d = 0.4928*dc
# heat capaticty gamma cp/cv
g = 1.01

a = 27*Rs*Rs*Tc*Tc/64/Pc
b = 1/8*Rs*Tc/Pc
# cv for PIG
# cv = CP.CoolProp.PropsSI('CVMASS','T', t, 'P',p*0.1,  fluidname) 
cv = Rs/(g-1)
Z = CP.CoolProp.PropsSI('Z','T', t, 'P',p*0.1,  fluidname) 
N = 2*cv/Rs
# N = 120
v = 1/d

#reference value 0
v0 = v*0.5
t0 = t*0.5
c = (1+2/N)*Rs*t*(v/(v-b))**2-2*a/v
s = Rs*np.log((v-b)/(v0-b)) + (Rs)/(g-1)*np.log(t/t0)

# new density with same entropy
v1 = 1.18*v0
d1 = 1/v1
# find t1 = t1(s,v1)
t2 = np.linspace(t*0.9,t*1.1,500) 
t2 = pd.Series(t2)
s1 = np.zeros(t2.size) 
for i in t2.index:
    s1[i] = Rs*np.log((v1-b)/(v0-b)) + Rs/(g-1)*np.log(t2[i]/t0)
print("index:",np.argmin(abs(s1-s)))    
t1 = t2[np.argmin(abs(s1-s))]
c1 = (1+2/N)*Rs*t1*(v1/(v1-b))**2-2*a/v1
Gamma = 1 + d/c*(c1-c)/(d1-d)
print("Gamma:", Gamma)
