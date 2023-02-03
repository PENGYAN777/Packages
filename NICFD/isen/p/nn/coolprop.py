#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute isentropic relatiomnships for PIG
input: total pressure and compressibiility factor, pt and zt

"""

from scipy.optimize import fsolve
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import CoolProp as CP
from newIOpairs import TGfromZP, PGfromZT, PTfromZG, ZPfromTG, ZTfromPG, ZGfromPT

"""
fluid property
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

"""
total conditions
"""

pt = 3e6 # total pressure
print("total pressure[Pa]:", pt)
zt = 0.75 # total compressibility factors
print("total compressibility factor", zt)
tt,gt = TGfromZP(zt, pt)
# tt = 500 #total temprature
print("total temperature[K]:", tt)
print("total Gamma:", gt)
dt = CP.CoolProp.PropsSI('Dmass','P|gas',pt,'T',tt,fluidname) 
ct = CP.CoolProp.PropsSI('Dmass','P|gas',pt,'T',tt,fluidname) 
s = CP.CoolProp.PropsSI('Smass','T',tt,'P',pt,fluidname)
print("entropy:", s)
ht = CP.CoolProp.PropsSI('Hmass','T',tt,'P',pt,fluidname)
print("total enthalpy:", ht)

"""
compute isentropic relationship
"""

p = np.linspace(pt*0.1,pt,100) # 
p = pd.Series(p)
Z = np.zeros(p.size)
h = np.zeros(p.size)
v = np.zeros(p.size)
c = np.zeros(p.size)
m = np.zeros(p.size)
d = np.zeros(p.size)
t = np.zeros(p.size)
g = np.zeros(p.size)
A = np.zeros(p.size)
for i in p.index:
    t[i] = CP.CoolProp.PropsSI('T','P|gas',p[i],'Smass',s,fluidname) 
    d[i] = CP.CoolProp.PropsSI('Dmass','P|gas',p[i],'T',t[i],fluidname) 
    Z[i] = CP.CoolProp.PropsSI('Z','P|gas',p[i],'Smass',s,fluidname) 
    g[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P|gas',p[i],'T',t[i],fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','Smass',s,'P|gas',p[i],fluidname)
    v[i] = math.sqrt(abs(2*(ht-h[i])))
    c[i] = CP.CoolProp.PropsSI('A','P|gas',p[i],'T',t[i],fluidname) 
    m[i] = v[i]/c[i]

"""
find sonic condition, assume A* = 1
"""
dstar = d[np.argmin(abs(m-1))]
vstar = v[np.argmin(abs(m-1))]
for i in p.index:
    A[i] = dstar*vstar/d[i]/v[i]


"""
write into csv file
"""    

p = p/pt
d = d/dt
t = t/tt
c = c/ct

pd.DataFrame(p).to_csv('z75.csv', index_label = "Index", header  = ['pressure']) 
data = pd.read_csv("z75.csv", ",")
# append new columns
d =pd.DataFrame({'density': d, 'temperature': t, 'soundspeed': c,'Mach': m, 'Z': Z,'g': g,'A':A})
newData = pd.concat([data, d], join = 'outer', axis = 1)
# save newData in csv file
# newData.to_csv("m4sh.csv")
newData.to_csv("z75.csv")

    
    
