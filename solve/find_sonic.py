#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute Prandtl-Meyer relatiomnships for CoolProp EOS
input: upstream P,Z,M and defletion angle theta.

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
0. fluid property
"""
fluidname = "MM"
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
1. input total conditions
"""


pt = 1.4e6 # total pressure
zt = 0.6

tt,gt = TGfromZP(zt,pt)
dt = CP.CoolProp.PropsSI('Dmass','P',pt,'T',tt,fluidname) 
s = CP.CoolProp.PropsSI('Smass','P',pt,'T',tt,fluidname) 
ht = CP.CoolProp.PropsSI('Hmass','P',pt,'T',tt,fluidname) 
"""
2. compute isentropic relationship
"""

p = np.linspace(pt*0.1,pt,1000) # P<Pc

p = pd.Series(p)
Z = np.zeros(p.size) # P/rho RT
h = np.zeros(p.size) # enthalpy
u = np.zeros(p.size) # velocity
v = np.zeros(p.size) # specific volume
c = np.zeros(p.size) # sound speed
m = np.zeros(p.size) # Mach number
d = np.zeros(p.size) # density
t = np.zeros(p.size) # temperature
g = np.zeros(p.size) # Gamma

for i in p.index:
    if abs(p[i]-Pc)<0.05*Pc:
        p[i] = 0.95*Pc
    d[i] = CP.CoolProp.PropsSI('Dmass','P',p[i],'Smass',s,fluidname) 
    v[i] = 1/d[i]
    t[i] = CP.CoolProp.PropsSI('T','P',p[i],'Smass',s,fluidname) 
    Z[i] = CP.CoolProp.PropsSI('Z','P',p[i],'Smass',s,fluidname) 
    g[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',p[i],'Smass',s,fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[i],fluidname)
    u[i] = math.sqrt(abs(2*(ht-h[i])))
    c[i] = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s,fluidname) 
    m[i] = u[i]/c[i]

"""
3. find sonic condition, assume A* = 1
"""
print("index for sonic condition:",np.argmin(abs(m-1)))
dstar = d[np.argmin(abs(m-1))]
vstar = 1/dstar
Tstar = t[np.argmin(abs(m-1))]
Pstar = p[np.argmin(abs(m-1))]



