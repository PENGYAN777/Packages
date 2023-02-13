#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute normal shock relatiomnships for CoolProp EOS
input: upstream P,T,M.

"""

from scipy.optimize import fsolve
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import CoolProp as CP


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
1. input pre-shock conditions
"""


P1 = Pc*0.4 # pre-shock pressure
T1 = Tc*1.1 # pre-shcok temperature
d1 = CP.CoolProp.PropsSI('Dmass','P',P1,'T',T1,fluidname) # pre-shock density
Z1 = CP.CoolProp.PropsSI('Z','P',P1,'T',T1,fluidname) 
G1 = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',P1,'T',T1,fluidname) 
print("Z,Gamma:", Z1,G1)
s1 = CP.CoolProp.PropsSI('Smass','P',P1,'T',T1,fluidname) 
cv = CP.CoolProp.PropsSI('Cvmass','P',P1,'T',T1,fluidname) 
cp = CP.CoolProp.PropsSI('Cpmass','P',P1,'T',T1,fluidname) 
print("gamma:", cp/cv)
M1 = 2
c1 = CP.CoolProp.PropsSI('A','P',P1,'T',T1,fluidname)
u1 = M1*c1 
h1 = CP.CoolProp.PropsSI('Hmass','P',P1,'T',T1,fluidname) 
ht1 = h1 + 0.5*u1*u1
print("ht:", ht1)
# compute total pressure/temperature
Pt1 = CP.CoolProp.PropsSI('P','Smass',s1,'Hmass',ht1,fluidname) 
Tt1 = CP.CoolProp.PropsSI('T','Smass',s1,'Hmass',ht1,fluidname) 
htotal1 = CP.CoolProp.PropsSI('Hmass','P',Pt1,'T',Tt1,fluidname) 
print("ht1:", htotal1)

"""
2. compute post-shock properites
"""

p2 = np.linspace(P1*2,P1*4.1,1000) # post-shock pressure 
p2 = pd.Series(p2)
u2 = np.zeros(p2.size) 
d2 = np.zeros(p2.size) 
diff = np.zeros(p2.size) 
 
for i in p2.index:
    if abs(p2[i]-Pc)<0.01*Pc:
        p2[i] = 0.99*Pc
    u2[i] = (P1+d1*u1*u1-p2[i])/d1/u1
    d2[i] =  d1*u1/u2[i]
    diff[i] = ht1 - 0.5*u2[i]*u2[i] - CP.CoolProp.PropsSI('Hmass','P',p2[i],'Dmass',d2[i],fluidname) 
print("min diff:", diff[np.argmin(abs(diff))])
P2 = p2[np.argmin(abs(diff))]
D2 = d2[np.argmin(abs(diff))]    
T2 = CP.CoolProp.PropsSI('T','P',P2,'Dmass',D2,fluidname) 
c2 = CP.CoolProp.PropsSI('A','P',P2,'T',T2,fluidname) 
U2 = u2[np.argmin(abs(diff))]    
M2 = U2/c2
h2 = CP.CoolProp.PropsSI('Hmass','P',P2,'Dmass',D2,fluidname) 
htotal2 = h2 + 0.5*U2*U2
print("ht2:", htotal2)
print("M2:", M2)
