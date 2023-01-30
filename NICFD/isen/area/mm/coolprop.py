#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

compute isentropic relatiomnships with constat area 
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

"""
1. input total conditions
"""

pt = 1.3e6 # total pressure
print("total pressure[Pa]:", pt)
zt = 0.6 # total compressibility factor
print("total compressibility factor", zt)
tt,gt = TGfromZP(zt, pt)
# tt = 500 #total temprature
print("total temperature[K]:", tt)
print("total Gamma:", gt)
dt = CP.CoolProp.PropsSI('Dmass','P',pt,'T',tt,fluidname) 
ct = CP.CoolProp.PropsSI('Dmass','P',pt,'T',tt,fluidname) 
s = CP.CoolProp.PropsSI('Smass','T',tt,'P',pt,fluidname)
print("entropy:", s)
ht = CP.CoolProp.PropsSI('Hmass','T',tt,'P',pt,fluidname)
print("total enthalpy:", ht)

"""
2. read area ratio
"""
wall = pd.read_csv("nozzle.csv", " ", skiprows=0)
A = wall.iloc[:,2]
Astar = A[np.argmin(abs(A))]
Ar = A/Astar # A/A^*, i.e., area ratio
x = A = wall.iloc[:,1]
    
"""
3. find conic conditions
"""    
# note that if we use more points, i.e., from 100 to 1000, the curve wil be more smooth.
p = np.linspace(pt*0.01,pt,500) # static pressure
p = pd.Series(p)
m = np.zeros(p.size) # Mach number
for i in p.index:
    c = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s,fluidname)
    h = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[i],fluidname)
    v = math.sqrt(abs(2*(ht-h)))
    m[i] = v/c
pstar = p[np.argmin(abs(m-1))]
tstar = CP.CoolProp.PropsSI('T','P',pstar,'Smass',s,fluidname)
dstar = CP.CoolProp.PropsSI('Dmass','P',pstar,'Smass',s,fluidname)
cstar = CP.CoolProp.PropsSI('A','P',pstar,'Smass',s,fluidname)
hstar = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',pstar,fluidname)
vstar = math.sqrt(abs(2*(ht-hstar)))

"""
4. compute thermodynamics properties at differenet area section.
! remember that at a given A/A^*, we can have subsonic and supersonic flow.
"""
Psub = np.zeros(Ar.size) # static pressure for subsonic flow
Tsub = np.zeros(Ar.size) # static temperature
Dsub = np.zeros(Ar.size) # static density
Csub = np.zeros(Ar.size) # static speed of sound
Msub = np.zeros(Ar.size) # Mach number
Psup = np.zeros(Ar.size) # static pressure for supersonic flow
Tsup = np.zeros(Ar.size) # static temperature
Dsup = np.zeros(Ar.size) # static density
Csup = np.zeros(Ar.size) # static speed of sound
Msup = np.zeros(Ar.size) # Mach number
for i in Ar.index:
    print("i= ",i)
    """
    for subsonic flow
    """
    r = np.zeros(p.size) # ratio: (rho*A*M*c)/()^* 
    for j in p.index:
        dsub = CP.CoolProp.PropsSI('Dmass','P',p[j],'Smass',s,fluidname)
        hsub = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[j],fluidname)
        vsub = math.sqrt(abs(2*(ht-hsub)))
        csub = CP.CoolProp.PropsSI('A','P',p[j],'Smass',s,fluidname)
        msub = vsub/csub
        if msub<1:
            r[j]= dsub*vsub*Ar[i]/dstar/vstar
    k = np.argmin(abs(r-1))
    Psub[i] =p[k]
    Tsub[i] = CP.CoolProp.PropsSI('T','P',p[k],'Smass',s,fluidname)
    Dsub[i] = CP.CoolProp.PropsSI('Dmass','P',p[k],'Smass',s,fluidname)
    Csub[i] = CP.CoolProp.PropsSI('A','P',p[k],'Smass',s,fluidname)
    h = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[k],fluidname)
    v = math.sqrt(abs(2*(ht-h)))
    Msub[i] = v/Csub[i]
    """
    for supersonic flow
    """
    r = np.zeros(p.size) # ratio: (rho*A*M*c)/()^* 
    for j in p.index:
        dsup = CP.CoolProp.PropsSI('Dmass','P',p[j],'Smass',s,fluidname)
        hsup = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[j],fluidname)
        vsup = math.sqrt(abs(2*(ht-hsup)))
        csup = CP.CoolProp.PropsSI('A','P',p[j],'Smass',s,fluidname)
        msup = vsup/csup
        if msup>1:
            r[j]= dsup*vsup*Ar[i]/dstar/vstar
    k = np.argmin(abs(r-1))
    Psup[i] =p[k]
    Tsup[i] = CP.CoolProp.PropsSI('T','P',p[k],'Smass',s,fluidname)
    Dsup[i] = CP.CoolProp.PropsSI('Dmass','P',p[k],'Smass',s,fluidname)
    Csup[i] = CP.CoolProp.PropsSI('A','P',p[k],'Smass',s,fluidname)
    h = CP.CoolProp.PropsSI('Hmass','Smass',s,'P',p[k],fluidname)
    v = math.sqrt(abs(2*(ht-h)))
    Msup[i] = v/Csup[i]
            
"""
5. compute by assuming flow from left to right
"""
P = np.zeros(Ar.size) # static pressure 
T = np.zeros(Ar.size) # static temperature
D = np.zeros(Ar.size) # static density
C = np.zeros(Ar.size) # static speed of sound
M = np.zeros(Ar.size) # Mach number
Z = np.zeros(Ar.size) # compressibility factor
G = np.zeros(Ar.size) # Gamma
for i in Ar.index:
    if i<np.argmin(Ar):
        P[i] = Psub[i]
        T[i] = Tsub[i]
        D[i] = Dsub[i]
        C[i] = Csub[i]
        M[i] = Msub[i]
        Z[i] = CP.CoolProp.PropsSI('Z','P',P[i],'Smass',s,fluidname) 
        G[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',P[i],'Smass',s,fluidname)
    else:
        P[i] = Psup[i]
        T[i] = Tsup[i]
        D[i] = Dsup[i]
        C[i] = Csup[i]
        M[i] = Msup[i]
        Z[i] = CP.CoolProp.PropsSI('Z','P',P[i],'Smass',s,fluidname) 
        G[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',P[i],'Smass',s,fluidname)
        
"""
6. write data into file
"""
P = P/pt
D = D/dt
T = T/tt
C = C/ct

pd.DataFrame(P).to_csv('z6.csv', index_label = "Index", header  = ['pressure']) 
data = pd.read_csv("z6.csv", ",")
# append new columns
d =pd.DataFrame({'density': D, 'temperature': T, 'soundspeed': C,'Mach': M, 'Z': Z,'g': G,'A':Ar, 'x':x})
newData = pd.concat([data, d], join = 'outer', axis = 1)
# save newData in csv file
# newData.to_csv("m4sh.csv")
newData.to_csv("z6.csv")





