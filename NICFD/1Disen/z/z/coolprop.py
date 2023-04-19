#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:35:29 2022

@author: yan

investigate effect of Z on quasi 1D isentropic flow

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
1. input 
"""
Ar = 2 # A/A^*
g = 1.4 # cp/cv
Z = np.arange(0.50, 1.00, 0.1).tolist()
Msup = np.zeros(len(Z)) 
Psup = np.zeros(len(Z)) 
Msub = np.zeros(len(Z)) 
Psub = np.zeros(len(Z)) 

for i in range(len(Z)):
    z = Z[i]
    print('-----------------------Zt----------------------:',z)
    """
    2. solve Mach number and P/Pt
    """
    
    """
    2.1 supersonic  
    """
    print('-----------------------Supersonic--------------------')
    n = 1000 # number of test points
    M = np.linspace(1,3,n)
    diff = np.zeros(n)
    for j in range(n):
        m = M[j]
        diff[j] = Ar*m - ((2+(g-1)*z)/(2+(g-1)*z*m*m))**((g+1)/(2-2*g)/z)
        diff[j] = abs(diff[j])
    
    Msup[i] = M[np.argmin(abs(diff))]
    print('Msup:',Msup[i])
    
    Psup[i] = (1+(g-1)/2*z*Msup[i]*Msup[i])**(g/(g-1)/z) #P/Pt
    Psup[i] = 1/Psup[i]
    print('P/Pt supersonic:', Psup[i])
    
    """
    2.2 subsonic  
    """
    print('-----------------------Supersonic--------------------')
    M = np.linspace(0.01,1,n)
    diff = np.zeros(n)
    for j in range(n):
        m = M[j]
        diff[j] = Ar*m - ((2+(g-1)*z)/(2+(g-1)*z*m*m))**((g+1)/(2-2*g)/z)
        diff[j] = abs(diff[j])
    
    Msub[i] = M[np.argmin(abs(diff))]
    print('Msub:',Msub[i])
    
    Psub[i] = (1+(g-1)/2*z*Msub[i]*Msub[i])**(g/(g-1)/z) #P/Pt
    Psub[i] = 1/Psub[i]
    print('P/Pt subsonic:', Psub[i])


"""
3. Write csv file 
"""
data = 'data.csv'
pd.DataFrame(Z).to_csv(data, index_label = "Index", header  = ['Z']) 
result = pd.read_csv(data, ",")
# append new columns
D =pd.DataFrame({'Msup': Msup, 'Psup': Psup, 'Msub': Msub, 'Psub': Psub, })
newData = pd.concat([result, D], join = 'outer', axis = 1)
# save newData in csv file
# newData.to_csv("m4sh.csv")
newData.to_csv(data)