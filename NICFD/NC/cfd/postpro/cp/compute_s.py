#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:51:54 2024

@author: yan
"""

import numpy as np
import CoolProp as CP
import pandas as pd

fluidname = "HEOS::MD4M"
data = pd.read_csv("m6.csv", ",")
P = data.iloc[:,6] 
T = data.iloc[:,8] 
D = data.iloc[:,0] 
print("size", P.index)

G = np.zeros(P.size)
Z = np.zeros(P.size)
# ht = np.zeros(P.size)
s = np.zeros(P.size)
# c = np.zeros(P.size)
for i in P.index:
    G[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics',
                                'T',T[i],'P',P[i],fluidname)
    Z[i] =  CP.CoolProp.PropsSI('Z','T',T[i],'P',P[i],fluidname)
    s =  CP.CoolProp.PropsSI('S','T',T[i],'P',P[i],fluidname)
    # c[i] =  CP.CoolProp.PropsSI('SPEED_OF_SOUND','T',T[i],'P',P[i],fluidname)
    # u = c[i]*data.iloc[i,1] 
    # ht[i] = h + 0.5*u*u
    
# append new columns
# shG =pd.DataFrame({'G':G, })
shG =pd.DataFrame({'G':G, 'Z':Z, 's':s,})
newData = pd.concat([data, shG], join = 'outer', axis = 1)
# save newData in csv file
# newData.to_csv("m4sh.csv")
newData.to_csv("m6new.csv")