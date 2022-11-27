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

fluidname = "nitrogen"
data = pd.read_csv("pr.csv", " ")
# get P,T from 2nd and 3rd column 
P = data.iloc[:,1] * 2e5
T = data.iloc[:,2] * 300
print("size", P.index)
Z = np.zeros(P.size)
Gamma = np.zeros(P.size)
for i in P.index:
    Z[i] =  CP.CoolProp.PropsSI('Z','T',T[i],'P',P[i],fluidname)
    Gamma[i] =  CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',T[i],'P',P[i],fluidname)

# append new columns
ZG =pd.DataFrame({'Z': Z,'Gamma': Gamma})
newData = pd.concat([data, ZG], join = 'outer', axis = 1)

# save newData in csv file
newData.to_csv("prNew.csv")
