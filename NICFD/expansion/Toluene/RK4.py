#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 17:10:25 2023

@author: yan
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import CoolProp as CP

fluidname = "Toluene"

def rk4(v0,vn,t0,m0,nu0,n):
    dv = (vn-v0)/n          # estep Delta v
    """
    initialization
    """
    v_arr = np.zeros(n + 1)   # create an array of zeros for spefific volume
    t_arr = np.zeros(n +1)    # create an array of zeros for Temperature
    m_arr = np.zeros(n + 1)   # create an array of zeros for Mach number
    nu_arr = np.zeros(n + 1)   # create an array of zeros for Prandtl Mayer function nu
    
    v_arr[0] = v0              # add initial value to array
    t_arr[0] = t0              # add initial value to array
    m_arr[0] = m0              # add initial value to array
    nu_arr[0] = nu0              # add initial value to array

    # Euler's method
    for i in range (1, n + 1):  
       v = v_arr[i-1]
       t = t_arr[i-1]
       m = m_arr[i-1]
       nu = nu_arr[i-1]
       g = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',t,'Dmass',1/v,fluidname ) 
       p = CP.CoolProp.PropsSI('P','T',t,'Dmass',1/v,fluidname ) 
       dnudv = math.sqrt(m*m-1)/v/m/m
       dmdv = -m/v*(1 - g - 1/m/m)
       dtdv = -t/CP.CoolProp.PropsSI('Cvmass','T',t,'Dmass',1/v,"Toluene")* CP.CoolProp.PropsSI('d(P)/d(T)|Dmass','P',p,'T',t,fluidname )
       v_arr[i] = v + dv
       t_arr[i] = t + dv*dtdv
       m_arr[i] = m + dv*dmdv
       nu_arr[i] = nu + dv*dnudv
       
    return v_arr,t_arr,m_arr,nu_arr
