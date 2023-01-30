#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 09:58:33 2022

@author: P.Yan

Enable new input and output pairs for CoolProp, including:
    (Z,P),(Z,T),(Z,Gamma),(Gamma,P),(Gamma,T)
"""
import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import pandas as pd


# give fluid name
fluidname = "MM"
# update fluid
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
Tmax =  CP.CoolProp.PropsSI("Tmax",fluidname)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
Pmax =  CP.CoolProp.PropsSI("pmax",fluidname)

# define input pairs (Z,P)
def TGfromZP(Z,P):
    """
    Parameters
    ----------
    P : double
        pressure.
    Z : double
        compressibility factor.
    Returns
    -------
    T : double
        temperature
    G : double
        fundamental derivative of gasdynamics Gamma.
    """
    print("fluid name is:", fluidname)
    print("input pairs Z,P[Pa]:" ,Z , P )
    # compute T for given Z,P
    Trange = np.linspace(Tc-100,Tmax,500)
    Trange = pd.Series(Trange)
    Z_error = np.zeros(Trange.size)
    for i in Trange.index:
        Z_error[i] = CP.CoolProp.PropsSI('Z','T',Trange[i],'P',P,fluidname) - Z
    Z_error = abs(Z_error)
    T = Trange[np.argmin(Z_error)]
    print("T[K] for given Z,P", T)
    #print("check: given P,T",P,T)
    print("Z is:", CP.CoolProp.PropsSI('Z','T',T,'P',P,fluidname))
    Gamma = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',T,'P',P,fluidname)
    print("Gamma is:", Gamma)
    return T,Gamma

# define input pairs (Z,P)
def PGfromZT(Z,T):
    print("fluid name is:", fluidname)
    print("input pairs Z,T[K]:" ,Z , T )
    # compute P for given Z,T
    Prange = np.linspace(Pc/10,Pc*3,1000)
    Prange = pd.Series(Prange)
    Z_error = np.zeros(Prange.size)
    for i in Prange.index:
        Z_error[i] = CP.CoolProp.PropsSI('Z','P',Prange[i],'T',T,fluidname) - Z
    Z_error = abs(Z_error)
    P = Prange[np.argmin(Z_error)]
    print("P[Pa] for given Z,T", P)
    print("Z is:", CP.CoolProp.PropsSI('Z','T',T,'P',P,fluidname))
    Gamma = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',T,'P',P,fluidname)
    print("Gamma is:", Gamma)
    return P,Gamma

# define input pairs (Z,Gamma)
def PTfromZG(Z,Gamma):
    #print("fluid name is:", fluidname)
    #print("input pairs Z,Gamma:" ,Z , Gamma )
    # compute P,T for given Z,Gamma
    Trange = np.linspace(Tc+30,Tmax/5,100)
    Trange = pd.Series(Trange)
    Prange = np.zeros(100)
    Gammat = np.zeros(100)
    # loop for P,T
    for i in Trange.index:
        Prange[i], Gammat[i] = PGfromZT(Z, Trange[i])
    Gammat = abs(Gammat - Gamma)
    P = Prange[np.argmin(Gammat)]
    T = Trange[np.argmin(Gammat)]
    #print("P[Pa],T[K]:", P, T)
    return P,T

# define input pairs (T,Gamma)
def ZPfromTG(T,Gamma):
    print("fluid name is:", fluidname)
    print("input pairs T,Gamma:" ,T , Gamma )
    # compute P for given t,Gamma
    Prange = np.linspace(Pc/5,Pc*5,500)
    Prange = pd.Series(Prange)
    G_error = np.zeros(Prange.size)
    for i in Prange.index:
        G_error[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','P',Prange[i],'T',T,fluidname) - Gamma
    G_error = abs(G_error)
    P = Prange[np.argmin(G_error)]
    print("P[Pa] for given T,Gamma", P)
    #print("Gamma is:", Gamma)
    Z = CP.CoolProp.PropsSI('Z','T',T,'P',P,fluidname)
    print("Z is:", Z)
    return Z,P

# define input pairs (P,Gamma)
def ZTfromPG(P,Gamma):
    print("fluid name is:", fluidname)
    print("input pairs P,Gamma:" ,P , Gamma )
    # compute P for given t,Gamma
    Trange = np.linspace(Tc-100,Tmax,500)
    Trange = pd.Series(Trange)
    G_error = np.zeros(Trange.size)
    for i in Trange.index:
        G_error[i] = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',Trange[i],'P',P,fluidname) - Gamma
    G_error = abs(G_error)
    T = Trange[np.argmin(G_error)]
    print("T[K] for given P,Gamma", T)
    Z = CP.CoolProp.PropsSI('Z','T',T,'P',P,fluidname)
    print("Z is:", Z)
    return Z,T

def ZGfromPT(P,T):
    Z = CP.CoolProp.PropsSI('Z','T',T,'P',P,fluidname)
    G = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',T,'P',P,fluidname)
    print("Z, Gamma: ", Z , G)