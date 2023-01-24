#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 11:18:57 2022

@author: P.Yan

main function of CoolProp extensions
"""
import os
import CoolProp as CP
#os.system('python PT_diagram.py')
#os.system('python plotContour.py')
# new input output pairs
from newIOpairs import TGfromZP, PGfromZT, PTfromZG, ZPfromTG, ZTfromPG, ZGfromPT

#T1,Gamma1 = TGfromZP(0.9, 1e6)

#P2,Gamma2 = PGfromZT(0.9,500)

#P3, T3 =    PTfromZG(0.9, 0.9)

#Z4, P4  = ZPfromTG(500, 0.9)

#Z5, T5  = ZTfromPG(8e6,0.8)

# check consistence of P,T, Z,Gamma
PGfromZT(0.5,160)

# compute active degree of freedom
print("------------compute N-----------")
fluidname = "MM"
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
w =  CP.CoolProp.PropsSI("acentric",fluidname)
print("caentric factor:", w)
cv =CP.CoolProp.PropsSI('Cvmass','T',Tc,'P',Pc/10,fluidname)
print("isochoric specific heat capacity : 	J/kg/K", cv)
N = 2*cv/Rs
print("active degree of freddom:", N)
Z = CP.CoolProp.PropsSI('Z','T',Tc,'P',Pc/10,fluidname)
print("Z:",Z)

Z = CP.CoolProp.PropsSI('Z','T',542.13,'P',9.04e5,"MDM")
print("Z:",Z)
G = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',542.13,'P',9.04e5,"MDM")
print("G:",G)

