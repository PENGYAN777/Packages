#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 11:18:57 2022

@author: P.Yan

main function of CoolProp extensions
"""
import os

#os.system('python PT_diagram.py')
os.system('python plotContour.py')
# new input output pairs
from newIOpairs import TGfromZP, PGfromZT, PTfromZG, ZPfromTG, ZTfromPG, ZGfromPT

T1,Gamma1 = TGfromZP(0.9, 1e6)

#P2,Gamma2 = PGfromZT(0.9,500)

#P3, T3 =    PTfromZG(0.9, 0.9)

#Z4, P4  = ZPfromTG(500, 0.9)

#Z5, T5  = ZTfromPG(8e6,0.8)

# check consistence of P,T, Z,Gamma
ZGfromPT(1e6,T1)