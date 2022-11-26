#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 11:36:47 2022

@author: yan

PLot contour of Z and Gamma for given fluid in (P,T) diagram
"""

from CoolProp.CoolProp import PropsSI
print(PropsSI('D', 'P', 1E5, 'T', 300, 'nitrogen'))
