#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 11:11:37 2022
@author: yan
Installing CoolProp
"""

import sys
import subprocess

# download CoolProp
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 
'CoolProp'])


