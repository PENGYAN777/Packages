#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:59:53 2023

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
"""
read csv file
"""
data= pd.read_csv("data.csv", ",", skiprows=0)


"""
plot 
"""
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
ax2 = axes.twinx()
axes.plot(data.iloc[:,2] , data.iloc[:,3] , 'ko', lw=lwh, label="Msup")
ax2.plot(data.iloc[:,2] , data.iloc[:,4] , 'k*', lw=lwh, label="${(P/P_t)}_{sup}$")

# ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax2.legend(loc=5 , prop={'size': 10}) 
ax2.set_ylabel('${(P/P_t)}_{sup}$',fontsize=12) 

# axes.set_ylim([1.9, 2.0])
axes.set_xlabel('$Z$',fontsize=12)
axes.set_ylabel('$M_{sup}$',fontsize=12) 
# axes.set_title('Mach number vs $Z$',fontsize=14)
axes.legend(loc=6 , prop={'size': 10}) #
fig1.savefig("files/M.pdf")

