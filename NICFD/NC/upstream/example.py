#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:01:16 2024

@author: yan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Given points
x = np.array([1,2,2,1])
y = np.array([2,2,4,4])
f = np.array([3,5,7,9])

# Generate a meshgrid for plotting
x_grid = np.linspace(min(x), max(x), 10)
y_grid = np.linspace(min(y), max(y), 10)
X, Y = np.meshgrid(x_grid, y_grid)

# Interpolate f values over the meshgrid
Z = griddata((x, y), f, (X, Y), method='cubic')

# Plot contour
plt.contourf(X, Y, Z, cmap='viridis')
plt.colorbar(label='f(x, y)')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour Plot of f(x, y)')
plt.scatter(x, y, c='red', label='Data Points')  # Plot data points for reference
plt.legend()
plt.show()
