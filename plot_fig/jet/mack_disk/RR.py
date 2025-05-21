import math
import numpy as np
from numpy import sin, cos, pi, linspace
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Circle


"""
1. plot the nozzle wall and centerline
"""

lwh = 1
# Define the wall points

x = [0, 0.1, 0.1]
y = [0.2, 0.2, 0.4]

# Plot size settings
fig1 = plt.figure(dpi=300)
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7])  # size of figure

# Create the plot
axes.plot(x, y, 'k-',linewidth = lwh)

# # add  hash line along the wall

for x in [i * 0.02 for i in range(5)]:
    axes.plot([x, x], [0.2, 0.21], 'k-', linewidth = lwh/2)

for i in range(8):
    y = 0.22 + i * 0.02
    axes.plot([0.09, 0.1], [y, y], 'k-', linewidth = lwh/2)



# plot centerline 
axes.plot([0,0.7], [0,0], 'k-.', linewidth = lwh)

"""
2. plot the flow features
"""

# constant pressure boundary
# Define center and point on the circle
# Given data
center = np.array([0.3, -0.2])
point1 = np.array([0.1, 0.2])
point2 = np.array([0.5, 0.2])

# Calculate radius
radius = np.linalg.norm(point1 - center)

# Calculate angles of points relative to center
angle1 = np.arctan2(point1[1] - center[1], point1[0] - center[0])
angle2 = np.arctan2(point2[1] - center[1], point2[0] - center[0])


angles = np.linspace(angle1, angle2, 100)

# Calculate arc points
arc_x = center[0] + radius * np.cos(angles)
arc_y = center[1] + radius * np.sin(angles)

axes.plot(arc_x, arc_y, 'b--', linewidth=lwh, label='Constant pressure boundary')
# expansion fans

axes.plot([0.1,0.16], [0.2,0.0], 'k--', linewidth = lwh, label='Expansion fans')
axes.plot([0.1,0.18], [0.2,0.0], 'k--', linewidth = lwh)
axes.plot([0.1,0.20], [0.2,0.0], 'k--', linewidth = lwh)

# relected expansion fans
axes.plot([0.16,0.20], [0.0,0.240], 'k--', linewidth = lwh)
axes.plot([0.18,0.24], [0.0,0.244], 'k--', linewidth = lwh)
axes.plot([0.20,0.30], [0.0,0.246], 'k--', linewidth = lwh)


# compression fans

axes.plot([0.20,0.40], [0.240,0.10], 'g--', linewidth = lwh,label='Compression fans')
axes.plot([0.24,0.40], [0.244,0.10], 'g--', linewidth = lwh)
axes.plot([0.30,0.40], [0.246,0.10], 'g--', linewidth = lwh)



# # interception and reflected shock, Mach disk

axes.plot([0.40,0.50], [0.10,0.00], 'r', linewidth = lwh, label='Interception/Reflected shock')
axes.plot([0.50,0.60], [0.00,0.10], 'r', linewidth = lwh)




# # Labels 

# Define circle parameters

center = (0.03, 0.10)
radius = 0.015  # Adjust as needed

# Create the circle patch
circle = Circle(center, radius, edgecolor='k', facecolor='white', linewidth=lwh)
axes.add_patch(circle)

plt.text(0.03, 0.095, '1',ha= 'center', size =6)
plt.text(0.09, 0.095, '$P_e>P_a$',ha= 'center', size =8)

######################################################################################

center = (0.16, 0.16)
radius = 0.015  # Adjust as needed

# Create the circle patch
circle = Circle(center, radius, edgecolor='k', facecolor='white', linewidth=lwh)
axes.add_patch(circle)

plt.text(0.16, 0.155, '2',ha= 'center', size =6)
plt.text(0.16, 0.195, '$P_a$',ha= 'center', size =8)

######################################################################################
center = (0.30, 0.10)
radius = 0.015  # Adjust as needed

# Create the circle patch
circle = Circle(center, radius, edgecolor='k', facecolor='white', linewidth=lwh)
axes.add_patch(circle)

plt.text(0.30, 0.095, '3',ha= 'center', size =6)
plt.text(0.30, 0.055, '$P_1<P_a$',ha= 'center', size =8)


######################################################################################

center = (0.50, 0.10)
radius = 0.015  # Adjust as needed

# Create the circle patch
circle = Circle(center, radius, edgecolor='k', facecolor='white', linewidth=lwh)
axes.add_patch(circle)

plt.text(0.50, 0.095, '4',ha= 'center', size =6)
plt.text(0.50, 0.155, '$P_a$',ha= 'center', size =8)


######################################################################################

center = (0.25, 0.30)
radius = 0.015  # Adjust as needed

# Create the circle patch
circle = Circle(center, radius, edgecolor='k', facecolor='white', linewidth=lwh)
axes.add_patch(circle)

plt.text(0.25, 0.295, '5',ha= 'center', size =6)
plt.text(0.30, 0.295, '$P_a$',ha= 'center', size =8)




# Set aspect ratio to 'equal' and box
axes.set_aspect('equal', 'box')

# Turn off the axes
plt.axis('off')

axes.legend(loc=0, prop={'size': 7})

# Save the figure as a .eps file
fig1.savefig("geo_RR.png")
