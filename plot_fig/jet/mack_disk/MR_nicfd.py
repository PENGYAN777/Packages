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

axes.plot([0.40,0.45], [0.10,0.05], 'r', linewidth = lwh, label='Interception/Reflected shock')
axes.plot([0.45,0.45], [0.05,0.00], 'k', linewidth = lwh, label='Mach Disk')
axes.plot([0.45,0.50], [0.05,0.10], 'r', linewidth = lwh)

# slip line

axes.plot([0.45,0.55], [0.05,0.05], 'b', linewidth = lwh, label='Slip Line')



# # Labels 

plt.text(0.30, 0.05, 'Ideal',ha= 'center', size =10)


""""
3. plot the mirrow part
"""

# Mirror nozzle wall and hash lines (y -> -y)
x_wall = np.array([0, 0.1, 0.1])
y_wall = np.array([0.2, 0.2, 0.4])

axes.plot(x_wall, -y_wall, 'k-', linewidth=lwh)  # mirrored wall

# mirrored hash lines along the wall
for xh in [i * 0.02 for i in range(5)]:
    axes.plot([xh, xh], [-0.2, -0.21], 'k-', linewidth=lwh / 2)

for i in range(8):
    y_hash = 0.22 + i * 0.02
    axes.plot([0.09, 0.1], [-y_hash, -y_hash], 'k-', linewidth=lwh / 2)

# mirrored constant pressure boundary arc
axes.plot(arc_x, -arc_y, 'b--', linewidth=lwh)

# mirrored expansion fans
axes.plot([0.1, 0.14], [-0.2, -0.0], 'k--', linewidth=lwh)
axes.plot([0.1, 0.16], [-0.2, -0.0], 'k--', linewidth=lwh)
axes.plot([0.1, 0.18], [-0.2, -0.0], 'k--', linewidth=lwh)

# mirrored reflected expansion fans
axes.plot([0.14, 0.18], [-0.0, -0.23], 'k--', linewidth=lwh)
axes.plot([0.16, 0.22], [-0.0, -0.244], 'k--', linewidth=lwh)
axes.plot([0.18, 0.28], [-0.0, -0.245], 'k--', linewidth=lwh)

# # mirrored compression fans
axes.plot([0.19, 0.35], [-0.240, -0.10], 'g--', linewidth=lwh)
axes.plot([0.23, 0.35], [-0.244, -0.10], 'g--', linewidth=lwh)
axes.plot([0.29, 0.35], [-0.246, -0.10], 'g--', linewidth=lwh)

# # mirrored interception/reflected shock and Mach disk
axes.plot([0.35, 0.40], [-0.10, -0.05], 'r', linewidth=lwh)
axes.plot([0.40, 0.40], [-0.05, -0.00], 'k', linewidth=lwh)
axes.plot([0.40, 0.45], [-0.05, -0.10], 'r', linewidth=lwh)

# # mirrored slip line
axes.plot([0.40, 0.50], [-0.05, -0.05], 'b', linewidth=lwh)

# Optional: Add symmetric labels, adjust position as needed

plt.text(0.30, -0.07, 'Nonideal', ha='center', size=10)



# Set aspect ratio to 'equal' and box
# axes.set_aspect('equal', 'box')

# Turn off the axes
plt.axis('off')

axes.legend(loc=0, prop={'size': 5})

# Save the figure as a .eps file
fig1.savefig("geo_MR_nicfd.png")
