#!/usr/bin/python

# This function plots the the triangulated surface obaine from the simulations
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from matplotlib import cm

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import



r0 = np.loadtxt("init_coords.dat", delimiter="\t") # Coordinates minim
rf = np.loadtxt("minim_coords.dat", delimiter="\t") # Coordinates minim
tri = np.loadtxt("mesh.dat", delimiter="\t") # Triangles

x0=r0[:,0]
y0=r0[:,1]
z0=r0[:,2]

xm=rf[:,0]
ym=rf[:,1]
zm=rf[:,2]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x0, y0, z0, triangles=tri, color=(0.5,0.5,0.5,0.5), edgecolor=(0.7,0.7,0.7, 1.0))
ax.plot_trisurf(xm, ym, zm, triangles=tri, color=(0.8,0.0,0.0,0.2), edgecolor=(0.5,0.0,0.0, 0.3))


plt.axis('off')
plt.show()
#plt.savefig("mesh.jpg", dpi=600)
plt.close()
