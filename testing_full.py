# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 23:39:27 2023

@author: Home-Pc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates
from lib import B_comp, grid
from field import dipolebetter

latitudes , longitudes = np.loadtxt('lat-lon.txt')
hs = np.pi/180
r = 696340
base_grid = grid(r, latitudes, longitudes, hs)

B_map = grid(r, latitudes, longitudes, hs)
m = np.asarray([0, 0, 10e20])

for latlon in B_map.latlon:
    r1 = coordinates(696340, latlon[0], latlon[1])
    B_map.set_value(dipolebetter(r1, m, returnBl=True), latlon[0], latlon[1])

r_high = 800000
high_grid = grid(r_high, latitudes, longitudes, hs)
B_map_high = grid(r_high, latitudes, longitudes, hs)
B_map_comp = grid(r_high, latitudes, longitudes, hs)

for latlon in B_map_high.latlon:
    r1 = coordinates(r_high, latlon[0], latlon[1], latlon=True)
    B_c = np.asarray(B_comp(r1, base_grid, B_map))
    B_map_comp.set_value(B_c, latlon[0], latlon[1], vector=True)
    B_d = np.asarray(dipolebetter(r1, m, returnBl=False, returnxyz=True))
    B_map_high.set_value(B_d, latlon[0], latlon[1], vector=True)
    print(B_c, B_d)

n = B_map_comp.num
xx = np.empty(n)
yy = np.empty(n)
uu = np.empty(n)
vv = np.empty(n)
cc = np.empty(n)
for i, (cs, val) in enumerate(zip(B_map_comp.coors_set, B_map_comp.valuesvector)):
    x, y = cs.project()
    u, z, v = val
    if abs(u) > 10 or abs(v) > 10:
        u, v = [0, 0]
    xx[i], yy[i], uu[i], vv[i] = x, y, u, v
    cc = math.sqrt(u**2 + v**2)

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))

n1 = B_map_high.num
xx1 = np.empty(n1)
yy1 = np.empty(n1)
uu1 = np.empty(n1)
vv1 = np.empty(n1)
cc1 = np.empty(n1)
for i, (cs, val) in enumerate(zip(B_map_high.coors_set, B_map_high.valuesvector)):
    x, y = cs.project()
    u, z, v = val
    xx1[i], yy1[i], uu1[i], vv1[i] = x, y, u, v
    cc1 = math.sqrt(u**2 + v**2)
    
ax2.quiver(xx1, yy1, uu1, vv1, np.arctan2(vv1, uu1))
plt.show()

