# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 23:39:27 2023

@author: Home-Pc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import mplstereonet
from coordinates import coordinates
from lib import B_comp, grid
from field import dipolebetter


latitudes, longitudes = np.loadtxt('lat-lon.txt')
hs = np.pi / 90
r = 800000

base_grid = grid(r, latitudes, longitudes, hs)
B_map = grid(r, latitudes, longitudes, hs)
m = np.asarray([0, 0, 10e6])

for (lat, lon) in B_map.latlon:
    r1 = coordinates(r, lat, lon, latlon=True)
    B_map.set_value(dipolebetter(r1, m, returnBl=True), lat, lon)

r_high = 900000
high_grid = grid(r_high, latitudes, longitudes, hs)
B_map_high = grid(r_high, latitudes, longitudes, hs)
B_map_comp = grid(r_high, latitudes, longitudes, hs)

base_grid = grid(r, latitudes, longitudes, hs)

for (lat, lon) in B_map_high.latlon:
    r1 = coordinates(r_high, lat, lon, latlon=True)
    B_c = B_comp(r1, base_grid, B_map)

    B_map_comp.set_value(B_c, lat, lon, vector=True)

    B_d = np.asarray(dipolebetter(r1, m, returnBl=False, returnxyz=True))

    B_map_high.set_value(B_d, lat, lon, vector=True)
    print(lat, lon)
    print(f"[{B_c[0]:.2}, {B_c[1]:.2}, {B_c[2]:.2}]",
          f'[{B_d[0]:.2}, {B_d[1]:.2}, {B_d[2]:.2}]')

n = B_map_comp.num
xx = np.empty(n)
yy = np.empty(n)
uu = np.empty(n)
vv = np.empty(n)
cc = np.empty(n)
sq = np.empty(n)
for i, (cs, val) in enumerate(zip(B_map_comp.coors_set,
                                  B_map_comp.valuesvector)):
    x, y = cs.project()
    u, v, z = val
    if abs(u)>1e-4 or abs(v)>1e-4:
        u, v = 0, 0
    xx[i], yy[i], uu[i], vv[i] = x, y, u, v
    cc[i] = math.sqrt(u**2 + v**2)

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))

n1 = B_map_high.num
xx1 = np.empty(n1)
yy1 = np.empty(n1)
uu1 = np.empty(n1)
vv1 = np.empty(n1)
cc1 = np.empty(n1)
for i, (cs, val) in enumerate(zip(B_map_high.coors_set,
                                  B_map_high.valuesvector)):
    x, y = cs.project()
    u, v, z = val
    xx1[i], yy1[i], uu1[i], vv1[i] = x, y, u, v
    cc1[i] = math.sqrt(u**2 + v**2)

ax2.quiver(xx1, yy1, uu1, vv1, np.arctan2(vv1, uu1))

plt.show()
