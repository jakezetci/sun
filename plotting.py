# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:52:18 2023

@author: cosbo
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from coordinates import coordinates
from lib import B_comp, grid
from field import B_dipole
import numpy as np
import math

latitudes , longitudes = np.loadtxt('lat-lon.txt')
hs = 90/15
r = 696340
base_grid = grid(r, latitudes, longitudes, hs)

B_map = grid(r, latitudes, longitudes, hs)

for latlon in B_map.latlon:
    r1 = coordinates(696340, latlon[0], latlon[1])
    B_map.set_value(B_dipole(r1, returnBl=True),
                    latlon[0], latlon[1])

r_high = 800000

B_map_high = grid(r_high, latitudes, longitudes, hs)
B_map_comp = grid(r_high, latitudes, longitudes, hs)

for latlon in B_map_high.latlon:
    r1 = coordinates(r_high, latlon[0], latlon[1], latlon=True)
    B_c = np.asarray(B_comp(r1, base_grid, B_map))
    B_map_comp.set_value(B_c, latlon[0], latlon[1], vector=True)

n = B_map_comp.num
xx = np.empty(n)
yy = np.empty(n)
uu = np.empty(n)
vv = np.empty(n)
cc = np.empty(n)
for i, (cs, val) in enumerate(zip(B_map_comp.coors_set, B_map_comp.valuesvector)):
    x, y = cs.project()
    u, z, v = val
    xx[i], yy[i], uu[i], vv[i] = x, y, u, v
    cc = math.sqrt(u**2 + v**2)

plt.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))
plt.show()
