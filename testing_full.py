# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 23:39:27 2023

@author: Home-Pc
"""

import numpy as np
from coordinates import coordinates
from lib import B_comp, grid
from field import B_dipole

latitudes , longitudes = np.loadtxt('lat-lon.txt')
hs = np.pi/180
r = 696340
base_grid = grid(r, latitudes, longitudes, hs)

B_map = grid(r, latitudes, longitudes, hs)


for latlon in B_map.latlon:
    r1 = coordinates(696340, latlon[0], latlon[1])
    B_map.set_value(B_dipole(r1), latlon[0], latlon[1])

r_high = 800000
high_grid = grid(r_high, latitudes, longitudes, hs)
B_map_high = grid(r_high, latitudes, longitudes, hs)
B_map_comp = grid(r_high, latitudes, longitudes, hs)

for latlon in B_map_high.latlon:
    r1 = coordinates(r_high, latlon[0], latlon[1], latlon=True)
    B_c = np.asarray(B_comp(r1, base_grid, B_map))
    B_d = np.asarray(B_dipole(r1, returnBl=False, returnxyz=True))
    print(B_c, B_d)


r = coordinates(r_high, 0, 0, latlon=True)
B_p = B_comp(r, base_grid, B_map)
print(B_p)
