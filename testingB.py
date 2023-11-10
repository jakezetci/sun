# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 22:47:33 2023

@author: Home-Pc
"""

import numpy as np
from coordinates import coordinates
from lib import B_comp_map, grid
from field import B_dipole

latitudes , longitudes = np.loadtxt('lat-lon.txt')
hs = 90/15
r = 696340
base_grid = grid(r, latitudes, longitudes, hs)

B_map = grid(r, latitudes, longitudes, hs)

a = (B_map.latlon)

for latlon in B_map.latlon:
    r1 = coordinates(696340, latlon[0], latlon[1])
    B_map.set_value(B_dipole(r1), latlon[0], latlon[1])


r = coordinates(700000, 0, 0, latlon=True)
B_p = B_comp_map(r, base_grid, B_map)
print(B_p)
