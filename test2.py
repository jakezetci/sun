# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 00:01:44 2023

@author: Home-Pc
"""
import numpy as np
from coordinates import coordinates
from lib import B_comp, grid
from field import B_dipole

latitudes , longitudes = np.loadtxt('lat-lon.txt')

a, b = np.asarray([3, 2, 1]), np.asarray([1, 1, 1])
x, y, z = a
print(a/4)