# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:16:44 2023

@author: cosbo
"""

import numpy as np
import math


a = np.asarray([1, 2, 3])
b = np.asarray([1, 4, 6])
x, y, z = b - a

try:
    print(math.sqrt(-1))
except ValueError:
    print("fafa")


latitudes , longitudes = np.loadtxt('lat-lon.txt')
latlon = list(map(list, zip(latitudes, longitudes)))
latlon = np.array(latlon)
rows = np.where(((latlon == [-90, -180]).all(1)))
print(latlon[rows])