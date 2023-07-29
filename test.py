# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 23:34:15 2023

@author: Home-Pc
"""

import numpy as np

n = 15
m = 30
lats = np.linspace(-90, 90, num=n)
lons = np.linspace(-180, 180, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon.txt", L)
