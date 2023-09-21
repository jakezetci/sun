# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 23:34:15 2023

@author: Home-Pc
"""

import numpy as np

n = 361
m = 361
"""
lats = np.hstack((np.linspace(-30, -15, num=n),
                 np.linspace(15, 30, num=n)))
lons = np.hstack((np.linspace(-30, -15, num=m),
                 np.linspace(15, 30, num=m)))
"""
lats = np.linspace(-90, 90, num=n)
lons = np.linspace(-90, 90, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_big.txt", L)


