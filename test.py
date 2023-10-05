# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 23:34:15 2023

@author: Home-Pc
"""

import numpy as np

n = 341
m = 361
"""
lats = np.hstack((np.linspace(-30, -15, num=n),
                 np.linspace(15, 30, num=n)))
lons = np.hstack((np.linspace(-30, -15, num=m),
                 np.linspace(15, 30, num=m)))
"""
lats = np.linspace(-85, 85, num=n)
lons = np.linspace(-90, 90, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_half_.txt", L)

n = 11
m = 11

lats = np.linspace(55, 65, num=n)
lons = np.linspace(25, 35, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_small.txt", L)

n = 3
m = 3

lats = np.linspace(59, 61, num=n)
lons = np.linspace(29, 31, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_extrasmall.txt", L)

n = 341
m = 721
lats = np.linspace(-85, 85, num=n)
lons = np.linspace(-180, 180, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_big__.txt", L)

n = 99
m = 99

lats = np.linspace(35, 65, num=n)
lons = np.linspace(15, 50, num=m)
latitudes = np.repeat(lats, m)
longitudes = np.tile(lons, n)
L = np.asarray([latitudes, longitudes])
np.savetxt("lat-lon_exl.txt", L)