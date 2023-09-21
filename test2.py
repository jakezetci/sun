# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 00:01:44 2023

@author: Home-Pc
"""
import numpy as np
from coordinates import coordinates
from lib import B_comp, grid
from field import B_dipole
import os



latitudes, longitudes = np.loadtxt('lat-lon_big.txt')
r = 700000
N = int(np.unique(longitudes).size)

lonmin, lonmax = longitudes.min(), longitudes.max()

hs = (np.pi / 180) * (lonmax - lonmin) / N
B_map = grid(r, latitudes, longitudes, hs)
B_map.save_pkl('bigmap')