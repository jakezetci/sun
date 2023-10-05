# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 00:01:44 2023

@author: Home-Pc
"""
import numpy as np
from coordinates import coordinates, ll2xyz
from lib import B_comp, grid, magneticline
from field import B_dipole, dipolebetter
import os
"""

latitudes, longitudes = np.loadtxt('lat-lon_big.txt')
N = int(np.size(np.unique(longitudes)))
latlims = [latitudes.min(), latitudes.max()]
lonlims = [longitudes.min(), longitudes.max()]

hs = (np.pi / 180) * (lonlims[1] - lonlims[0]) / N
r = 696340
B_mapempty = grid(r, latitudes, longitudes, hs)
B_mapempty.save_pkl('bigempty')

latitudes, longitudes = np.loadtxt('lat-lon_small.txt')
N = int(np.size(np.unique(longitudes)))
latlims2 = [latitudes.min(), latitudes.max()]
lonlims2 = [longitudes.min(), longitudes.max()]

hs = (np.pi / 180) * (lonlims2[1] - lonlims2[0]) / N
r = 800000
B_mapempty = grid(r, latitudes, longitudes, hs)
B_mapempty.save_pkl('smallempty')

latitudes, longitudes = np.loadtxt('lat-lon_extrasmall.txt')
N = int(np.size(np.unique(longitudes)))
latlims2 = [latitudes.min(), latitudes.max()]
lonlims2 = [longitudes.min(), longitudes.max()]
r = 700000
hs = (np.pi / 180) * (lonlims2[1] - lonlims2[0]) / N
B_mapempty = grid(r, latitudes, longitudes, hs)
B_mapempty.save_pkl('exsmallempty')

latitudes, longitudes = np.loadtxt('lat-lon_half_.txt')
N = int(np.size(np.unique(longitudes)))
latlims2 = [latitudes.min(), latitudes.max()]
lonlims2 = [longitudes.min(), longitudes.max()]
r = 696340
hs = (np.pi / 180) * (lonlims2[1] - lonlims2[0]) / N
B_mapempty = grid(r, latitudes, longitudes, hs)
B_mapempty.save_pkl('halfempty_')
"""
latitudes, longitudes = np.loadtxt('lat-lon_exl.txt')
N = int(np.size(np.unique(longitudes)))
latlims2 = [latitudes.min(), latitudes.max()]
lonlims2 = [longitudes.min(), longitudes.max()]
r = 696340
hs = (np.pi / 180) * (lonlims2[1] - lonlims2[0]) / N
B_mapempty = grid(r, latitudes, longitudes, hs)
B_mapempty.save_pkl('empty_exl')
