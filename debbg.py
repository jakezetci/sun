# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 22:20:58 2023

@author: cosbo
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates, ll2xyz
from lib import B_comp, grid
from field import dipolebetter
from plots import sphere
from plotting import plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
import pickle
import copy


def create_grid(latlim, lonlim, N):
    # n - количество ячеек на градус
    n = (latlim[1] - latlim[0]) * N + 1
    m = (lonlim[1] - lonlim[0]) * N + 1
    lats = np.linspace(latlim[0], latlim[1], num=n)
    lons = np.linspace(lonlim[0], lonlim[1], num=m)
    latitudes = np.repeat(lats, m)
    longitudes = np.tile(lons, n)
    L = np.asarray([latitudes, longitudes])
    latitudes, longitudes = L
    r = 696340
    hs = (np.pi / 180) * (lonlim[1] - lonlim[0]) / m
    B_mapempty = grid(r, latitudes, longitudes, hs)
    return B_mapempty

m = np.asarray(ll2xyz(1, 10, 10)) * 1e12
pos = coordinates(600000, 10, 10, latlon=True)
pos = pos.vector


latlims = [-30, 30]
lonlims = [-20, 14]
N = 5
empty = create_grid(latlims, lonlims, N)

mapL = model_grid(empty, dipole=m, dipolepos=pos,
                  vector=False, name='Lexc', returnobj=True)



r = coordinates(696340+11000, 5, -12, latlon=True)
B, debug = B_comp(r, mapL, debug=True, change=False)
B2 = B_comp(r, mapL, debug=False, change=True)

print(B)
B_mod = dipolebetter(r, m=m, rdipole=pos, returnxyz=True)
print(B_mod)
debug = np.array(debug)

map_debug = copy.deepcopy(mapL)
map_debug2 = copy.deepcopy(mapL)
N = map_debug.num
for i in range(N):
    lat, lon = map_debug.latlon[i]
    map_debug.set_value(abs(debug[i][1]), lat, lon, index=i, vector=False)
    map_debug2.set_value(debug[i][1], lat, lon, index=i, vector=False)
j = np.argmax(debug[:,1])
print(map_debug.values[j], map_debug.latlon[j])
j = np.argmax(debug[:,0])
print(map_debug.values[j], map_debug.latlon[j])

plotmap(map_debug, lines=100)
plotmap(map_debug2, lines=100)
plotmap(mapL, lines=20)

plt.figure()
plt.plot(mapL.lat, debug[:, 1], 'o', ms=0.2)
yy = np.linspace(np.min(debug[:,1 ]), np.max(debug[:,1]))
xxdot = np.full_like(yy, 60)
xxdot2 = np.full_like(yy, 30)
plt.plot(xxdot, yy)

plt.title('lats')
plt.figure()
plt.plot(mapL.lon, debug[:, 1], 'o', ms=0.2)
plt.plot(xxdot2, yy)    
plt.title('lons')
plt.figure()
plt.plot(mapL.lat, mapL.values, 'o', ms=0.2)
plt.title('lats')
plt.figure()
plt.plot(mapL.lon, mapL.values, 'o', ms=0.2)
plt.title('lons')