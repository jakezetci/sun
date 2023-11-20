# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:16:44 2023

@author: cosbo
"""

"""
import cartopy.crs as ccrs
ax =plt.axes(projection=ccrs.Orthographic())
plt.pcolormesh(lons, lats,val, edgecolors='k', 
               linewidths=1, transform=ccrs.PlateCarree())

ax.coastlines()
ax.gridlines()

plt.show() """


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import Coordinates
from lib import B_comp_map, Grid, ll2xyz
from field import dipolebetter
from plots import sphere

from computing import model_grid
from plots import disk, plotmap

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


m = np.asarray(ll2xyz(60, 30, 1)) * 1e8
pos = Coordinates(500000, 60, 30, latlon=True)


with open("checkpoint1 testBIG.pkl", "rb") as f:
    map1 = pickle.load(f)
