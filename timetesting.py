# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:37:25 2023

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
    n = (latlim[1] - latlim[0]) * N
    m = (lonlim[1] - lonlim[0]) * N
    lats = np.linspace(latlim[0], latlim[1], num=n)
    lons = np.linspace(lonlim[0], lonlim[1], num=m)
    latitudes = np.repeat(lats, m)
    longitudes = np.tile(lons, n)
    L = np.asarray([latitudes, longitudes])
    latitudes, longitudes = L
    r = 696340
    hs = (np.pi / 180) * (lonlim[1] - lonlim[0]) / n
    B_mapempty = grid(r, latitudes, longitudes, hs)
    return B_mapempty




latlimits = [0, 87]
lonlimits = [0, 90]

for i in range(latlims):
    