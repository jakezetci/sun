# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:52:18 2023

@author: cosbo
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from coordinates import coordinates
from lib import B_comp, grid
from field import B_dipole
import numpy as np
import math
from plots import sphere, disk


def plotmap(B_map, mode=disk):
    """
    Args:
        B_map (grid): a grid with either values of valuesvector.

    Returns:
        A plt.figure.
    """
    N = int(np.sqrt(B_map.lon.size))
    latlims = [B_map.lat.min(), B_map.lat.max()]
    lonlims = [B_map.lon.min(), B_map.lon.max()]
    fig, ax = plt.subplots(1, 1)
    n = B_map.num
    if np.linalg.norm(B_map.valuesvector[0]) == 0:
        n = B_map.num
        xx = np.empty(n)
        yy = np.empty(n)
        sq = np.empty(n)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.values)):
            x, y = cs.x, cs.y
            s = val
            xx[i], yy[i], sq[i] = x, y, s
        ax.plot(xx, yy, c=s)
    else:
        xx = np.empty(n)
        yy = np.empty(n)
        uu = np.empty(n)
        vv = np.empty(n)
        cc = np.empty(n)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.valuesvector)):
            x, y = cs.x, cs.y
            u, v, z = val
            xx[i], yy[i], uu[i], vv[i] = x, y, u, v
            cc[i] = math.sqrt(u**2 + v**2)
        ax.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))
    mode(ax, latlims, lonlims, n=N, r=B_map.r)
    plt.show()
