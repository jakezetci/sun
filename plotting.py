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
import pickle


def plotmap(B_map, mode=disk, limit=False, every=1, lines=1):
    """
    Args:
        B_map (grid): a grid with either values of valuesvector.

    Returns:
        A plt.figure.
    """
    N = int(np.size(np.unique(B_map.lon)))
    latlims = [B_map.lat.min(), B_map.lat.max()]
    lonlims = [B_map.lon.min(), B_map.lon.max()]
    fig, ax = plt.subplots(1, 1)
    n = B_map.num
    if np.linalg.norm(B_map.valuesvector[0]) == 0:
        n = B_map.num
        xx = np.zeros(n//every + 1)
        yy = np.zeros(n//every + 1)
        sq = np.zeros(n//every + 1)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.values)):
            if i % every == 0:
                
                x, y = cs.x, cs.y
                s = val
                xx[i//every], yy[i//every], sq[i//every] = x, y, s
        ax.scatter(xx, yy, c=sq, s=0.5, cmap='inferno')
    else:
        xx = np.zeros(n)
        yy = np.zeros(n)
        uu = np.zeros(n)
        vv = np.zeros(n)
        cc = np.zeros(n)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.valuesvector)):
            x, y = cs.x, cs.y
            u, v, z = val
            if limit is not False:
                if abs(u)>limit or abs(v)>limit:
                    u, v = 0, 0
            xx[i], yy[i], uu[i], vv[i] = x, y, u, v
            cc[i] = math.sqrt(u**2 + v**2)
        ax.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))
    mode(ax, latlims, lonlims, n=N//lines, r=B_map.r)
    fig.show()

if __name__ == "__main__":
    """
    with open('maglinemodel.pkl', 'rb') as fff:
        line1 = pickle.load(fff)
    with open('maglinecomp2.pkl', 'rb') as fff:
        line2 = pickle.load(fff)
    fig2, ax3 = plt.subplots(1, 1)
    ax3.set_title("расчёты...", fontsize=20, fontname='HSE Slab')
    x1, y1 = np.array(line1.pointsxyz)[:,0], np.array(line1.pointsxyz)[:,1]
    points2 = np.array(line2.pointsxyz)[:,:2]


    x2, y2 = np.array(line2.pointsxyz)[:,0], np.array(line2.pointsxyz)[:,1]
    ax3.plot(x1, y1,'o', lw=0.5)
    ax3.plot(x2, y2, 'o', lw=0.5)
    """
    with open('BIGL.pkl', 'rb') as fff:
        bigmap = pickle.load(fff)
        bigmap.change_coors()
        bigmap.save_pkl('BIGL')

    
    plotmap(bigmap, every=2)