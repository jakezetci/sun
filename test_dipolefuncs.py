# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 16:10:27 2023

@author: cosbo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates, ll2xyz, ll2pt, pt2xyz
from lib import B_comp, Grid, create_grid, Magneticline
from field import dipolebetter, twomonopoles
from plots import sphere, disk, plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
from computing import create_model_plotmap, alert_bot
from plots import config
import pickle


alerts = False
dipolelat, dipolelon = 60, 30

phi, theta = ll2pt(60, 30)

m = np.asarray([np.cos(-phi), 0, np.sin(phi)]) * 1e12


pos = coordinates(600000, dipolelat, dipolelon, latlon=True)
pos = pos.vector
pos = [0, 0, 0]
M = np.asarray([0, 2, 0])

q1, q2 = 1e18, -1e18
rq1, rq2 = [0, 4, 0], [0, -4, 0]

relativelat1, relativelon1 = -2.857142857142854, +2.8571428571428577
relativelat, relativelon = -30, -15
n = 50
m = 50
xs = np.linspace(-1000, 1000, num=n)
ys = np.linspace(-1000, 1000, num=m)
pointsx = np.repeat(xs, m)
pointsy = np.tile(ys, n)
coor_set = []
for x, y in zip(pointsx, pointsy):
    coor_set.append(coordinates(x, y, 0))

uu, vv, ww = [], [], []
for coor in coor_set:
    u, v, w = dipolebetter(coor, M, returnxyz=True)
    uu.append(u)
    vv.append(v)
    ww.append(w)

fig, ax = plt.subplots(1, 1)

ax.quiver(pointsx, pointsy, uu, vv, np.arctan2(uu, vv))

point1 = coordinates(5, 5, 0)
in_value = dipolebetter(point1, M, returnxyz=True)
line = Magneticline(point1, in_value, step=1)
steps = 100000
line = model_magneticline(line, M, pos, returnobj=True,
                          name='model newdipole v6_5', maxsteps=steps,
                          timestamp=10000, alert=alerts)

xx_model, yy_model, zz_model = np.array(line.pointsxyz).T

ax.plot(xx_model, yy_model, 'o', color='pink', label='model', ms=0.2)
