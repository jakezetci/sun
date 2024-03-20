# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 20:26:15 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import Coordinates, ll2xyz, ll2pt, pt2xyz
from lib import B_comp_map, Grid, create_grid, Magneticline
from field import dipolebetter
from plots import sphere, disk, plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
from computing import create_model_plotmap, alert_bot
from plots import config
import pickle

dipolelat, dipolelon = 60, 30

phi, theta = ll2pt(60, 30)

alerts = False
computed = True

m = np.asarray([np.cos(-phi), 0, np.sin(phi)]) * 1e5*1e27
pos = Coordinates(600000*1e3, dipolelat, dipolelon, latlon=True)
pos = pos.vector

borderlat, borderlon = 25, 30

latlims = [dipolelat - borderlat, dipolelat + borderlat]
lonlims = [dipolelon - dipolelon, dipolelon + borderlon]
N = 2
if computed == False:

    empty_grid = create_grid(latlims, lonlims, N, name="perfectempty")
    Lmap = model_grid(empty_grid, dipole=m, dipolepos=pos,
                      vector=False, name='BIGMAP1010_', returnobj=True)
else:
    pass

steps = 10000
xmin, xmax = -8.74*1e7, 4.448e5*1e8
ymin, ymax = 5.12*1e8, 7.22*1e8
computed = True
if computed == True:
    with open('Lmaps/Диполь 60, 30 closeup22.pkl', 'rb') as fmap:
        picturemap = pickle.load(fmap)
    fig, ax = plotmap(picturemap, n_linesx=19, n_linesy=9, alpha=0.6, lw=0.6, xlabel='X, m', ylabel='Y, m',
                      title='Magnetic field lines of a dipole at (60,30)', ignoretop=True,
                      ignorecorners=5,
                      #xlimit=[xmin, xmax], ylimit=[ymin, ymax],
                      dpi=140,
                      figsize=(5.6, 5.6),
                      grid=False, 
                      )
else:
    latmax, latmin, lonmax, lonmin = 80, 40, 80, 10
    fig, ax = create_model_plotmap([49, 89], [-10, 80], N=10, M=m,
                                   dipolepos=pos, name='Диполь 60, 30 closeup22', n_linesx=17, n_linesy=9, alpha=0.6, lw=1.5,
                                   vector=False, xlabel='X, km', ylabel='Y, km',
                                   title='Magnetic lines of a dipole at (60,30)')


points = [[60, 36], [59, 38], [58, 32]]
colors = ['#029C63',  '#F6C3C3', '#11A0D7']


computed = True
if computed == False:

    for i, (lat, lon) in enumerate(points):
        point = Coordinates(696440*1e3, lat, lon, latlon=True)
        in_value = dipolebetter(point, dipolemoment=m,
                                rdipole=pos, returnxyz=True)
        line_model = Magneticline(point, in_value, step=30*1e3)
        line_comp = Magneticline(point, in_value, step=60*1e3)
        line_model = model_magneticline(line_model, m, pos, returnobj=True,
                                        name=f'presentable model line {lat} {lon}', maxsteps=steps,
                                        timestamp=1000, alert=alerts, stoppoint=696000*1e3)
        line_comp = comp_magneticline(line_comp, Lmap.values, Lmap.xyz, Lmap.area, returnobj=True,
                                      name=f'presentable comp line {lat} {lon}', maxsteps=steps,
                                      timestamp=100, alert=alerts, stoppoint=696000*1e3)

        xx_model, yy_model, zz_model = np.array(line_model.points).T
        xx_comp, yy_comp, zz_comp = np.array(line_comp.points).T
        ax.plot(xx_model, yy_model, '--', color=colors[i],
                label=f'model line {i+1}', lw=2)
        ax.plot(xx_comp, yy_comp, '-', color=colors[i],
                label=f'computed line {i+1}', lw=2)
else:
    for i, (lat, lon) in enumerate(points):
        with open(f'maglines/presentable model line {lat} {lon}.pkl', 'rb') as f_model:
            line_model = pickle.load(f_model)
        with open(f'maglines/presentable comp line {lat} {lon}.pkl', 'rb') as f_comp:
            line_comp = pickle.load(f_comp)
        xx_model, yy_model, zz_model = np.array(line_model.points).T
        xx_comp, yy_comp, zz_comp = np.array(line_comp.points).T
        ax.plot(xx_model, yy_model, '--', color=colors[i],
                lw=1.5)
        ax.plot(xx_comp, yy_comp, '-', color=colors[i],
                lw=1.5)
        xmax_line, ymax_line = np.max(np.hstack([xx_model, xx_comp])), np.max(
            np.hstack([yy_model, yy_comp]))
        xmin_line, ymin_line = np.min(np.hstack([xx_model, xx_comp])), np.min(
            np.hstack([yy_model, yy_comp]))
        if ymax_line > ymax:
            ymax = ymax_line
        if xmax_line > xmax:
            xmax = xmax_line
        if ymin_line < ymin:
            ymin = ymin_line
        if xmin_line < xmin:
            xmin = xmin_line
ax.plot([], [], '--', label='model lines of a dipole', color='black')
ax.plot([], [], '-', label='computed lines of a dipole', color='black')


ax.legend(loc='upper right', fontsize='medium')

if alerts:
    rand = np.random.randint(10, 99)
    fig.savefig(f'pics/lines{rand}.png')
    alert_bot('вот картинка..', imagepath=f'pics/lines{rand}.png')

plt.show()