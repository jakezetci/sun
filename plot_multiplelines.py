# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 20:26:15 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates, ll2xyz, ll2pt, pt2xyz
from lib import B_comp_map, Grid, create_grid, Magneticline
from field import dipolebetter
from plots import sphere, disk, plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
from computing import create_model_plotmap, alert_bot
from plots import config
import pickle

dipolelat, dipolelon = 60, 30

phi, theta = ll2pt(60, 30)

alerts = True
computed = True

m = np.asarray([np.cos(-phi), 0, np.sin(phi)]) * 1e5
pos = coordinates(600000, dipolelat, dipolelon, latlon=True)
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

steps = 5000
xmin, xmax = -1e4, 5.5e5
ymin, ymax = 4.4e5, 7.4e5
if computed == True:
    with open('Lmaps/Диполь 60, 30 closeup.pkl', 'rb') as fmap:
        picturemap = pickle.load(fmap)
    fig, ax = plotmap(picturemap, lines=21, alpha=0.6, lw=1.5, xlabel='X, km', ylabel='Y, km',
                      title='Magnetic lines of a dipole at (60,30)', ignoretop=True,
                      xlimit=[xmin, xmax], ylimit=[ymin, ymax])
else:
    latmax, latmin, lonmax, lonmin = 80, 40, 85, 5
    fig, ax = create_model_plotmap([latmin, latmax], [lonmin, lonmax], N=10, M=m,
                                   dipolepos=pos, name='Диполь 60, 30 closeup', lines=20, alpha=0.6, lw=1.5,
                                   vector=False, xlabel='X, km', ylabel='Y, km',
                                   title='Magnetic lines of a dipole at (60,30)')


points = [[60, 36], [59, 38], [58, 32]]
colors = ['green', 'pink', 'cyan']

if computed == False:

    for i, (lat, lon) in enumerate(points):
        point = coordinates(696340+100, lat, lon, latlon=True)
        in_value = dipolebetter(point, dipolemoment=m,
                                rdipole=pos, returnxyz=True)
        line_model = Magneticline(point, in_value, step=300)
        line_comp = Magneticline(point, in_value, step=600)
        line_model = model_magneticline(line_model, m, pos, returnobj=True,
                                        name=f'presentable model line {lat} {lon}', maxsteps=steps,
                                        timestamp=1000, alert=alerts, stoppoint=696000)
        line_comp = comp_magneticline(line_comp, Lmap, returnobj=True,
                                      name=f'presentable comp line {lat} {lon}', maxsteps=steps,
                                      timestamp=100, alert=alerts, stoppoint=696000)

        xx_model, yy_model, zz_model = np.array(line_model.pointsxyz).T
        xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
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
        xx_model, yy_model, zz_model = np.array(line_model.pointsxyz).T
        xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
        ax.plot(xx_model, yy_model, '--', color=colors[i],
                label=f'model line {i+1}', lw=2)
        ax.plot(xx_comp, yy_comp, '-', color=colors[i],
                label=f'computed line {i+1}', lw=2)
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


ax.legend(loc='best', fontsize='x-large')

if alerts:
    rand = np.random.randint(10, 99)
    fig.savefig(f'pics/lines{rand}.png')
    alert_bot('вот картинка..', imagepath=f'pics/lines{rand}.png')
