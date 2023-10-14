# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:46:53 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates, ll2xyz, xy2ll
from lib import B_comp, Grid, create_grid, Magneticline
from field import dipolebetter
from plots import sphere, disk, plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
from computing import create_model_plotmap
from plots import config
import pickle

dipolelat, dipolelon = 60, 30
m = np.asarray(ll2xyz(10, 10, 1)) * 1e5
pos = coordinates(600000, dipolelat, dipolelon, latlon=True)
pos = pos.vector
relativelat1, relativelon1 = 2.857142857142854, -2.8571428571428577
relativelat, relativelon = -8, -4
borderlat, borderlon = 22, 30
latlims = [dipolelat-borderlat+relativelat, dipolelat+borderlat+relativelat]
lonlims = [dipolelon-borderlon+relativelon, dipolelon+borderlon+relativelon]
N = 2
empty_grid = create_grid(latlims, lonlims, N, name="perfectempty")
Lmap = model_grid(empty_grid, dipole=m, dipolepos=pos,
                    vector=False, name='BIGMAP1010_', returnobj=True)

print('Lmapdone')

perfect_point = coordinates(696340+100, dipolelat+relativelat, dipolelon+relativelon, latlon=True)
in_value = dipolebetter(perfect_point, m=m, rdipole=pos, returnxyz=True)
Bc = in_value
Bd = B_comp(perfect_point, Lmap, debug=False, change=False)
print(Bc, Bd)
"""
line_model = Magneticline(perfect_point, in_value, step=100)

line_comp =  Magneticline(perfect_point, in_value, step=100)
steps = 500

line_model = model_magneticline(line_model, m, pos, returnobj=True, 
                                name='picture_model7', maxsteps=steps,
                                timestamp=100)

line_comp = comp_magneticline(line_comp, Lmap, returnobj=True,
                              name='picture_comp7', maxsteps=steps,
                              timestamp=50)
"""


with open('maglines/picture_model2.pkl', 'rb') as f_model:
    line_model = pickle.load(f_model)
with open('maglines/picture_comp2.pkl', 'rb') as f_comp:
    line_comp = pickle.load(f_comp)



xx_model, yy_model, zz_model = np.array(line_model.pointsxyz).T
xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
xmax, ymax = np.max(np.hstack([xx_model, xx_comp])), np.max(np.hstack([yy_model, yy_comp]))
xmin, ymin = np.min(np.hstack([xx_model, xx_comp])), np.min(np.hstack([yy_model, yy_comp]))
"""
latmax, lonmax = xy2ll(xmax, ymax)
latmin, lonmin = xy2ll(xmin, ymin)
"""
"""
latmax, latmin, lonmax, lonmin = 80, 10, 70, -10

fig, ax = create_model_plotmap([latmin, latmax], [lonmin, lonmax], N=4, M=m,
                               dipolepos=pos, name='Диполь 60, 30', n=5, alpha=0.6, lw=1.5)
"""
with open('Lmaps/Диполь 60, 30.pkl', 'rb') as fmap:
    picturemap = pickle.load(fmap)

fig, ax = plotmap(picturemap, lines=5, alpha=0.6, lw=1)

"""

fig, ax = config(grid=False)

disk(ax, [latmin, latmax], [lonmin, lonmax], n=10)
"""

ax.plot(xx_model, yy_model, '-',
        color='pink', label='model', lw=4)
ax.plot(xx_comp, yy_comp, '--',
        color='green', label='comp', lw=4)
ax.legend(loc='best', fontsize='x-large')






"""
with open('smallempty.pkl', 'rb') as ff:
    empty_small = pickle.load(ff)
    empty_small.change_coors()


map_comp = comp_grid(empty_small, Lmap, vector=True,
                     name='small_comp_test3', returnobj=True, timestamp=True)




with open('smallempty.pkl', 'rb') as ff:
    empty_small = pickle.load(ff)
    empty_small.change_coors()

map_model = model_grid(empty_small, dipole=m, dipolepos=pos,
                       vector=True, name='small_model3', returnobj=True)
"""

"""
with open('small_model3.pkl', 'rb') as ff:
    map_model = pickle.load(ff)

with open('small_comp_test3.pkl', 'rb') as fff:
    map_comp = pickle.load(fff)

values= np.array(map_comp.valuesvector)
j = np.argmax(values[:,1])
print(values[j], map_comp.latlon[j])


with open('BIG_bigL.pkl', 'rb') as f:
    B_map_full = pickle.load(f)
    B_map_full.change_coors()
    B_map_full.save_pkl('BIG_bigL')


with open('smallempty.pkl', 'rb') as ff:
    empty_small = pickle.load(ff)

with open('small_model2.pkl', 'rb') as ff:
    map_model = pickle.load(ff)

map_comp = comp_grid(empty_small, newmap, vector=True,
                     name='small_comp_test', returnobj=True, timestamp=True)

with open('small_comp_test2.pkl', 'rb') as fff:
    map_comp = pickle.load(fff)

plotmap(map_comp, limit=0.1)
plotmap(map_model)
"""