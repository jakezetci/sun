# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:46:53 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates, ll2xyz, ll2pt, pt2xyz
from lib import B_comp, Grid, create_grid, Magneticline
from field import dipolebetter
from plots import sphere, disk, plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
from computing import create_model_plotmap, alert_bot
from plots import config
import pickle

dipolelat, dipolelon = 60, 30

phi, theta = ll2pt(60, 30)
alerts = True

m = np.asarray(ll2xyz(60, 30, 1)) * 1e5
pos = coordinates(600000, dipolelat, dipolelon, latlon=True)
pos = pos.vector
relativelat1, relativelon1 = 6, 8
relativelat, relativelon = -17, -15
borderlat, borderlon = 22, 30
latlims = [dipolelat-borderlat+relativelat1, dipolelat+borderlat+relativelat1]
lonlims = [dipolelon-borderlon+relativelon1, dipolelon+borderlon+relativelon1]
N = 2
empty_grid = create_grid(latlims, lonlims, N, name="perfectempty")
Lmap = model_grid(empty_grid, dipole=m, dipolepos=pos,
                  vector=False, name='BIGMAP1010_', returnobj=True)

print('Lmapdone')

perfect_point = coordinates(
    696340+100, dipolelat+relativelat1, dipolelon+relativelon1, latlon=True)
in_value = dipolebetter(perfect_point, dipolemoment=m,
                        rdipole=pos, returnxyz=True)
Bc = in_value
Bd = B_comp(perfect_point, Lmap, debug=False, change=False)
print(Bc, Bd)


steps = 2000
line_model = Magneticline(perfect_point, in_value, step=100)

line_comp = Magneticline(perfect_point, in_value, step=300)


line_model = model_magneticline(line_model, m, pos, returnobj=True,
                                name='model good dipole 60 30 start: 66 38', maxsteps=int(steps*3),
                                timestamp=200, alert=alerts, stoppoint=696000)

line_comp = comp_magneticline(line_comp, Lmap, returnobj=True,
                              name='comp good dipole 60 30 start: 66 38', maxsteps=steps,
                              timestamp=100, alert=alerts, stoppoint=696000)


"""




with open('maglines/model newdipole v6.pkl', 'rb') as f_model:
    line_model = pickle.load(f_model)
with open('maglines/comp newdipole v6.pkl', 'rb') as f_comp:
    line_comp = pickle.load(f_comp)


"""

xx_model, yy_model, zz_model = np.array(line_model.pointsxyz).T
xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
xmax, ymax = np.max(np.hstack([xx_model, xx_comp])), np.max(
    np.hstack([yy_model, yy_comp]))
xmin, ymin = np.min(np.hstack([xx_model, xx_comp])), np.min(
    np.hstack([yy_model, yy_comp]))
"""
latmax, lonmax = xy2ll(xmax, ymax)
latmin, lonmin = xy2ll(xmin, ymin)
"""


with open('Lmaps/Диполь 60, 30.pkl', 'rb') as fmap:
    picturemap = pickle.load(fmap)
fig, ax = plotmap(picturemap, lines=20, alpha=0.6, lw=1)
"""

latmax, latmin, lonmax, lonmin = 80, 10, 70, -10
fig, ax = create_model_plotmap([latmin, latmax], [lonmin, lonmax], N=5, M=m,
                               dipolepos=pos, name='Диполь 60, 30', lines=10, alpha=0.6, lw=1.5,
                               vector=False)


fig, ax = create_model_plotmap([latmin, latmax], [lonmin, lonmax], N=5, M=m,
                               dipolepos=pos, name='Диполь 60, 30 new', n=5, alpha=0.6, lw=1.5,
                               vector=False)
"""
"""

fig, ax = config(grid=False)

disk(ax, [latmin, latmax], [lonmin, lonmax], n=10)
"""

ax.plot(xx_model, yy_model, '--',
        color='pink', label='model', lw=2)
ax.plot(xx_comp, yy_comp, '-',
        color='green', label='comp', lw=2)
ax.legend(loc='best', fontsize='x-large')

if alerts:
    rand = np.random.randint(10, 99)
    fig.savefig('lines{rand}.png')
    alert_bot('вот картинка..', imagepath='lines2.png')
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
