# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 22:56:50 2023

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
alerts = False

m = np.asarray([np.cos(-phi), 0, np.sin(phi)]) * 1e5
#m = np.asarray(ll2xyz(60, 30, 1)) * 1e5
pos = coordinates(600000, dipolelat, dipolelon, latlon=True)
pos = pos.vector
relativelat1, relativelon1 = 6, 8

borderlat, borderlon = 22, 30
"""
latlims = [dipolelat-borderlat+relativelat1, dipolelat+borderlat+relativelat1]
lonlims = [dipolelon-borderlon+relativelon1, dipolelon+borderlon+relativelon1]
"""

pointlat, pointlon = 52, 20

perfect_point = coordinates(
    696340+100, pointlat, pointlon, latlon=True)
latlims = [pointlat - borderlat, pointlat + borderlat]
lonlims = [pointlon - pointlon, pointlon + borderlon]


steps = 2000
with open('Lmaps/Диполь 60, 30.pkl', 'rb') as fmap:
    picturemap = pickle.load(fmap)
fig, ax = plotmap(picturemap, lines=20, alpha=0.6, lw=1, ignoretop=True)

lat_array = np.linspace(45, 61, num=17)
pointlat, pointlon = 52, 30
points = [[60, 36], [59, 38], [58, 32]]
for lat, lon in points:
    perfect_point = coordinates(
        696340+100, lat, lon, latlon=True)
    in_value = dipolebetter(perfect_point, dipolemoment=m,
                            rdipole=pos, returnxyz=True)
    line_model = Magneticline(perfect_point, in_value, step=100)
    line_model = model_magneticline(line_model, m, pos, returnobj=True,
                                    name='model test', maxsteps=int(steps*3),
                                    timestamp=1000, alert=alerts, stoppoint=696000)

    xx_model, yy_model, zz_model = np.array(line_model.pointsxyz).T

    ax.plot(xx_model, yy_model, '--',
            label=f'model {lat} {lon}', lw=2)


ax.legend(loc='best', fontsize='x-large')

if alerts:
    rand = np.random.randint(10, 99)
    fig.savefig(f'lines{rand}.png')
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
