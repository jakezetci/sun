# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:46:53 2023

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


m = np.asarray(ll2xyz(1, 10, 10)) * 1e12
pos = coordinates(600000, 10, 10, latlon=True)
pos = pos.vector


with open('bigempty_.pkl', 'rb') as f:
    empty = pickle.load(f)
    empty.change_coors()

Lmap = model_grid(empty, dipole=m, dipolepos=pos,
                    vector=False, name='BIGMAP1010_', returnobj=True)

print('Lmapdone')
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