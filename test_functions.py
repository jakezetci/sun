# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:46:53 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import coordinates
from lib import B_comp, grid, ll2xyz
from field import dipolebetter
from plots import sphere
from plotting import plotmap
from computing import model_grid

import pickle


m = np.asarray(ll2xyz(1, 60, 30)) * 1e8
pos = coordinates(500000, 60, 30, latlon=True)
with open('bigmap.pkl', 'rb') as f:
    B_map = pickle.load(f)

newmap = model_grid(B_map, dipole=m,
                   dipolepos=pos, vector=True, name='testBIG', returnobj=True)
