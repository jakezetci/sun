# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 21:22:29 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import mplstereonet
from coordinates import coordinates
from lib import B_comp, grid
from field import dipolebetter
from plots import sphere
import os
import time

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


def crtname(name):
    nme = f'{np.random.randint(1,999)}'
    if (os.path.exists(f'checkpoint2 {name}.pkl') or
            os.path.exists(f'checkpoint1 {name}.pkl') or
            os.path.exists(f'{name}.pkl')):
        return crtname()
    else:
        return nme


def model_grid(B_map, dipole, dipolepos, vector=False, name=False,
               returnobj=False, checkpoints=721):
    """
    Model dipole grid computation with checkpoints

    Args:
        func (TYPE): DESCRIPTION.
        latlon (TYPE): DESCRIPTION.
        r (TYPE): DESCRIPTION.
        dipole (TYPE): DESCRIPTION.
        dipolepos (TYPE): DESCRIPTION.
        vector (TYPE, optional): DESCRIPTION. Defaults to False.
        name (TYPE, optional): DESCRIPTION. Defaults to False.
        returnobj (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        None.

    """

    if name is False:
        name = crtname(name)

    checkpointsnum = int(checkpoints)
    B_map.save_pkl(name=f'checkpoint2 {name}')
    B_map.save_pkl(name=f'checkpoint1 {name}')
    r = B_map.r
    for i in range(B_map.progress, B_map.num):
        lat, lon = B_map.latlon[i]
        r1 = coordinates(r, lat, lon, latlon=True)
        B_map.set_value(dipolebetter(r1, m=dipole, rdipole=dipolepos,
                                     returnxyz=vector, returnBl=not vector),
                        lat, lon, index=i, vector=vector)
        B_map.progress1()
        if i % checkpointsnum == 0:
            os.remove(f'checkpoint1 {name}.pkl')
            os.replace(f'checkpoint2 {name}.pkl', f'checkpoint1 {name}.pkl')
            B_map.save_pkl(name=f'checkpoint2 {name}')
    B_map.save_pkl(name=name)
    if returnobj:
        return B_map


def comp_grid(Grid, B_map, vector=False, name=False,
              returnobj=False, checkpoints=721, timestamp=False, debug=False):

    if name is False:
        name = crtname(name)
    checkpointsnum = int(checkpoints)
    Grid.save_pkl(name=f'checkpoint2 {name}')
    Grid.save_pkl(name=f'checkpoint1 {name}')
    R = Grid.r
    for i in range(Grid.progress, Grid.num):
        tic = time.perf_counter()
        lat, lon = Grid.latlon[i]
        r1 = coordinates(R, lat, lon, latlon=True)
        value = B_comp(r1, B_map, debug=debug)
        Grid.set_value(value, lat, lon, index=i, vector=vector)
        Grid.progress1()
        if i % checkpointsnum == 0:
            os.remove(f'checkpoint1 {name}.pkl')
            os.replace(f'checkpoint2 {name}.pkl', f'checkpoint1 {name}.pkl')
            Grid.save_pkl(name=f'checkpoint2 {name}')
        if timestamp:
            toc = time.perf_counter()
            print(f'value {i} done in {toc - tic:0.2f} seconds')
    Grid.save_pkl(name=name)
    if returnobj:
        return Grid


def model_magneticline(magline, dipole, dipolepos, name=False, returnobj=False,
                       maxsteps=200):

    if name is False:
        name = crtname(name)

    magline.save_pkl(name=f'checkpoint {name}')

    for i in range(maxsteps):
        magline.add_value(dipolebetter, dipole, dipolepos)
        magline.save_pkl(name=f'checkpoint1 {name}')
    magline.save_pkl(name)
    if returnobj:
        return magline


def comp_magneticline(magline, B_map, name=False, returnobj=False,
                      maxsteps=200, timestamp=False):

    if name is False:
        name = crtname(name)

    magline.save_pkl(name=f'checkpoint {name}')
    start = magline.progress
    for i in range(start, maxsteps):
        tic = time.perf_counter()
        magline.add_value_comp(B_map)
        magline.save_pkl(name=f'checkpoint1 {name}')
        if timestamp:
            toc = time.perf_counter()
            print(f'value {i} done in {toc - tic:0.2f} seconds')
    magline.save_pkl(name)
    if returnobj:
        return magline


if __name__ == "__main__":
    pass
