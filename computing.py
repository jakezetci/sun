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
from lib import B_comp, Grid, create_grid
from field import dipolebetter
from plots import sphere, disk, plotmap
import os
import time
import telebot
import alertbot

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

call_id = 412289661

def alert_bot(status_message, image=False,
              API_TOKEN='6559135670:AAHH8GpwVNX5k90FdFdWArasp9sc05fSLpI'):
    bot = telebot.TeleBot('6559135670:AAHH8GpwVNX5k90FdFdWArasp9sc05fSLpI')
    bot.send_message(call_id, status_message)
    if image is not False:
        bot.send_photo(call_id, image)


def crtname(name):
    nme = f'{np.random.randint(1,999)}'
    if (os.path.exists(f'checkpoint2 {name}.pkl') or
            os.path.exists(f'checkpoint1 {name}.pkl') or
            os.path.exists(f'{name}.pkl')):
        return crtname()
    else:
        return nme


def model_grid(B_map, dipole, dipolepos, vector=False, name=False,
               returnobj=False, checkpoints=721,
               alert=True):
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
    if vector is True:
        folder = 'vectormaps'
    else:
        folder = 'Lmaps'
    if name is False:
        name = crtname(name)

    checkpointsnum = int(checkpoints)
    B_map.save_pkl(name=f'checkpoint2 {name}', empty=False, vector=vector)
    B_map.save_pkl(name=f'checkpoint1 {name}', empty=False, vector=vector)
    r = B_map.r
    for i in range(B_map.progress, B_map.num):
        lat, lon = B_map.latlon[i]
        r1 = coordinates(r, lat, lon, latlon=True)
        B_map.set_value(dipolebetter(r1, m=dipole, rdipole=dipolepos,
                                     returnxyz=vector, returnBl=not vector),
                        lat, lon, index=i, vector=vector)
        B_map.progress1()
        if i % checkpointsnum == 0:
            os.remove(f'{folder}/checkpoint1 {name}.pkl')
            os.replace(f'{folder}/checkpoint2 {name}.pkl',
                       f'{folder}/checkpoint1 {name}.pkl')
            B_map.save_pkl(name=f'checkpoint2 {name}')
    B_map.save_pkl(name=name)
    if alert is True:
        alert_bot('модельная сетка посчиталась...')
    if returnobj:
        return B_map


def comp_grid(grid, B_map, vector=False, name=False,
              returnobj=False, checkpoints=721, timestamp=False, debug=False,
              alert=True):
    if vector is True:
        folder = 'vectormaps'
    else:
        folder = 'Lmaps'
    if name is False:
        name = crtname(name)
    checkpointsnum = int(checkpoints)
    grid.save_pkl(name=f'checkpoint2 {name}')
    grid.save_pkl(name=f'checkpoint1 {name}')
    R = grid.r
    tic = time.perf_counter()
    for i in range(Grid.progress, Grid.num):

        lat, lon = Grid.latlon[i]
        r1 = coordinates(R, lat, lon, latlon=True)
        value = B_comp(r1, B_map, debug=debug)
        grid.set_value(value, lat, lon, index=i, vector=vector)
        grid.progress1()
        if i % checkpointsnum == 0:
            os.remove(f'{folder}/checkpoint1 {name}.pkl')
            os.replace(f'{folder}/checkpoint2 {name}.pkl',
                       f'{folder}/checkpoint1 {name}.pkl')
            grid.save_pkl(name=f'checkpoint2 {name}')
        if timestamp:
            toc = time.perf_counter()
            print(
                f'values {i-checkpointsnum} - {i} done in {toc - tic:0.2f} seconds')
            tic = time.perf_counter()
    grid.save_pkl(name=name)
    if alert is True:
        alert_bot('расчётная сетка посчиталась...')
    if returnobj:
        return grid


def model_magneticline(magline, dipole, dipolepos, name=False, returnobj=False,
                       maxsteps=200, timestamp=False, alert=True):

    if name is False:
        name = f'{magline.initial_point.r:.2}-{magline.initial_value:.2} dipole {dipole}'

    magline.save_pkl(name=f'checkpoint {name}')
    start = magline.progress
    tic = time.perf_counter()
    if timestamp is not False:
        for i in range(start, maxsteps):
            magline.add_value(dipole, dipolepos)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f'values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds')
                tic = time.perf_counter()
                magline.save_pkl(name)
            else:
                magline.save_pkl(name=f'checkpoint1 {name}')
    else:
        for i in range(start, maxsteps):
            tic = time.perf_counter()
            magline.add_value(dipolebetter, dipole, dipolepos)
            magline.save_pkl(name=f'checkpoint1 {name}')

    magline.save_pkl(name)
    if alert is True:
        alert_bot('модельная магнитная линия посчиталась...')
    if returnobj:
        return magline


def comp_magneticline(magline, B_map, name=False, returnobj=False,
                      maxsteps=200, timestamp=False, alert=True):

    if name is False:
        name = crtname(name)

    magline.save_pkl(name=f'checkpoint {name}')
    start = magline.progress
    tic = time.perf_counter()
    if timestamp is not False:
        for i in range(start, maxsteps):

            magline.add_value_comp(B_map)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f'values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds')
                tic = time.perf_counter()
                magline.save_pkl(name)
            else:
                magline.save_pkl(name=f'checkpoint1 {name}')
    else:
        for i in range(start, maxsteps):
            tic = time.perf_counter()
            magline.add_value_comp(B_map)
            magline.save_pkl(name=f'checkpoint1 {name}')
    magline.save_pkl(name)
    if alert is True:
        alert_bot('расчётная магнитная линия посчиталась...')

    if returnobj:
        return magline
    


def create_model_plotmap(latlim, lonlim, N, M, dipolepos, vector=False,
                         name=False, n=10, alpha=0.7, lw=0.8):
    grid = create_grid(latlim, lonlim, N)
    computed = model_grid(grid, M, dipolepos, name=name, returnobj=True,
                          vector=vector)
    return plotmap(computed, lines=n, alpha=alpha, lw=lw)


if __name__ == "__main__":
    pass
