# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 21:22:29 2023

@author: cosbo
"""

import numpy as np


import os
import time
import telebot

try:
    from coordinates import coordinates
    from lib import B_comp_map, Grid, create_grid, Magneticline, B_comp
    from field import dipolebetter
    from plots import sphere, disk, plotmap
except ModuleNotFoundError:
    from sun.coordinates import coordinates
    from sun.lib import B_comp_map, Grid, create_grid, Magneticline, B_comp
    from sun.field import dipolebetter
    from sun.plots import sphere, disk, plotmap

import pickle

with open('callid.txt') as f:
    call_id = int(f.read())

with open('TOKEN.txt') as ff:
    API_TOKEN = ff.read()


def create_grid(latlim: tuple, lonlim: tuple, N, r=696340 * 1000, name=False):
    # N - количество ячеек на градус
    n = int((latlim[1] - latlim[0]) * N)
    m = int((lonlim[1] - lonlim[0]) * N)
    lats = np.linspace(latlim[0], latlim[1], num=n)
    lons = np.linspace(lonlim[0], lonlim[1], num=m)
    latitudes = np.repeat(lats, m)
    longitudes = np.tile(lons, n)
    L = np.asarray([latitudes, longitudes])
    latitudes, longitudes = L
    hs = (np.pi / 180) * (lonlim[1] - lonlim[0]) / m
    B_mapempty = Grid(r, latitudes, longitudes, hs)
    if name is not False:
        B_mapempty.save_pkl(name, empty=True)
    return B_mapempty


def alert_bot(
    status_message,
    imagepath=False
):
    bot = telebot.TeleBot(API_TOKEN)
    bot.send_message(call_id, status_message)
    if imagepath is not False:
        bot.send_photo(call_id, telebot.types.InputFile(imagepath))


def crtname(name):
    nme = f"{np.random.randint(1,999)}"
    if (
        os.path.exists(f"checkpoint2 {name}.pkl")
        or os.path.exists(f"checkpoint1 {name}.pkl")
        or os.path.exists(f"{name}.pkl")
    ):
        return crtname()
    else:
        return nme


def model_grid(
    B_map,
    dipole,
    dipolepos,
    vector=False,
    name=False,
    returnobj=False,
    checkpoints=10000,
    alert=True,
):
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
        folder = "vectormaps"
    else:
        folder = "Lmaps"
    if name is False:
        name = crtname(name)

    checkpointsnum = int(checkpoints)
    B_map.save_pkl(name=f"checkpoint2 {name}", empty=False, vector=vector)
    B_map.save_pkl(name=f"checkpoint1 {name}", empty=False, vector=vector)
    r = B_map.r
    for i in range(B_map.progress, B_map.num):
        lat, lon = B_map.latlon[i]
        r1 = coordinates(r, lat, lon, latlon=True)
        B_map.set_value(
            dipolebetter(
                r1,
                dipolemoment=dipole,
                rdipole=dipolepos,
                returnxyz=vector,
                returnBl=not vector,
            ),
            lat,
            lon,
            index=i,
            vector=vector,
        )
        B_map.progress1()
        if i % checkpointsnum == 0:
            os.remove(f"{folder}/checkpoint1 {name}.pkl")
            os.replace(
                f"{folder}/checkpoint2 {name}.pkl", f"{folder}/checkpoint1 {name}.pkl"
            )
            B_map.save_pkl(name=f"checkpoint2 {name}")
    B_map.save_pkl(name=name)
    if alert is True:
        alert_bot("модельная сетка посчиталась...")
    if returnobj:
        return B_map


def comp_grid(
    grid: Grid,
    B_map: Grid,
    vector=False,
    name=False,
    returnobj=False,
    checkpoints=721,
    timestamp=False,
    debug=False,
    alert=True,
):
    if vector is True:
        folder = "vectormaps"
    else:
        folder = "Lmaps"
    if name is False:
        name = crtname(name)
    checkpointsnum = int(checkpoints)
    grid.save_pkl(name=f"checkpoint2 {name}")
    grid.save_pkl(name=f"checkpoint1 {name}")
    R = grid.r
    tic = time.perf_counter()
    for i in range(grid.progress, grid.num):
        lat, lon = grid.latlon[i]
        r1 = grid.coors_set[i]
        value = B_comp_map(r1, B_map, debug=debug)
        grid.set_value(value, index=i, vector=vector)
        grid.progress1()
        if i % checkpointsnum == 0:
            os.remove(f"{folder}/checkpoint1 {name}.pkl")
            os.replace(
                f"{folder}/checkpoint2 {name}.pkl", f"{folder}/checkpoint1 {name}.pkl"
            )
            grid.save_pkl(name=f"checkpoint2 {name}")
        if timestamp:
            toc = time.perf_counter()
            print(
                f"values {i-checkpointsnum} - {i} done in {toc - tic:0.2f} seconds")
            tic = time.perf_counter()
    grid.save_pkl(name=name)
    if alert is True:
        alert_bot("расчётная сетка посчиталась...")
    if returnobj:
        return grid


def model_magneticline(
    magline: Magneticline,
    dipole,
    dipolepos,
    name=False,
    returnobj=False,
    maxsteps=2000,
    timestamp=10000,
    stoppoint=None,
    alert=True,
):
    def stoppoint_model():
        start = magline.progress
        tic = time.perf_counter()
        for i in range(start, maxsteps):
            magline.add_value(dipole, dipolepos)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f"values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
                magline.save_pkl(name)
            if magline.points[-1].r < stoppoint:
                break
        return magline

    def steps_model():
        start = magline.progress
        tic = time.perf_counter()
        for i in range(start, maxsteps):
            magline.add_value(dipole, dipolepos)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f"values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
                magline.save_pkl(name)
        return magline

    if name is False:
        name = (
            f"{magline.initial_point.r:.2}-{magline.initial_value:.2} dipole {dipole}"
        )

    start = magline.progress
    tic = time.perf_counter()
    if stoppoint is None:
        magline = steps_model()
    else:
        magline = stoppoint_model()
    magline.save_pkl(name)
    if alert is True:
        alert_bot("модельная магнитная линия посчиталась...")
    if returnobj:
        return magline


def comp_magneticline(
    magline: Magneticline,
    B_map,
    name=False,
    returnobj=False,
    maxsteps=2000,
    timestamp=1000,
    alert=True,
    stoppoint=False,
):
    def steps_comp():
        start = magline.progress
        tic = time.perf_counter()
        for i in range(start, maxsteps):
            magline.add_value_comp(B_map)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f"values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
                magline.save_pkl(name)
        return magline

    def stoppoint_comp():
        start = magline.progress
        tic = time.perf_counter()
        for i in range(start, maxsteps):
            magline.add_value_comp(B_map)
            if i % timestamp == 0:
                toc = time.perf_counter()
                print(
                    f"values {i-timestamp}-{i} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
                magline.save_pkl(name)
            if magline.points[-1].r < stoppoint:
                break
        return magline

    if name is False:
        name = crtname(name)

    if stoppoint is None:
        magline = steps_comp()
    else:
        magline = stoppoint_comp()
    magline.save_pkl(name)
    if alert is True:
        alert_bot("расчётная магнитная линия посчиталась...")

    if returnobj:
        return magline


def create_model_plotmap(
    latlim,
    lonlim,
    N,
    M,
    dipolepos,
    vector=False,
    name=False,
    lines=10,
    alpha=0.7,
    lw=0.8,
    title="",
    xlabel="xlabel",
    ylabel="ylabel",
):
    grid = create_grid(latlim, lonlim, N)
    computed = model_grid(grid, M, dipolepos, name=name,
                          returnobj=True, vector=vector)
    return plotmap(
        computed,
        lines=lines,
        alpha=alpha,
        lw=lw,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )


def comp_grid_points(
    grid: Grid,
    points: np.array,
    values: np.array,
    name=False,
    returnobj=True,
    checkpoints=721,
    timestamp=False,
    debug=False,
    alert=True,
):
    if name is False:
        name = crtname(name)
    checkpointsnum = int(checkpoints)
    tic = time.perf_counter()
    for i in range(grid.progress, grid.num):
        r1 = grid.coors_set[i]
        value = B_comp(r1, values, points)
        grid.set_value(value, index=i, vector=True)
        grid.progress1()
        if i % checkpointsnum == 0:

            if timestamp:
                toc = time.perf_counter()
                print(
                    f"values {i-checkpointsnum} - {i} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
    grid.save_pkl(name=name)
    if alert is True:
        alert_bot("расчётная сетка посчиталась...")
    if returnobj:
        return grid


def compute_grid_energy(grid: Grid):
    energy = 0
    for value, area in zip(grid.valuesvector, grid.area):
        energy = energy + np.linalg.norm(value) * area
    return energy


if __name__ == "__main__":
    pass
