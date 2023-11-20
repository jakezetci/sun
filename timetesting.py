# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:37:25 2023

@author: cosbo
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from coordinates import Coordinates, ll2xyz
from lib import B_comp_map, Grid, create_grid
from field import dipolebetter
from plots import sphere, subplots, config
from plotting import plotmap
from computing import model_grid, model_magneticline, comp_magneticline, comp_grid
import pickle
import copy
from sklearn.metrics import mean_squared_error as mse
from sklearn.metrics import mean_absolute_percentage_error as mape
import time


def test_dipolelatlon(
    lattest,
    lontest,
    iters=10,
    dipoleR=600000,
    latlim=[-30, 30],
    lonlim=[-30, 30],
    N=2,
    M=1e12,
    Rrelative=np.array([10000, 4, -4]),
    relative=True,
    timestamps=10,
):
    NN = np.sqrt(iters / ((lattest[1] - lattest[0]) * (lontest[1] - lontest[0])))
    n = int(NN * (lattest[1] - lattest[0]))
    m = int(NN * (lontest[1] - lontest[0]))
    iters = n * m
    latstest = np.linspace(lattest[0], lattest[1], num=n)
    lonstest = np.linspace(lontest[0], lontest[1], num=m)
    hs = (np.pi / 180) * (lontest[1] - lontest[0]) / m
    latitudes = np.repeat(latstest, m)
    longitudes = np.tile(lonstest, n)
    L = np.asarray([latitudes, longitudes]).T
    B_c_array, B_d_array = np.zeros((iters, 3)), np.zeros((iters, 3))
    mse_array, mse_xyz_array = np.zeros(iters), np.zeros((iters, 3))
    for i in range(iters):
        if i % timestamps == 0:
            tic = time.perf_counter()
        dipolelat, dipolelon = L[i]
        mm = np.asarray(ll2xyz(dipolelat, dipolelon, 1)) * M
        pos = Coordinates(dipoleR, dipolelat, dipolelon, latlon=True)
        pos = pos.vector
        r_point = np.array([696340, dipolelat, dipolelon]) + Rrelative
        r_point = Coordinates(*r_point, latlon=True)
        grid = create_grid(
            [r_point.lat + latlim[0], r_point.lat + latlim[1]],
            [r_point.lon + lonlim[0], r_point.lon + lonlim[1]],
            N=N,
        )
        mapL = model_grid(
            grid,
            dipole=mm,
            dipolepos=pos,
            vector=False,
            name="Ltemptest",
            returnobj=True,
        )
        B_c = B_comp_map(r_point, mapL)
        B_d = dipolebetter(r_point, m=mm, rdipole=pos, returnxyz=True)
        B_c_array[i], B_d_array[i] = B_c, B_d
        if relative:
            mse_array[i] = mape(B_c, B_d)
            mse_xyz_array[i] = np.abs((B_c - B_d) / B_d)
        else:
            mse_array[i], mse_xyz_array[i] = mse(B_c, B_d), np.sqrt(np.abs(B_c - B_d))
        if i % timestamps == 0:
            toc = time.perf_counter()
            print(f"value {i} / {iters} done in {toc - tic:0.2f} seconds")
    index_min = np.argmin(mse_array)
    lat_min, lon_min = latitudes[index_min], longitudes[index_min]
    if relative:
        fig, (error_lat, error_lon) = subplots(
            2,
            1,
            xlabel=["dipole latitude", "dipole longitude"],
            ylabel=["relative error", "relative error"],
            title=["MAPE of latitudes", "MAPE of longitudes"],
        )
    else:
        fig, (error_lat, error_lon) = subplots(
            2,
            1,
            xlabel=["dipole latitude", "dipole longitude"],
            ylabel=["mean square error", "mean square error"],
            title=["MSE of latitudes", "MSE of longitudes"],
        )
    points_lat = np.where(longitudes == lon_min)
    error_lat.plot(
        latitudes[points_lat],
        mse_array[points_lat],
        "o",
        ms=8,
        label=f"longitude ={lon_min:.2f}",
    )
    error_lat.legend(loc="best", fontsize="x-large")
    points_lon = np.where(latitudes == lat_min)
    error_lon.plot(
        longitudes[points_lon],
        mse_array[points_lon],
        "o",
        ms=8,
        label=f"latitude ={lat_min:.2f}",
    )
    error_lon.legend(loc="best", fontsize="x-large")

    map_error = Grid(600000, latitudes, longitudes, hs=hs)
    for j, val in enumerate(mse_array):
        map_error.set_value(val, 0, 0, vector=False, index=j)
    plotmap(map_error, lines=1, ms=40)
    return (
        latitudes[index_min],
        longitudes[index_min],
        mse_array[index_min],
        mse_xyz_array[index_min],
        B_c_array[index_min],
        B_d_array[index_min],
        index_min,
    )


def test_pointlatlon(
    lattest,
    lontest,
    iters=10,
    dipoleR=600000,
    latlim=[-30, 30],
    lonlim=[-30, 30],
    N=2,
    M=1e12,
    dipolelat=10,
    dipolelon=17,
    relative=True,
    rRelative=10000,
    timestamps=10,
):
    NN = np.sqrt(iters / ((lattest[1] - lattest[0]) * (lontest[1] - lontest[0])))
    n = int(NN * (lattest[1] - lattest[0]))
    m = int(NN * (lontest[1] - lontest[0]))
    iters = n * m
    latstest = np.linspace(lattest[0], lattest[1], num=n)
    lonstest = np.linspace(lontest[0], lontest[1], num=m)
    hs = (np.pi / 180) * (lontest[1] - lontest[0]) / m
    latitudes = np.repeat(latstest, m)
    longitudes = np.tile(lonstest, n)
    L = np.asarray([latitudes, longitudes]).T
    B_c_array, B_d_array = np.zeros((iters, 3)), np.zeros((iters, 3))
    mse_array, mse_xyz_array = np.zeros(iters), np.zeros((iters, 3))
    for i in range(iters):
        if i % timestamps == 0:
            tic = time.perf_counter()
        Rrelative = np.array([rRelative, *L[i]])
        mm = np.asarray(ll2xyz(dipolelat, dipolelon, 1)) * M
        pos = Coordinates(dipoleR, dipolelat, dipolelon, latlon=True)
        pos = pos.vector
        r_point = np.array([696340, dipolelat, dipolelon]) + Rrelative
        r_point = Coordinates(*r_point, latlon=True)
        grid = create_grid(
            [r_point.lat + latlim[0], r_point.lat + latlim[1]],
            [r_point.lon + lonlim[0], r_point.lon + lonlim[1]],
            N=N,
        )
        mapL = model_grid(
            grid,
            dipole=mm,
            dipolepos=pos,
            vector=False,
            name="Ltemptest",
            returnobj=True,
        )
        B_c = B_comp_map(r_point, mapL)
        B_d = dipolebetter(r_point, m=mm, rdipole=pos, returnxyz=True)
        B_c_array[i], B_d_array[i] = B_c, B_d
        if relative:
            mse_array[i] = mape(B_c, B_d)
            mse_xyz_array[i] = np.abs((B_c - B_d) / B_d)
        else:
            mse_array[i], mse_xyz_array[i] = mse(B_c, B_d), np.sqrt(np.abs(B_c - B_d))
        if i % timestamps == 0:
            toc = time.perf_counter()
            print(f"value {i} / {iters} done in {toc - tic:0.2f} seconds")
    index_min = np.argmin(mse_array)
    lat_min, lon_min = latitudes[index_min], longitudes[index_min]
    if relative:
        fig, (error_lat, error_lon) = subplots(
            2,
            1,
            xlabel=["point latitude", "point longitude"],
            ylabel=["relative error", "relative error"],
            title=["MAPE of latitudes", "MAPE of longitudes"],
        )
    else:
        fig, (error_lat, error_lon) = subplots(
            2,
            1,
            xlabel=["point latitude", "point longitude"],
            ylabel=["mean square error", "mean square error"],
            title=["MSE of latitudes", "MSE of longitudes"],
        )
    points_lat = np.where(longitudes == lon_min)
    error_lat.plot(
        latitudes[points_lat],
        mse_array[points_lat],
        "o",
        ms=8,
        label=f"longitude ={lon_min:.2f}",
    )
    points_lon = np.where(latitudes == lat_min)
    error_lon.plot(
        longitudes[points_lon],
        mse_array[points_lon],
        "o",
        ms=8,
        label=f"latitude ={lat_min:.2f}",
    )
    error_lat.legend(loc="upper right", fontsize="x-large")

    error_lon.legend(loc="upper right", fontsize="x-large")
    map_error = Grid(600000, latitudes, longitudes, hs=hs)
    for j, val in enumerate(mse_array):
        map_error.set_value(val, 0, 0, vector=False, index=j)
    plotmap(map_error, lines=1, ms=40)
    return (
        latitudes[index_min],
        longitudes[index_min],
        mse_array[index_min],
        mse_xyz_array[index_min],
        B_c_array[index_min],
        B_d_array[index_min],
        index_min,
    )


def test_borders(
    sizelims,
    iters=100,
    dipoleR=600000,
    N=2,
    M=1e12,
    Rrelative=np.array([10000, -12, 12]),
    relative=True,
    dipolelat=10,
    dipolelon=17,
    timestamps=10,
):
    size_array = np.linspace(sizelims[0], sizelims[1], num=iters)
    B_c_array, B_d_array = np.zeros((iters, 3)), np.zeros((iters, 3))
    mse_array, mse_xyz_array = np.zeros(iters), np.zeros((iters, 3))
    for i in range(iters):
        if i % timestamps == 0:
            tic = time.perf_counter()
        sz = size_array[i]
        latlim, lonlim = [-sz, +sz], [-sz, sz]
        mm = np.asarray(ll2xyz(dipolelat, dipolelon, 1)) * M
        pos = Coordinates(dipoleR, dipolelat, dipolelon, latlon=True)
        pos = pos.vector
        r_point = np.array([696340, dipolelat, dipolelon]) + Rrelative
        r_point = Coordinates(*r_point, latlon=True)
        grid = create_grid(
            [r_point.lat + latlim[0], r_point.lat + latlim[1]],
            [r_point.lon + lonlim[0], r_point.lon + lonlim[1]],
            N=N,
        )
        mapL = model_grid(
            grid,
            dipole=mm,
            dipolepos=pos,
            vector=False,
            name="Ltemptest",
            returnobj=True,
        )
        B_c = B_comp_map(r_point, mapL)
        B_d = dipolebetter(r_point, m=mm, rdipole=pos, returnxyz=True)
        B_c_array[i], B_d_array[i] = B_c, B_d
        if relative:
            mse_array[i] = mape(B_c, B_d)
            mse_xyz_array[i] = np.abs((B_c - B_d) / B_d)
        else:
            mse_array[i], mse_xyz_array[i] = mse(B_c, B_d), np.sqrt(np.abs(B_c - B_d))
        if i % timestamps == 0:
            toc = time.perf_counter()
            print(f"value {i} / {iters} done in {toc - tic:0.2f} seconds")
    index_min = np.argmin(mse_array)
    sizemin = size_array[index_min]
    if relative:
        fig, error_size = config(
            None,
            None,
            "grid size in degrees",
            "relative error",
            "MAPE of grid size",
            logscaley=True,
            gapx=1.0,
            gapy=0.5,
        )
    else:
        fig, error_size = config(
            None,
            None,
            "grid size in degrees",
            "absolute error",
            "MSE of grid size",
            logscaley=True,
            gapx=1.0,
            gapy=0.5,
        )
    error_size.plot(size_array, mse_array, "o", ms=6, label=f"size min ={sizemin:.2f}")
    error_size.legend(loc="best", fontsize="x-large")

    return (
        size_array[index_min],
        mse_array[index_min],
        mse_xyz_array[index_min],
        B_c_array[index_min],
        B_d_array[index_min],
        index_min,
    )


def test_density(
    Nlims,
    iters=10,
    dipoleR=600000,
    latlim=[-10, 10],
    lonlim=[-10, 10],
    M=1e12,
    Rrelative=np.array([10000, 4, -4]),
    relative=True,
    dipolelat=10,
    dipolelon=17,
    timestamps=10,
):
    N_array = np.linspace(Nlims[0], Nlims[1], num=iters)
    B_c_array, B_d_array = np.zeros((iters, 3)), np.zeros((iters, 3))
    mse_array, mse_xyz_array = np.zeros(iters), np.zeros((iters, 3))
    for i in range(iters):
        if i % timestamps == 0:
            tic = time.perf_counter()
        NN = N_array[i]
        mm = np.asarray(ll2xyz(dipolelat, dipolelon, 1)) * M
        pos = Coordinates(dipoleR, dipolelat, dipolelon, latlon=True)
        pos = pos.vector
        r_point = np.array([696340, dipolelat, dipolelon]) + Rrelative
        r_point = Coordinates(*r_point, latlon=True)
        grid = create_grid(
            [r_point.lat + latlim[0], r_point.lat + latlim[1]],
            [r_point.lon + lonlim[0], r_point.lon + lonlim[1]],
            N=NN,
        )
        mapL = model_grid(
            grid,
            dipole=mm,
            dipolepos=pos,
            vector=False,
            name="Ltemptest",
            returnobj=True,
        )
        B_c = B_comp_map(r_point, mapL)
        B_d = dipolebetter(r_point, m=mm, rdipole=pos, returnxyz=True)
        B_c_array[i], B_d_array[i] = B_c, B_d
        if relative:
            mse_array[i] = mape(B_c, B_d)
            mse_xyz_array[i] = np.abs((B_c - B_d) / B_d)
        else:
            mse_array[i], mse_xyz_array[i] = mse(B_c, B_d), np.sqrt(np.abs(B_c - B_d))
        if i % timestamps == 0:
            toc = time.perf_counter()
            print(f"value {i} / {iters} done in {toc - tic:0.2f} seconds")
    index_min = np.argmin(mse_array)
    Nmin = N_array[index_min]
    if relative:
        fig, error_N = config(
            None,
            None,
            "grid N in degree",
            "relative error (log scale)",
            "MAPE of grid N",
            logscalex=False,
            logscaley=True,
            gapx=1.0,
            gapy=0.1,
        )
    else:
        fig, error_N = config(
            None,
            None,
            "grid N in degrees",
            "absolute error (log scale) ",
            "MSE of grid N",
            logscalex=False,
            logscaley=True,
            gapx=1.0,
            gapy=0.1,
        )
    error_N.plot(N_array, mse_array, "o", ms=6, label=f"N min ={Nmin:.2f}")
    error_N.legend(loc="best", fontsize="x-large")
    return (
        N_array[index_min],
        mse_array[index_min],
        mse_xyz_array[index_min],
        B_c_array[index_min],
        B_d_array[index_min],
        index_min,
    )


if __name__ == "__main__":
    latlimits = [-50, 50]
    lonlimits = [-50, 50]
    r1, r2, r3, r4, r5, r6, r7 = test_dipolelatlon(
        latlimits,
        lonlimits,
        iters=500,
        latlim=[-10, 10],
        lonlim=[-10, 10],
        relative=True,
    )

    l1, l2 = [-20, 20], [-20, 20]
    rr1, rr2, rr3, rr4, rr5, rr6, rr7 = test_pointlatlon(
        l1,
        l2,
        iters=500,
        latlim=[-10, 10],
        lonlim=[-10, 10],
        relative=True,
        dipolelat=r1,
        dipolelon=r2,
    )

    sizelims = [0.5, 45]
    rrr1, rrr2, rrr3, rrr4, rrr5, rrr6 = test_borders(
        sizelims,
        iters=100,
        Rrelative=np.array([10000, rr1, rr2]),
        relative=True,
        dipolelat=r1,
        dipolelon=r2,
        timestamps=10,
    )

    perfsize = [-rrr1, +rrr1]
    N_lims = [0.2, 10]
    t1, t2, t3, t4, t5, t6 = test_density(
        N_lims,
        iters=50,
        Rrelative=np.array([10000, rr1, rr2]),
        relative=True,
        dipolelat=r1,
        dipolelon=r2,
        latlim=perfsize,
        lonlim=perfsize,
        timestamps=10,
    )

    print(
        f"best parameters: dipole coors: {r1:0.2f, r2:0.2f},",
        f"point relative coors: {rr1:0.2f, rr2:0.2f},",
        f"relative borders:0 [{-rrr1, rrr2}]",
        f"grid density: {t1}",
    )
