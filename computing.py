# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 21:22:29 2023

@author: cosbo
"""

import numpy as np

from astropy.io import fits
import astropy.units as u
import os
import time
import telebot
import cpp_module as cpp
import multiprocessing as mp
from collections import namedtuple




from coordinates import Coordinates, xyR2xyz
from lib import B_comp_map, Grid, create_grid, Magneticline, B_comp, Grid3D
from field import dipolebetter
from plots import sphere, disk, plotmap
from pipeline import download_map_and_harp, bitmaps_to_points, arcsecs_to_radian, bitmaps_to_points_slow

try:
    with open("callid.txt") as f:
        call_id = int(f.read())

    with open("TOKEN.txt") as ff:
        API_TOKEN = ff.read()
except FileNotFoundError:
    pass

shared_data_class = namedtuple(
    'shared_data_class', ['values', 'points', 'areas'])

def __mponecube(vector):
    global values
    global points
    global areas

    B = cpp.b_comp(vector, values, points, areas)

    
    return np.inner(B, B)

def __mp_init(shared_data):
    global values
    global points
    global areas
    
    values = shared_data.values
    points = shared_data.points
    areas = shared_data.areas
    pass


def mp_energy(bitmap_path, magnetogram_path, density=5,
           onlyactive=True, threads=mp.cpu_count(), mode='default',
           track_loc=False):
    """calculates energy of a single active region provided
    !utilising multiprocessing! (about 6 times faster)

    Args:
        bitmap_path: _description_
        magnetogram_path: _description_
        density: _description_. Defaults to 5.
        onlyactive: _description_. Defaults to True.
        threads: _description_. Defaults to mp.cpu_count().
    """
    values, points, areas, hdrs, (cX, cY) = bitmaps_to_points(TIME=False, downloaded=True,
                                                              bitmaps=bitmap_path,
                                                              magnetogram=magnetogram_path,
                                                              onlyactive=onlyactive,
                                                              returnhdr=True)

    xyz, r, basic_volume = create_3Dgrid(hdrs[0], density, cX, cY, mode=mode)
    energy = 0.0
    a = 697000 * 1e3
    valid_indeces = np.where(r > a)
    valid_points = xyz[valid_indeces] 
    shared_data = shared_data_class(values, points, areas)

    pool = mp.Pool(processes=threads, 
                   initializer=__mp_init, initargs=(shared_data, ))
    energy = 0.0
    from tqdm import tqdm

    with tqdm(total=np.shape(valid_indeces)[1], maxinterval=0.1) as pbar:
        for res in pool.imap(__mponecube, valid_points):
            if np.isnan(energy):
                pass
            else:
                energy = energy + res
            pbar.update(1)
        pool.close()
        pool.join()
        pbar.update()
        time.sleep(2.0)
    if track_loc:

        return energy * basic_volume/(8*np.pi), (cX, cY)
    else:
        return energy * basic_volume/(8*np.pi)

            


    

def single_bitmap_energy(bitmap_path, magnetogram_path, density=5,
                         onlyactive=True):
    """calculates energy of one single active region provided
    MAIN FUNCTION 

    Args:
        bitmap_path (_type_): _description_
        magnetogram_path (_type_): _description_
        density (int, optional): _description_. Defaults to 5.
        gap (float, optional): _description_. Defaults to 0.0.
        onlyactive (bool, optional): _description_. Defaults to True.
        timestamp (int, optional): _description_. Defaults to 100.
    """

    values, points, areas, hdrs, (cX, cY) = bitmaps_to_points(TIME=False, downloaded=True,
                                                              bitmaps=bitmap_path,
                                                              magnetogram=magnetogram_path,
                                                              onlyactive=onlyactive,
                                                              returnhdr=True)

    xyz, r, basic_volume = create_3Dgrid(hdrs[0], density, cX, cY)
    energy = 0.0
    tic = time.perf_counter()

    a = 697000 * 1e3
    valid_indeces = np.where(r > a)
    valid_points = xyz[valid_indeces]
    print(f'total values = {np.shape(valid_indeces)[1]}')
    from tqdm import tqdm

    with tqdm(total=np.shape(valid_indeces)[1], maxinterval=0.1) as pbar:
        for i, vector in enumerate(valid_points[::-1]):
            B = cpp.b_comp(vector, values, points, areas)
            if np.isnan(B).any():
                print(f'alert!! B is nan, i={i}')
            else:
                energy = energy + (np.linalg.norm(B)**2 * basic_volume)
            pbar.update(1)
        pbar.update()
        time.sleep(2.0)
    return energy/(8*np.pi)

def create_3Dgrid(hdr, density, cX, cY):

    pxsizeX, pxsizeY = hdr["CDELT1"], hdr["CDELT2"]
    r_sun = hdr["RSUN_REF"]
    d_pixel = np.mean([pxsizeX, pxsizeY])
    dOBS = hdr["DSUN_OBS"]
    d_pixel = arcsecs_to_radian(d_pixel) * dOBS
    mapsizeX, mapsizeY = hdr['CRSIZE1'], hdr['CRSIZE2']
    ref1, ref2 = hdr["CRPIX1"], hdr["CRPIX2"]
    bitmapcenterX, bitmapcenterY = ref1 + mapsizeX/2, ref2 + mapsizeY/2
    pixel_xs = np.linspace(start=ref1, stop=ref1 +
                           mapsizeX, num=int(mapsizeX/density))
    pixel_ys = np.linspace(start=ref2, stop=ref2 +
                           mapsizeY, num=int(mapsizeY/density))
    xs_unique = -(pixel_xs-cX) * d_pixel
    ys_unique = -(pixel_ys-cY) * d_pixel
    z_num = np.max([mapsizeX, mapsizeY])
    z_size = z_num*d_pixel
    __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                          r_sun)
    zs_unique = np.linspace(z, z+z_size, num=int(z_num//(density*5)))
    xs, ys, zs = np.meshgrid(xs_unique, ys_unique, zs_unique)
    xs, ys, zs = xs.flatten(), ys.flatten(), zs.flatten()
    basic_volume = ((xs_unique[1]-xs_unique[0]) *
                    (ys_unique[1]-ys_unique[0])*(zs_unique[1]-zs_unique[0]))
    grid = Grid3D(xs, ys, zs)
    return grid.xyz, grid.r, basic_volume * 1e6

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


def alert_bot(status_message, imagepath=False):
    try:
        bot = telebot.TeleBot(API_TOKEN)
        bot.send_message(call_id, status_message)
        if imagepath is not False:
            bot.send_photo(call_id, telebot.types.InputFile(imagepath))
    except ConnectionError:
        print('No connection to internet to alert bot')
        pass


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
    alert=False,
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
    #B_map.save_pkl(name=f"checkpoint2 {name}", empty=False, vector=vector)
    #B_map.save_pkl(name=f"checkpoint1 {name}", empty=False, vector=vector)
    r = B_map.r
    for i in range(B_map.progress, B_map.num):
        lat, lon = B_map.latlon[i]
        r1 = Coordinates(r, lat, lon, latlon=True)
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
        """if i % checkpointsnum == 0:
            os.remove(f"{folder}/checkpoint1 {name}.pkl")
            os.replace(
                f"{folder}/checkpoint2 {name}.pkl", f"{folder}/checkpoint1 {name}.pkl"
            )
            B_map.save_pkl(name=f"checkpoint2 {name}")
            """
    # B_map.save_pkl(name=name)
    if alert is True:
        alert_bot("модельная сетка посчиталась...")
    if returnobj:
        return B_map

def create_3Dgrid(hdr, density, cX, cY, mode='default'):



    pxsizeX, pxsizeY = hdr["CDELT1"], hdr["CDELT2"]
    r_sun = hdr["RSUN_REF"] 
    d_pixel = np.mean([pxsizeX, pxsizeY])
    dOBS = hdr["DSUN_OBS"]
    d_pixel = arcsecs_to_radian(d_pixel) * dOBS
    r_sun = r_sun + d_pixel # safe call
    mapsizeX, mapsizeY = hdr['CRSIZE1'], hdr['CRSIZE2']
    ref1, ref2 = hdr["CRPIX1"], hdr["CRPIX2"]
    bitmapcenterX, bitmapcenterY = ref1 + mapsizeX/2, ref2 + mapsizeY/2
    pixel_xs = np.linspace(start=ref1, stop=ref1 +
                           mapsizeX, num=int(mapsizeX/density))
    pixel_ys = np.linspace(start=ref2, stop=ref2 +
                           mapsizeY, num=int(mapsizeY/density))
    xs_unique = -(pixel_xs-cX) * d_pixel
    ys_unique = -(pixel_ys-cY) * d_pixel
    if mode == 'default':

        z_num = np.min([mapsizeX, mapsizeY])
        z_size = z_num*d_pixel
        __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                              r_sun)
        zs_unique = np.linspace(z, z+z_size, num=int(z_num//density))
    elif mode == 'opti':
        z_num = np.min([mapsizeX, mapsizeY])//3
        z_size = z_num*d_pixel
        __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                              r_sun)
        zs_unique = np.linspace(z, z+z_size, num=int(z_num//(density*3)))
    xs, ys, zs = np.meshgrid(xs_unique, ys_unique, zs_unique)
    xs, ys, zs = xs.flatten(), ys.flatten(), zs.flatten()
    basic_volume = ((xs_unique[1]-xs_unique[0]) *
                    (ys_unique[1]-ys_unique[0])*(zs_unique[1]-zs_unique[0]))
    grid = Grid3D(xs, ys, zs)
    return grid.xyz, grid.r, basic_volume * 1e6

def comp_grid(
    grid: Grid,
    B_map: Grid,
    vector=False,
    name=False,
    returnobj=False,
    checkpoints=721,
    timestamp=False,
    debug=False,
    alert=False,
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
    values, points, areas,
    name=False,
    returnobj=False,
    maxsteps=2000,
    timestamp=1000,
    alert=True,
    stoppoint=False,
    sign=+1
):
    def steps_comp():
        start = magline.progress
        tic = time.perf_counter()
        for i in range(start, maxsteps):
            magline.add_value_comp(values, points, areas, sign=sign)
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
            magline.add_value_comp(values, points, areas, sign=sign)
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
        n_lines=lines,
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
    alert=False,
):
    if name is False:
        name = crtname(name)
    checkpointsnum = int(checkpoints)
    tic = time.perf_counter()
    for i in range(grid.progress, grid.num):
        r1 = grid.coors_set[i]
        value = cpp.b_comp(r1, values, points)
        grid.set_value(value, index=i, vector=True)
        grid.progress1()
        if i % checkpointsnum == 0:
            if timestamp:
                toc = time.perf_counter()
                print(
                    f"values {i-checkpointsnum} - {i} done in {toc - tic:0.2f} seconds"
                )
                tic = time.perf_counter()
    grid.save_pkl(name=name)
    if alert is True:
        alert_bot("расчётная сетка посчиталась...")
    if returnobj:
        return grid


def smart_compute(time, onlyactive=True):
    def create_3Dgrid(header):
        return Grid3D()
    magnetogram, bitmap_paths = download_map_and_harp(time, time)

    for bitmap_path in bitmap_paths:
        bitmap = fits.open(bitmap_path)

        databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
        ref1, ref2 = hdrbitmap["CRPIX1"], hdrbitmap["CRPIX2"]
        active_indeces = np.argwhere(databitmap == 34)
        if onlyactive is False:
            quiet_indeces = np.argwhere(databitmap == 33)
            active_indeces = np.vstack([active_indeces, quiet_indeces])
        correction = np.full_like(active_indeces, [ref1, ref2])
        active_onmap = active_indeces + correction
        for xindex, yindex in active_onmap:
            """
            B = dataMap[xindex, yindex]
            x, y = -d_pixel * (xindex - centerX), -d_pixel * (yindex - centerY)
            points.append(xyR2xyz(x, y, r_sun))
            values.append(B)
            """



def cpp_bitmap_energy(bitmap_path, magnetogram_path, density=5, gap=0.0,
                      onlyactive=True):
    """actually slower pls do not use it

    Args:
        bitmap_path (_type_): _description_
        magnetogram_path (_type_): _description_
        density (int, optional): _description_. Defaults to 5.
        gap (float, optional): _description_. Defaults to 0.0.
        onlyactive (bool, optional): _description_. Defaults to True.
        timestamp (int, optional): _description_. Defaults to 100.
    """
    

    values, points, areas, hdrs, (cX, cY) = bitmaps_to_points(TIME=False, downloaded=True,
                                                              bitmaps=bitmap_path,
                                                              magnetogram=magnetogram_path,
                                                              onlyactive=onlyactive,
                                                              returnhdr=True)

    xyz, r, basic_volume = create_3Dgrid(hdrs[0], density, cX, cY)

    tic = time.perf_counter()
    print(f'total values = {grid.num}')
    energy = cpp.energy(xyz, basic_volume, values,
                        points, areas, 40000)
    return energy


def compute_grid_energy(grid: Grid):
    energy = 0
    for value, area in zip(grid.valuesvector, grid.area):
        energy = energy + np.linalg.norm(value) * area
    return energy


def model_grid3D():
    pass


if __name__ == "__main__":
    pass
