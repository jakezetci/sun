# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 21:22:29 2023

@author: cosbo
"""

import numpy as np

import os
import time
import telebot
import cpp_module as cpp
import multiprocessing as mp
from collections import namedtuple
import matplotlib.pyplot as plt
from collections.abc import Iterable
import matplotlib
import time
from requests.exceptions import ConnectionError

from coordinates import Coordinates, xyR2xyz
from lib import B_comp_map, Grid, create_grid, Magneticline
from field import dipolebetter
from plots import plotmap
from pipeline import bitmaps_to_points, arcsecs_to_radian

try:
    with open("callid.txt") as f:
        call_id = int(f.read())

    with open("TOKEN.txt") as ff:
        API_TOKEN = ff.read()
except FileNotFoundError:
    pass

shared_data_class = namedtuple(
    'shared_data_class', ['values', 'points', 'areas', 'y0'])

Grid_nt = namedtuple('Grid3D',
                     ['xyz', 'r', 'num', 'loc_x', 'loc_y', 'basic_volume'])


def __mponecube(vector):
    global values
    global points
    global areas
    global y0

    B = cpp.b_comp(vector, values, points, areas)

    if vector[1] == y0:
        return np.inner(B, B), vector
    return np.inner(B, B)


def __mponecube_adaptive(arg_array):
    global values
    global points
    global areas
    global y0

    x, y, z, volume = arg_array
    vector = np.array([x, y, z])
    B = cpp.b_comp(vector, values, points, areas)
    if vector[1] == y0:
        return np.inner(B, B)*volume, vector

    return np.inner(B, B) * volume


def __mp_init(shared_data):
    global values
    global points
    global areas
    global y0

    values = shared_data.values
    points = shared_data.points
    areas = shared_data.areas
    y0 = shared_data.y0
    pass


def mp_energy(bitmap_path, magnetogram_path, density=5,
              onlyactive=True, threads=mp.cpu_count(), mode='default',
              follow_flux=False):
    """calculates energy of a single active region provided
    ! utilising multiprocessing ! (about 6 times faster)

    Args:
        bitmap_path: bitmap file
        magnetogram_path: magnetogram file
        density: pixels in one cube. Defaults to 5.
        onlyactive: if True, only accounts for bitmap=34. if False, also takes bitmap=33. Defaults to True.
        threads: number of CPU threads to use. Defaults to all threads.
        mode: if 'opti' - boosts up the speed by lowering the number of cubes, 
                if 'track_loc' - also returns the position of the active region
    """
    if mode == 'plot':
        plt.close()
        matplotlib.rcParams.update({'font.size': 21})
        matplotlib.rcParams.update({'font.family': 'HSE Sans'})

        """
        fig, (ax1, ax2) = plt.subplots(2, 1,
                                       squeeze=True,
                                       figsize=(11,12),
                                       sharex=True,
                                       gridspec_kw={'height_ratios':[1,1],
                                                    })
        """
        fig, ax2 = plt.subplots(1, 1, figsize=(16, 9))
        ax2.set_xlabel('X, m')
        ax2.set_ylabel('Z, m')
        ax2.set_title('Проекция расчётного куба X-Z')
        ax1 = False

        # ax1, ax2 = False, False

    else:
        ax1, ax2 = False, False

    values, points, areas, hdrs, (cX, cY) = bitmaps_to_points(TIME=False, downloaded=True,
                                                              bitmaps=bitmap_path,
                                                              magnetogram=magnetogram_path,
                                                              onlyactive=onlyactive,
                                                              returnhdr=True,
                                                              plot=ax1)

    area = (696000000**2 / (2024*2024)) * np.pi

    # areas = np.full_like(areas, area)

    if follow_flux:
        return np.inner(np.abs(values), areas)
    grid = create_3Dgrid(hdrs[0], density, cX, cY, mode=mode)
    energy = 0.0
    a = 0
    valid_indeces = np.where(grid.r > a)
    valid_points = grid.xyz[valid_indeces]

    ys_unique = np.unique(valid_points[:, 1])
    y0 = ys_unique[np.argsort(ys_unique)[len(ys_unique)//2]]
    # N_y0 = np.count_nonzero(valid_points[:, 1] == y0)
    points_to_display, vectors = [], []
    # nan_start = 3620541

    shared_data = shared_data_class(values, points, areas, y0)

    pool = mp.Pool(processes=threads,
                   initializer=__mp_init, initargs=(shared_data, ))
    energy = 0.0
    func = __mponecube
    if mode == 'adaptive' or mode == 'fineZ':
        func = __mponecube_adaptive
        v = np.array([grid.basic_volume])
        valid_points = np.concatenate([grid.xyz, v.T], axis=1)
        # valid_points = valid_points[nan_start:]
        # r = grid.r[nan_start:]

    from tqdm import tqdm
    with tqdm(total=np.shape(valid_points)[0], maxinterval=0.1) as pbar:
        for res in pool.imap(func, valid_points):
            if isinstance(res, Iterable):
                vector = res[1]
                res = res[0]
                points_to_display.append(res)
                vectors.append([vector[0], vector[2]])
            if np.isnan(res):
                print('nan happened')
                pass
            else:
                energy = energy + res
            pbar.update(1)
        pool.close()
        pool.join()
        pbar.update()
        time.sleep(1.0)

    if mode == 'track_loc':
        return energy * grid.basic_volume/(8*np.pi), grid.loc_x, grid.loc_y
    elif mode == 'adaptive':
        return energy / (8*np.pi), grid.loc_x, grid.loc_y
    elif mode == 'fineZ':

        # alert_bot('посчитана картинка:', 'temp.png')
        return energy / (8*np.pi), grid.loc_x, grid.loc_y

    elif mode == 'plot':
        x, y = np.asarray(vectors).T
        x = -x
        tt = ax2.scatter(x, y, c=points_to_display, alpha=0.6,
                         norm=matplotlib.colors.LogNorm())
        cbar2 = plt.colorbar(tt, ax=ax2, location='right')
        cbar2.ax.set_ylabel("локальная плотность энергии",
                            size='medium')
        z_r = np.sqrt(696000000**2 - y0**2 - x**2)
        plt.plot(x, z_r, label='поверхность фотосферы', lw=2)
        plt.legend(loc='lower right')
        fig.savefig('temp.png')

        plt.show()

        return energy * grid.basic_volume/(8*np.pi), grid.loc_x, grid.loc_y

    else:
        return energy * grid.basic_volume/(8*np.pi), grid.loc_x, grid.loc_y


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

    grid = create_3Dgrid(hdrs[0], density, cX, cY)
    energy = 0.0

    a = 697000 * 1e3
    valid_indeces = np.where(grid.r > 0)
    valid_points = grid.xyz[valid_indeces]
    print(f'total values = {np.shape(valid_indeces)[1]}')
    from tqdm import tqdm

    with tqdm(total=np.shape(valid_indeces)[1], maxinterval=0.1) as pbar:
        for i, vector in enumerate(valid_points[::-1]):
            B = cpp.b_comp(vector, values, points, areas)
            if np.isnan(B).any():
                print(f'alert!! B is nan, i={i}')
            else:
                energy = energy + (np.linalg.norm(B)**2 * grid.basic_volume)
            pbar.update(1)
        pbar.update()
        time.sleep(2.0)
    return energy/(8*np.pi)


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
        print('No connection to the Internet to alert the bot')
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
    # B_map.save_pkl(name=f"checkpoint2 {name}", empty=False, vector=vector)
    # B_map.save_pkl(name=f"checkpoint1 {name}", empty=False, vector=vector)
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
    B_map.save_pkl(name=name)
    if alert is True:
        alert_bot("модельная сетка посчиталась...")
    if returnobj:
        return B_map


def create_3Dgrid(hdr, density, cX, cY, mode='default'):
    def tilt(a):

        return -a/np.sqrt(4096**2 - a**2)
    pxsizeX, pxsizeY = hdr["CDELT1"], hdr["CDELT2"]
    r_sun = hdr["RSUN_REF"]
    d_pixel = np.mean([pxsizeX, pxsizeY])
    dOBS = hdr["DSUN_OBS"]
    d_pixel = arcsecs_to_radian(d_pixel) * dOBS
    r_sun = r_sun + d_pixel/np.pi  # for good measure
    mapsizeX, mapsizeY = hdr['CRSIZE1'], hdr['CRSIZE2']
    ref1, ref2 = hdr["CRPIX1"], hdr["CRPIX2"]
    bitmapcenterX, bitmapcenterY = ref1 + mapsizeX/2, ref2 + mapsizeY/2

    extra_n = np.min((mapsizeX, mapsizeY))
    tilt_x, tilt_y = tilt(bitmapcenterX-cX), tilt(bitmapcenterY-cY)

    pixel_xs = np.linspace(start=ref1, stop=ref1 +
                           mapsizeX, num=int(mapsizeX//density))
    pixel_ys = np.linspace(start=ref2, stop=ref2 +
                           mapsizeY, num=int(mapsizeY//density))

    if tilt_x < 0:
        appendix = np.linspace(ref1+mapsizeX, ref1+mapsizeX+extra_n *
                               np.abs(tilt_x), num=int(np.abs(tilt_x)*extra_n//density))
        pixel_xs = np.append(pixel_xs, appendix)
    else:
        appendix = np.linspace(
            ref1, ref1-extra_n*np.abs(tilt_x), num=int(np.abs(tilt_x)*extra_n//density))
        pixel_xs = np.append(pixel_xs, appendix)

    if tilt_y < 0:

        appendix = np.linspace(ref2+mapsizeY, ref2+mapsizeY *
                               np.abs(tilt_y), num=int(np.abs(tilt_y)*extra_n//density))
        pixel_ys = np.append(pixel_ys, appendix)
    else:
        appendix = np.linspace(
            ref2, ref2-extra_n*np.abs(tilt_y), num=int(np.abs(tilt_y)*extra_n//density))
        pixel_ys = np.append(pixel_ys, appendix)

    loc_x, loc_y = -(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel
    xs_unique = -(pixel_xs-cX) * d_pixel
    ys_unique = -(pixel_ys-cY) * d_pixel

    # xs_unique = xs_unique + d_pixel/np.pi #shift to avoid nans
    # ys_unique = ys_unique + d_pixel/np.pi #shift

    if mode == 'default':
        z_num = np.min([mapsizeX, mapsizeY])
        z_size = z_num*d_pixel
        __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                              r_sun)
        zs_unique = np.linspace(z, z+z_size, num=int(z_num//density))
    elif mode == 'opti' or mode == 'track_loc' or mode == 'plot':
        z_num = np.min([mapsizeX, mapsizeY])
        z_size = z_num*d_pixel
        __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                              r_sun)
        zs_unique = np.linspace(z, z+z_size, num=int(z_num//(density)))
    elif mode == 'adaptive':
        z_num = np.max([mapsizeX, mapsizeY])
        z_size = z_num*d_pixel
        __x, __y, z = xyR2xyz(-(bitmapcenterX-cX)*d_pixel, -(bitmapcenterY-cY) * d_pixel,
                              r_sun)
        zs_unique = np.geomspace(z, z+z_size, num=int(z_num//(density*3)))
        xs, ys, zs = np.meshgrid(xs_unique, ys_unique, zs_unique)
        xs, ys, zs = xs.flatten(), ys.flatten(), zs.flatten()
        basic_volume = np.abs(np.roll(zs_unique, -1) - zs_unique)
        basic_volume[-1] = basic_volume[-2]
        basic_volume = basic_volume * \
            (xs_unique[1]-xs_unique[0]) * (ys_unique[1]-ys_unique[0])
    elif mode == 'fineZ':
        xs, ys = np.meshgrid(xs_unique, ys_unique)
        xs, ys = xs.flatten(), ys.flatten()
        z_num = int(np.max([mapsizeX, mapsizeY])//(density * 2) * 2)

        z_size = z_num*d_pixel
        xyz = np.zeros((z_num * len(xs), 3))
        basic_volume = np.zeros(z_num*len(xs))
        for i, (_x, _y) in enumerate(zip(xs, ys)):

            __x__, __y__, z = xyR2xyz(_x, _y, r_sun)
            # zs = np.array([np.linspace(z, z+z_size, num=int(z_num))])
            zs_small = np.linspace(z, z+z_size/3, num=int(z_num/2))
            zs_big = np.linspace(z+z_size/3, z+z_size, num=int(z_num/2))

            zs = np.array([np.hstack((zs_small, zs_big))])
            a = np.full((z_num, 2), [_x, _y])
            xyz[i*z_num:z_num*(i+1)] = np.concatenate((a, zs.T), axis=1)
            basic_volume[i*z_num:z_num*(i+1)] = np.abs(np.roll(zs, -1) - zs)
            basic_volume[z_num*(i+1)-1] = basic_volume[z_num*(i+1)-2]

        r = np.linalg.norm(xyz, axis=1)
        num = np.shape(r)[0]
        basic_volume = (xs_unique[1]-xs_unique[0]) * \
            (ys_unique[1]-ys_unique[0])*basic_volume

        grid = Grid_nt(xyz, r, num, loc_x, loc_y, basic_volume*1e6)

        return grid

    if mode != 'adaptive':
        xs, ys, zs = np.meshgrid(xs_unique, ys_unique, zs_unique)
        xs, ys, zs = xs.flatten(), ys.flatten(), zs.flatten()
        basic_volume = ((xs_unique[1]-xs_unique[0]) *
                        (ys_unique[1]-ys_unique[0])*(zs_unique[1]-zs_unique[0]))

    xyz = np.asarray([xs, ys, zs]).T
    r = np.linalg.norm(xyz, axis=1)
    num = np.shape(r)[0]

    if mode == 'track_loc' or 'plot':
        grid = Grid_nt(xyz, r, num, loc_x, loc_y, basic_volume*1e6)

        return grid
    else:
        grid = Grid_nt(xyz, r, num, loc_x, loc_y, basic_volume*1e6)
        return grid


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
            if np.linalg.norm(magline.points[-1]) < stoppoint:
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
            if np.linalg.norm(magline.points[-1]) < stoppoint:
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
    n_linesx=10,
    n_linesy=10,
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
        n_linesx=n_linesx,
        n_linesy=n_linesy,
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

    print(f'total values = {np.shape(xyz)[1]}')
    energy = cpp.energy(xyz, basic_volume, values,
                        points, areas, 40000)
    return energy


def compute_grid_energy(grid: Grid):
    energy = 0
    for value, area in zip(grid.valuesvector, grid.area):
        energy = energy + np.linalg.norm(value) * area
    return energy


if __name__ == "__main__":
    pass
