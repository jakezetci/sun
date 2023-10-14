# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 00:58:05 2023

@author: cosbo
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as tck
from matplotlib import cm
import numpy as np
import math


def config(xlimit=None, ylimit=None, xlabel="xlabel", ylabel="ylabel",
           title='title', logscalex=False, logscaley=False, gapx=0.0, gapy=0.0,
           grid=True, figsize=(16, 9), dpi=100):

    fig = plt.figure(layout='constrained')
    if xlimit is None is ylimit is None:
        ax = plt.axes()
    else:
        xlimit = [xlimit[0]-gapx, xlimit[1]+gapx]
        ylimit = [ylimit[0]-gapy, ylimit[1]+gapy]
        ax = plt.axes(
            xlim=xlimit,
            ylim=ylimit
            )

    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])
    fig.set_dpi(dpi)

    if grid:
        ax.grid(which='major', color='k')
        ax.grid(which='minor',
                color='gray',
                linestyle=':')

    COLOR = '#003087'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR

    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', labelsize='large')
    ax.tick_params(which='minor', labelsize='large', color='c')

    plt.xticks(fontname='HSE Slab')
    plt.yticks(fontname='HSE Slab')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_ylabel(ylabel, fontsize='xx-large', fontname='HSE Slab')
    ax.set_xlabel(xlabel, fontsize='xx-large', fontname='HSE Slab')
    ax.set_title(title, fontsize='x-large', fontname='HSE Slab')
    if logscalex:
        ax.set_xscale('log')
    if logscaley:
        ax.set_yscale('log')
    return fig, ax


def subplots(nrows=1, ncols=1, xlimit=None, ylimit=None, xlabel=("xlabel"),
             ylabel=("ylabel"), title=('title'), logscale=False):
    if logscale is False:
        logscale = np.full(nrows*ncols, False)
    fig, axes = plt.subplots(nrows, ncols)

    fig.set_figwidth(16*ncols)
    fig.set_figheight(6*nrows)
    fig.set_dpi(50)

    COLOR = '#003087'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR
    if xlimit is None:
        xlimit = np.full(nrows*ncols, None)
    if ylimit is None:
        ylimit = np.full(nrows*ncols, None)

    for ax, xl, yl, xx, yy, ttl, lg in zip(axes, xlimit, ylimit,
                                           xlabel, ylabel, title, logscale):
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

        ax.tick_params(which='both', width=2)
        ax.tick_params(which='major', labelsize=24, length=12)
        ax.tick_params(which='minor', labelsize=12, length=6, color='c')
        ax.set_ylabel(yy, fontsize=32, fontname='HSE Slab')
        ax.set_xlabel(xx, fontsize=32, fontname='HSE Slab')
        ax.set_title(ttl, fontsize=32, fontname='HSE Slab')
        ax.grid(which='major', color='k')
        ax.grid(which='minor',
                color='gray',
                linestyle=':')

        if lg is True:
            ax.set_yscale('log')
    fig.tight_layout()
    return fig, axes


def sphere(ax, latlim=(-90, 90), lonlim=(-180, 180), n=5,
           r=696340, color='grey', lw=0.2):
    """
    Generates a latitude - longitude grid on a subplot

    """
    lat_set = np.linspace(latlim[0], latlim[1], n)
    lon_set = np.linspace(lonlim[0], lonlim[1], n)
    lat_big = np.radians(90 - np.linspace(latlim[0], latlim[1]))
    lon_big = np.radians(np.linspace(lonlim[0], lonlim[1]))
    for one in lat_set:
        lats = np.full_like(lon_big, np.radians(90-one))

        phi = lon_big
        theta = lats

        xx = np.sin(theta) * np.cos(phi)/(1 - np.cos(theta))
        yy = np.sin(theta) * np.sin(phi)/(1 - np.cos(theta))

        ax.plot(xx, yy, linewidth=0.2, color='grey')
    for one in lon_set:
        lons = np.full_like(lat_big, np.radians(one))

        phi = lons
        theta = lat_big

        xx = np.sin(theta) * np.cos(phi)/(1 - np.cos(theta))
        yy = np.sin(theta) * np.sin(phi)/(1 - np.cos(theta))

        ax.plot(xx, yy, linewidth=lw, color=color)


def disk(ax, latlim=(-90, 90), lonlim=(-90, 90), n=5,
         r=696340, color='dimgrey', lw=0.8):
    lat_set = np.linspace(latlim[0], latlim[1], n)
    lon_set = np.linspace(lonlim[0], lonlim[1], n)
    lat_big = np.radians(90 - np.linspace(latlim[0], latlim[1]))
    lon_big = np.radians(np.linspace(lonlim[0], lonlim[1]))

    for i, one in enumerate(lon_set):
        lons = np.full_like(lat_big, np.radians(one))

        phi = lons
        theta = lat_big

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, '-.', linewidth=lw, color=color)
        arg_min, arg_max = np.argmin(yy), np.argmax(yy)
        if i == 0:
            rotationnum = math.atan((yy[2]-yy[0])/(xx[2] - xx[0]))
            rotationnum = math.degrees(rotationnum)+180
        elif i == n-1:
            ax.annotate(f'{one:.1f}', xy=(xx[arg_min], yy[arg_min]),
                        fontsize='large', color=color)

        else:
            ax.annotate(f'{one:.1f}', xy=(xx[arg_min], yy[arg_min]),
                        fontsize='large', color=color, xytext=(0,-1),
                        textcoords='offset fontsize')
            ax.annotate(f'{one:.1f}', xy=(xx[arg_max], yy[arg_max]),
                        color=color, fontsize='large', xytext=(0, 0.5),
                        textcoords='offset fontsize')
    for i, one in enumerate(lat_set):
        lats = np.full_like(lon_big, np.radians(90-one))

        phi = lon_big
        theta = lats

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, '-.', linewidth=lw, color=color)
        arg_min, arg_max = np.argmin(xx), np.argmax(xx)
        if i == 0:
            ax.annotate('longitudes', xy=(xx[arg_min], yy[arg_min]),
                        fontsize='x-large', color=color, xytext=(0, -1),
                        textcoords='offset fontsize')
            ax.annotate('latitudes', xy=(xx[arg_min], yy[arg_min]),
                        color=color,
                        fontsize='x-large', ha='right', rotation=rotationnum,
                        xytext=(-0.5, 0), textcoords='offset fontsize')

        else:
            ax.annotate(f'{one:.1f}', xy=(xx[arg_min], yy[arg_min]),
                        fontsize='large', color=color, xytext=(-2, 0),
                        textcoords='offset fontsize')
            ax.annotate(f'{one:.1f}', xy=(xx[arg_max], yy[arg_max]),
                        color=color, fontsize='large',  xytext=(1,0),
                        textcoords='offset fontsize')


def plotmap(B_map, mode=disk, limit=False, every=1, lines=1, ms=0.5, title='',
            alpha=1, lw=0.8):
    """
    Args:
        B_map (Grid): a grid with either values of valuesvector.

    Returns:
        A plt.figure.
    """
    N = int(np.size(np.unique(B_map.lon)))
    latlims = [B_map.lat.min(), B_map.lat.max()]
    lonlims = [B_map.lon.min(), B_map.lon.max()]
    fig, ax = config(title=title, grid=False)
    n = B_map.num
    if np.linalg.norm(B_map.valuesvector[0]) == 0:
        n = B_map.num
        xx = np.zeros(n//every + 1)
        yy = np.zeros(n//every + 1)
        sq = np.zeros(n//every + 1)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.values)):
            if i % every == 0:
                x, y = cs.x, cs.y
                s = val
                xx[i//every], yy[i//every], sq[i//every] = x, y, s
        ax.scatter(xx, yy, c=sq, s=ms, cmap='inferno', alpha=alpha)
    else:
        xx = np.zeros(n)
        yy = np.zeros(n)
        uu = np.zeros(n)
        vv = np.zeros(n)
        cc = np.zeros(n)
        for i, (cs, val) in enumerate(zip(B_map.coors_set,
                                          B_map.valuesvector)):
            x, y = cs.x, cs.y
            u, v, z = val
            if limit is not False:
                if abs(u) > limit or abs(v) > limit:
                    u, v = 0, 0
            xx[i], yy[i], uu[i], vv[i] = x, y, u, v
            cc[i] = math.sqrt(u**2 + v**2)
        ax.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))
    mode(ax, latlims, lonlims, n=lines, r=B_map.r, lw=lw)
    return fig, ax