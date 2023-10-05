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


def config(xlimit=None, ylimit=None, xlabel="xlabel", ylabel="ylabel",
           title='title', logscale=False):

    fig = plt.figure(layout='constrained')
    if xlimit is None is ylimit == None:
        ax = plt.axes()
    else:
        ax = plt.axes(
            xlim=xlimit,
            ylim=ylimit
            )

    fig.set_figwidth(16)
    fig.set_figheight(9)
    fig.set_dpi(50)

    ax.grid(which='major', color = 'k')
    ax.grid(which='minor',
            color = 'gray',
            linestyle = ':')

    COLOR = '#003087'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR

    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', labelsize=24, length=12)
    ax.tick_params(which='minor', labelsize=12, length=6, color='c')

    plt.xticks(fontname='HSE Slab')
    plt.yticks(fontname='HSE Slab')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_ylabel(ylabel, fontsize=32, fontname='HSE Slab')
    ax.set_xlabel(xlabel, fontsize=32, fontname='HSE Slab')
    ax.set_title(title, fontsize=32, fontname='HSE Slab')

    return fig, ax


def subplots(nrows=1, ncols=1, xlimit=None, ylimit=None, xlabel=("xlabel"),
             ylabel=("ylabel"), title=('title'), logscale=False):
    fig, axes = plt.subplots(nrows, ncols)
    fig.set_figwidth(8*ncols)
    fig.set_figheight(9*nrows)
    fig.set_dpi(50)

    COLOR = '#003087'
    mpl.rcParams['text.color'] = COLOR
    mpl.rcParams['axes.labelcolor'] = COLOR
    mpl.rcParams['xtick.color'] = COLOR
    mpl.rcParams['ytick.color'] = COLOR


    for ax, xl, yl, xx, yy, ttl in zip(axes, xlimit, ylimit,
                                       xlabel, ylabel, title):
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

        ax.tick_params(which='both', width=2)
        ax.tick_params(which='major', labelsize=24, length=12)
        ax.tick_params(which='minor', labelsize=12, length=6, color='c')
        ax.set_ylabel(yy, fontsize=32, fontname='HSE Slab')
        ax.set_xlabel(xx, fontsize=32, fontname='HSE Slab')
        ax.set_title(ttl, fontsize=32, fontname='HSE Slab')
        ax.set_xlim(xl)
        ax.set_ylim(yl)
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
         r=696340, color='grey', lw=0.2):
    lat_set = np.linspace(latlim[0], latlim[1], n)
    lon_set = np.linspace(lonlim[0], lonlim[1], n)
    lat_big = np.radians(90 - np.linspace(latlim[0], latlim[1]))
    lon_big = np.radians(np.linspace(lonlim[0], lonlim[1]))
    for one in lat_set:
        lats = np.full_like(lon_big, np.radians(90-one))

        phi = lon_big
        theta = lats

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, linewidth=0.2, color=color)

    for one in lon_set:
        lons = np.full_like(lat_big, np.radians(one))

        phi = lons
        theta = lat_big

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, linewidth=lw, color=color)
