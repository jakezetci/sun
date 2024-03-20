# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 00:58:05 2023

@author: cosbo
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as tck
import numpy as np
import math


def config(
    xlimit=None,
    ylimit=None,
    xlabel="X",
    ylabel="Y",
    title="title",
    logscalex=False,
    logscaley=False,
    gapx=0.0,
    gapy=0.0,
    grid=True,
    figsize=(12, 6),
    dpi=100,
):
    """
    creates a nice standardized base for a matplotlib plot,
    works even if no parameters are provided
    Parameters
    ----------
    xlimit, ylimit, optional
        limits tuples in a format of (min, max)
    xlabel, ylabel, optional
        the label text for x/y-axis, by default "X"/"Y"
    title, optional
        figure title text, by default 'title'
    logscalex, logscaley, optional
        True if you you need x/y-axis logscaled, by default False
    gapx, gapy, optional
        gap between min/max or and the figure borders, by default 0.0
    grid, optional
        True if you need a grid, False if not, by default True
    figsize, optional
        figure size in inches, by default (16, 9)
    dpi, optional
        pixels per inch for a figure, by default 100

    Returns
    -------
        matplotlib.figure.Figure, matplotlib.axes.Axes
    """

    fig = plt.figure(layout="constrained")
    if xlimit is None is ylimit is None:
        ax = plt.axes()
    else:
        xlimit = [xlimit[0] - gapx, xlimit[1] + gapx]
        ylimit = [ylimit[0] - gapy, ylimit[1] + gapy]
        ax = plt.axes(xlim=xlimit, ylim=ylimit)

    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])
    fig.set_dpi(dpi)

    if grid:
        ax.grid(which="major", color="k")
        ax.grid(which="minor", color="gray", linestyle=":")
    elif grid == "minor":
        #ax.grid(which="major", color="k", linestyle=':')
        ax.grid(which="minor", color="gray", linestyle=":")

    COLOR = "k"
    mpl.rcParams["text.color"] = COLOR
    mpl.rcParams["axes.labelcolor"] = COLOR
    mpl.rcParams["xtick.color"] = COLOR
    mpl.rcParams["ytick.color"] = COLOR

    ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
    ax.tick_params(which="both", width=2)
    ax.tick_params(which="major", labelsize="large")
    ax.tick_params(which="minor", labelsize="medium", color="c")

    plt.xticks(fontname="Arial")
    plt.yticks(fontname="Arial")

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")

    ax.set_ylabel(ylabel, fontsize="xx-large", fontname="Arial")
    ax.set_xlabel(xlabel, fontsize="xx-large", fontname="Arial",)
    ax.set_title(title, fontsize="xx-large", fontname="Arial")
    if logscalex:
        ax.set_xscale("log")
    if logscaley:
        ax.set_yscale("log")
    return fig, ax


def subplots(
    nrows=1,
    ncols=1,
    xlimit=None,
    ylimit=None,
    xlabel=("xlabel"),
    ylabel=("ylabel"),
    title=("title"),
    logscale=False,
):
    """
    creates a nice standardized base for matplotlib subplots,
    works even if no parameters are provided

    Parameters
    ----------
    nrows, ncols, optional
        number or rows/columns of the subplot grid, by default 1
    xlimit, ylimit, optional
        2-d array of limits in a format of [(min1, max1), (min2, max2), ..]
    xlabel, ylabel, optional
        the label text for x/y-axis, by default "X"/"Y"
    title, optional
        figure title text, by default 'title'

    Returns
    -------
        matplotlib.figure.Figure, (nrows*ncols of matplotlib.axes.Axes)
    """
    if logscale is False:
        logscale = np.full(nrows * ncols, False)
    fig, axes = plt.subplots(nrows, ncols)

    fig.set_figwidth(16 * ncols)
    fig.set_figheight(6 * nrows)
    fig.set_dpi(50)

    COLOR = "#003087"
    mpl.rcParams["text.color"] = COLOR
    mpl.rcParams["axes.labelcolor"] = COLOR
    mpl.rcParams["xtick.color"] = COLOR
    mpl.rcParams["ytick.color"] = COLOR
    if xlimit is None:
        xlimit = np.full(nrows * ncols, None)
    if ylimit is None:
        ylimit = np.full(nrows * ncols, None)

    for ax, xl, yl, xx, yy, ttl, lg in zip(
        axes, xlimit, ylimit, xlabel, ylabel, title, logscale
    ):
        ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())

        ax.tick_params(which="both", width=2)
        ax.tick_params(which="major", labelsize=24, length=12)
        ax.tick_params(which="minor", labelsize=12, length=6, color="c")
        ax.set_ylabel(yy, fontsize=32, fontname="Arial")
        ax.set_xlabel(xx, fontsize=32, fontname="Arial")
        ax.set_title(ttl, fontsize=32, fontname="Arial")
        ax.grid(which="major", color="k")
        ax.grid(which="minor", color="gray", linestyle=":")

        if lg is True:
            ax.set_yscale("log")
    fig.tight_layout()
    return fig, axes


def sphere(
    ax, latlim=(-90, 90), lonlim=(-180, 180), n=5, r=696340, color="grey", lw=0.2
):
    """
    Generates a latitude - longitude grid on a subplot

    """
    lat_set = np.linspace(latlim[0], latlim[1], n)
    lon_set = np.linspace(lonlim[0], lonlim[1], n)
    lat_big = np.radians(90 - np.linspace(latlim[0], latlim[1]))
    lon_big = np.radians(np.linspace(lonlim[0], lonlim[1]))
    for one in lat_set:
        lats = np.full_like(lon_big, np.radians(90 - one))

        phi = lon_big
        theta = lats

        xx = np.sin(theta) * np.cos(phi) / (1 - np.cos(theta))
        yy = np.sin(theta) * np.sin(phi) / (1 - np.cos(theta))

        ax.plot(xx, yy, linewidth=0.2, color="grey")
    for one in lon_set:
        lons = np.full_like(lat_big, np.radians(one))

        phi = lons
        theta = lat_big

        xx = np.sin(theta) * np.cos(phi) / (1 - np.cos(theta))
        yy = np.sin(theta) * np.sin(phi) / (1 - np.cos(theta))

        ax.plot(xx, yy, linewidth=lw, color=color)


def disk(
    ax,
    latlim=(-90, 90),
    lonlim=(-90, 90),
    n_linesx=5,
    n_linesy=5,
    r=696340 * 1e3,
    color="dimgrey",
    lw=0.8,
    ignoretop=False,
    ignorecorners=0
):
    """
    plots a latitude-longitude grid of a sphere on axes
    (it looks like a portion of a disk, depending on your
    latitude and longitude limits)

    Parameters
    ----------
    ax
        Axes to which the disk is added
    latlim, lonlims, optional
        latitude/longitude limits of your projection, by default (-90, 90)
    n_lines, optional
        number of different lines by one axis,
        the actual number of lines would be 2*n_lines, by default 5
    r, optional
        radius of the grid being projected in meters, by default 696340*1e3
    color, optional
        color of grid lines, by default 'dimgrey'
    lw, optional
        line width of the grid, by default 0.8
    ignoretop, optional
        True if you need to remove annotations at top, by default False
    """

    lat_set = np.linspace(latlim[0], latlim[1], n_linesy)
    lon_set = np.linspace(lonlim[0], lonlim[1], n_linesx)
    lat_big = np.radians(90 - np.linspace(latlim[0], latlim[1]))
    lon_big = np.radians(np.linspace(lonlim[0], lonlim[1]))

    for i, one in enumerate(lon_set):
        lons = np.full_like(lat_big, np.radians(one))

        phi = lons
        theta = lat_big

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, "-.", linewidth=lw, color=color)
        arg_min, arg_max = np.argmin(yy), np.argmax(yy)
        if i == 0:
            rotationnum = math.atan((yy[1] - yy[0]) / (xx[1] - xx[0]))
            rotationnum = math.degrees(rotationnum) + 18
        elif i == n_linesx - 1:
            # ax.annotate(f"{one:.1f}", xy=(xx[arg_min], yy[arg_min]),
            #   fontsize="medium",
            #  color=color,
            # )
            pass

        elif i > (n_linesx - ignorecorners-1):
            pass
        else:
            ax.annotate(
                f"{one:.1f}",
                xy=(xx[arg_min], yy[arg_min]),
                fontsize="x-small",
                color=color,
                xytext=(0, -1),
                textcoords="offset fontsize",
            )
            if not ignoretop:
                ax.annotate(
                    f"{one:.1f}",
                    xy=(xx[arg_max], yy[arg_max]),
                    color=color,
                    fontsize="x-small",
                    xytext=(0, 0.5),
                    textcoords="offset fontsize",
                )
    for i, one in enumerate(lat_set):
        lats = np.full_like(lon_big, np.radians(90 - one))

        phi = lon_big
        theta = lats

        xx = r * np.sin(theta) * np.sin(phi)
        yy = r * np.cos(theta)

        ax.plot(xx, yy, "-.", linewidth=lw, color=color)
        arg_min, arg_max = np.argmin(xx), np.argmax(xx)
        if i == 0:
            ax.annotate(
                "longitudes",
                xy=(xx[arg_min], yy[arg_min]),
                fontsize="x-small",
                color=color,
                xytext=(0, -1),
                textcoords="offset fontsize",
            )
            ax.annotate(
                "latitudes",
                xy=(xx[arg_min], yy[arg_min]),
                color=color,
                fontsize="x-small",
                ha="right",
                rotation=90,
                xytext=(-0.5, 0),
                textcoords="offset fontsize",
            )

        else:
            ax.annotate(
                f"{one:.1f}",
                xy=(xx[arg_min], yy[arg_min]),
                fontsize="x-small",
                color=color,
                xytext=(-3, 0),
                textcoords="offset fontsize",
            )
            ax.annotate(
                f"{one:.1f}",
                xy=(xx[arg_max], yy[arg_max]),
                color=color,
                fontsize="x-small",
                xytext=(1, 0),
                textcoords="offset fontsize",
            )


def plotmap(
    B_map,
    mode=disk,
    limit=False,
    every=1,
    ms=0.5,
    alpha=1,
    n_linesx=1,
    n_linesy=1,
    lw=0.8,
    ignoretop=False,
    ignorecorners=3,
    **configargs,
):
    """
    plots values of the Grid

    Parameters
    ----------
    B_map: lib.Grid
        the Grid to plot
    mode, optional
        either disk or sphere, a mode of projection, by default disk
    limit, optional
        discards every value above limit if provided, by default False
    every, optional
        every-th value is printed (for stacked maps), by default 1
    ms, optional
        marker size of values plotted, by default 0.5
    alpha, optional
        transparency of markers, by default 1

    **disk parameters**

    n_lines, optional
        number of different lines by one axis, by default 1
    lw, optional
        line width of the grid, by default 0.8
    ignoretop, optional
        True if you need to remove annotations at top, by default False
    **configargs, optional
        parameters to be parsed to config()

    Returns
    -------
        matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    N = int(np.size(np.unique(B_map.lon)))
    if n_linesx > N or n_linesy > N:
        n_linesx = N
        n_linesy = N

    latlims = [B_map.lat.min(), B_map.lat.max()]
    lonlims = [B_map.lon.min(), B_map.lon.max()]
    fig, ax = config(**configargs)
    n = B_map.num
    if np.linalg.norm(B_map.valuesvector[0]) == 0:
        n = B_map.num
        xx = np.zeros(n // every + 1)
        yy = np.zeros(n // every + 1)
        sq = np.zeros(n // every + 1)
        for i, (cs, val) in enumerate(zip(B_map.coors_set, B_map.values)):
            if i % every == 0:
                x, y = cs.x, cs.y
                s = val
                xx[i // every], yy[i // every], sq[i // every] = x, y, s
        sc = ax.scatter(xx, yy, c=sq, s=ms, cmap="inferno", alpha=alpha)
        cbar = fig.colorbar(sc, location='bottom', shrink=0.6, aspect=50,
                            pad=0.0001)
        cbar.ax.set_xlabel("line-of-sight component of magnetic field, Mx/cm2",
                           size='medium')
        cbar.ax.tick_params(labelsize='small', direction='in')
    else:
        xx = np.zeros(n)
        yy = np.zeros(n)
        uu = np.zeros(n)
        vv = np.zeros(n)
        cc = np.zeros(n)
        for i, (cs, val) in enumerate(zip(B_map.coors_set, B_map.valuesvector)):
            x, y = cs.x, cs.y
            u, v, z = val
            if limit is not False:
                if abs(u) > limit or abs(v) > limit:
                    u, v = 0, 0
            xx[i], yy[i], uu[i], vv[i] = x, y, u, v
            cc[i] = math.sqrt(u**2 + v**2)
        ax.quiver(xx, yy, uu, vv, np.arctan2(vv, uu))
    mode(ax, latlims, lonlims, n_linesx=n_linesx, n_linesy=n_linesy,
         r=B_map.r, lw=lw, ignoretop=ignoretop, ignorecorners=ignorecorners)
    return fig, ax
