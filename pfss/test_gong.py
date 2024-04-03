# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 15:10:37 2024

@author: cosbo
"""


import pfsspy
import pfsspy.utils

import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a


from astropy.io import fits


def set_axes_lims(ax):
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)


if __name__ == "__main__":
    file = fits.open("mrzqs110215t1754c2107_313.fits")
    gong_map = sunpy.map.Map("mrzqs110215t1754c2107_313.fits")
    nrho = 35
    rss = 2.5
    pfss_in = pfsspy.Input(gong_map, nrho, rss)

    m = pfss_in.map
    fig = plt.figure()
    ax = plt.subplot(projection=m)
    m.plot()
    plt.colorbar()
    ax.set_title('Input field')
    set_axes_lims(ax)
    pfss_out = pfsspy.pfss(pfss_in)

    ss_br = pfss_out.source_surface_br
    # Create the figure and axes
    fig = plt.figure()
    ax = plt.subplot(projection=ss_br)

    # Plot the source surface map
    ss_br.plot()
    # Plot the polarity inversion line
    ax.plot_coord(pfss_out.source_surface_pils[0])
    # Plot formatting
    plt.colorbar()
    ax.set_title('Source surface magnetic field')
    set_axes_lims(ax)
