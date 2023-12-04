# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 16:32:42 2023

@author: cosbo
"""

from computing import comp_grid_points, compute_grid_energy, create_grid
from pipeline import bitmaps_to_points, compute_harp_MEnergy, download_map_and_harp
from plots import config
import numpy as np
import pandas as pd
from computing import alert_bot

"""
magnetogram = [
    "C:/Users/cosbo/sunpy/data/hmi.m_720s.20230808_001200_TAI.3.magnetogram.fits"
]
bitmaps = [
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9915.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9914.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9913.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9911.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9909.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9907.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9905.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9904.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9902.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9900.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9899.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9893.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9892.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9890.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9882.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9879.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9875.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9872.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9864.20230808_001200_TAI.bitmap.fits",
    "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9862.20230808_001200_TAI.bitmap.fits",
]
"""
latlim = (-30, 30)
lonlim = (-30, 30)
N = 0.25


dates = pd.date_range(start="2023-08-08 00:12:00", freq="6H", periods=4).values
energys_low = np.zeros(4)
energys_high = np.zeros(4)
for i, date in enumerate(dates):
    magnetogram, bitmaps = download_map_and_harp(date, date)
    energy_low = np.sum(
        compute_harp_MEnergy(
            date, date, downloaded=True, magnetogram=magnetogram, bitmaps=bitmaps
        )
    )
    values, points = bitmaps_to_points(
        date, downloaded=True, magnetogram=magnetogram, bitmaps=bitmaps
    )
    grid_empty = create_grid(latlim, lonlim, N)
    ts = pd.to_datetime(str(date))
    name = ts.strftime("%Y.%m.%d %I:%M%p")
    grid_high = comp_grid_points(
        grid_empty,points,values,
        checkpoints=15,
        timestamp=True,
        alert=True,
        name=name,
    )
    energy_high = compute_grid_energy(grid_high)
    energys_low[i] = energy_low
    energys_high[i] = energy_high
np.savetxt("HIGH.txt", energys_high)
np.savetxt("LOW.txt", energys_low)
fig, ax = config()
ax.plot(dates, energys_low, "o", label="low", color="pink")
ax.plot(dates, energys_high, "o", label="high", color="green")
ax.legend(loc="best")

if True:
    rand = np.random.randint(10, 99)
    fig.savefig("energy.png")
    alert_bot("вот картинка..", imagepath="energy.png")
