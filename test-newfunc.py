# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:51:01 2023

@author: cosbo
"""

import computing
import pandas as pd
import numpy as np
import cpp_module as cpp
import pipeline
import matplotlib.pyplot as plt
import plots
bitmap_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits']
magnetogram_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits']

density = 4.7
day = 13
N = 4
dates = pd.date_range(start=f"2011-02-{day} 12:00:00",
                      freq="3H", periods=N).values

np.savetxt(f'dates_11158_{day}_dens{density}_slow_3.txt', dates, fmt='%s')
try:
    energys = np.loadtxt(
        f'energys_11158_density={density}, day={day}_slow_3.txt')
except:
    energys = []
print(energys)
start = len([energys])
for i in range(start, N):
    date = dates[i]
    magnetogram_path, bitmap_path = computing.download_map_and_harp(
        date, date, NOAA_AR=11158)

    energy = computing.single_bitmap_energy(bitmap_path, magnetogram_path, density=density,
                                            timestamp=100000, onlyactive=True)
    energys = np.append(energys, energy)

    np.savetxt(
        f'energys_11158_density={density}, day={day}_slow_3.txt', energys)

    date_as_str = np.datetime_as_string(date)
    print(f'{date} - {energy}')

fig, ax = plots.config(xlabel='date', ylabel='enegry, erg')
ax.plot(dates, energy, 'o'
        )
