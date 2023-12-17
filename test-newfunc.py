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
bitmap_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits']
magnetogram_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits']


dates = pd.date_range(start="2011-02-13 00:00:00",
                      freq="4H", periods=5).values

np.savetxt('dates_11158_13_dens2.txt', dates, fmt='%s')
energys = []

for date in dates:
    magnetogram_path, bitmap_path = computing.download_map_and_harp(
        date, date, NOAA_AR=11158)

    energy = computing.single_bitmap_energy(bitmap_path, magnetogram_path, density=2,
                                            timestamp=500000, onlyactive=True)
    energys.append(energy)

    np.savetxt('energys_0213_dens2_all.txt', energys)

    date_as_str = np.datetime_as_string(date)
    print(f'{date_as_str} - {energy}')
