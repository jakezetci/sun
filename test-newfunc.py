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

N = 72 # количество временных точек
density = 3
day = 12

noaa_ar = 11158
"""
missing - 2011-02-12T12:00:00
2011-02-14T00:00:00
2011-02-14T01:00:00
2011-02-14T00:00:00
2011-02-14T17:00:00

"""
if __name__ == '__main__':

    dates = pd.date_range(start=f"2011-02-{day} 12:00:00",
                          freq="1h", periods=N).values

    np.savetxt(f'dates_{noaa_ar}_{day}_dens{density}.txt', dates, fmt='%s')
    try:
        energys = np.loadtxt(
            f'energys_{noaa_ar}_density={density}, day={day}.txt')
    except:
        energys = []
    print(energys)
    start = len(energys)
    dates = np.datetime_as_string(dates, unit='s')

    __magnetogram_path, __bitmap_path = pipeline.download_map_and_harp(
            dates[0], dates[-1], NOAA_AR=noaa_ar)
    #быстрее выходит скачивать сразу много файлов

    failed_counter = 0
    for i in range(start+1, N):
        date = dates[i]
        print(date)
        try:
            bitmap_path = [pipeline.file_name(date, 'hmi.mharp_720s')]
            magnetogram_path = [pipeline.file_name(date, 'hmi.m_720s')]
        except FileNotFoundError:
            dts = np.loadtxt(f'dates_{noaa_ar}_{day}_dens{density}.txt')
            np.savetxt(f'dates_{noaa_ar}_{day}_dens{density}.txt', np.delete(dts, i - failed_counter))
            failed_counter += 1
            continue

        energy, x, y = computing.mp_energy(bitmap_path, magnetogram_path, density=density,
                                     onlyactive=True, threads=10, mode='track_loc')
        print(x,y)
        energys = np.append(energys, energy)

        np.savetxt(
            f'energys_{noaa_ar}_density={density}, day={day}.txt', energys)

        print(f'{date} - {energy}')
    print(f'{failed_counter} downloads failed')
    fig, ax = plots.config(xlabel='date', ylabel='enegry, erg')
    ax.plot(dates, energy, 'o'
            )
    plt.show()