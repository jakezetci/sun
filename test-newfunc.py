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

'''2011-02-11T00:00:00.000000000
'''
"04.09.2017 - 08.09.2017"
N = 36  # количество временных точек
density = 3  # пикселей на один объёмный куб
day = 12  # начальная дата
year = 2011
month = '02'
noaa_ar = 11158  # номер активной области
frequency = '4h'  # частота расчётов h - hour min - minute; данные приходят раз в 12 минут
start_time = '00:00:00'


if __name__ == '__main__':

    dates = pd.date_range(start=f"{year}-{month}-{day} {start_time}",
                          freq=frequency, periods=N).values

    np.savetxt(
        f'dates_{noaa_ar}_{day}_dens{density}_curv.txt', dates, fmt='%s')
    try:
        energys = np.loadtxt(
            f'energys_{noaa_ar}_density={density}, day={day}_curv.txt')
    except:
        energys = []
    print(energys)
    start = 2
    dates = np.datetime_as_string(dates, unit='s')

    #__magnetogram_path, __bitmap_path = pipeline.download_map_and_harp(
    #        dates[0], dates[-1], NOAA_AR=noaa_ar)
    # быстрее выходит скачивать сразу много файлов

    failed_counter = 0
    xs, ys = [], []
    for i in range(start, N):
        date = dates[i]
        print(date)
        try:
            bitmap_path = [pipeline.file_name(date, 'hmi.mharp_720s', general_path=r"C:/Users/Professional/sunpy/data")]
            magnetogram_path = [pipeline.file_name(date, 'hmi.m_720s', general_path=r"C:/Users/Professional/sunpy/data")]
        except FileNotFoundError:
            dts = np.loadtxt(
                f'dates_{noaa_ar}_{day}_dens{density}_curv.txt', dtype=np.datetime64)
            np.savetxt(f'dates_{noaa_ar}_{day}_dens{density}_curv.txt', np.delete(
                dts, i - failed_counter), fmt='%s')
            failed_counter += 1
            continue

        energy, x, y = computing.mp_energy(bitmap_path, magnetogram_path, density=density,
                                           threads=6, onlyactive=True, mode='fineZ')
        computing.alert_bot(f'{date} - {energy:.3}')
        print(f'{date} - {energy}')

        energys = np.append(energys, energy)
        xs.append(x)
        ys.append(y)
        xy = np.vstack([xs, ys]).T

        np.savetxt(
            f'energys_{noaa_ar}_density={density}, day={day}_curv.txt', energys)
        np.savetxt(
            f'locs_{noaa_ar}_density={density}, day={day}_curv.txt', xy)

    print(f'{failed_counter} downloads failed')
    fig, ax = plots.config(xlabel='date', ylabel='enegry, erg')
    dts = np.loadtxt(
        f'dates_{noaa_ar}_{day}_dens{density}_curv.txt', dtype=np.datetime64)

    ax.plot(dts, energys, 'o'
            )
    plt.show()
