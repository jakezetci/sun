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
import os

'''2011-02-11T00:00:00.000000000
'''
"04.09.2017 - 08.09.2017"
alerts = True


N = 16  # количество временных точек
density = 4  # пикселей на один объёмный куб
day = 13  # начальная дата
year = 2009
month = '02'
noaa_ar = 11012  # номер активной области
frequency = '96min'  # частота расчётов h - hour min - minute; данные приходят раз в 12 минут
start_time = '00:00:00'
home_path = os.path.expanduser("~") + '\\sunpy\\data'

energy_fname = f'data/energy_{noaa_ar}_dens={density}_const.txt'
date_fname = f'data/dates_{noaa_ar}_dens{density}_const.txt'
loc_fname = f'data/locs_{noaa_ar}_dens{density}_const.txt'


if __name__ == '__main__':

    #dates_total = pd.date_range(start=f"{year}-{month}-{day} {start_time}",
    #                      freq=frequency, periods=N).values
    #
    #dates_old = np.loadtxt(
    #            f'dates_{noaa_ar}_{day}_dens{density}_const.txt', dtype=np.datetime64)
    #
    #dates = np.setdiff1d(dates_total, dates_old)
    dates = pd.date_range(start=f"{year}-{month}-{day} {start_time}",
                          freq=frequency, periods=N).values

    np.savetxt(
        date_fname, dates, fmt='%s')
    try:
        energys = np.loadtxt(
            energy_fname)
    except:
        energys = []
    print(energys)
    try:
        start = len(energys)
    except TypeError:
        start = 1

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
            bitmap_path = [pipeline.file_name(
                date, 'hmi.mharp_720s', general_path=home_path)]
            magnetogram_path = [pipeline.file_name(
                date, 'hmi.m_720s', general_path=home_path)]
        except FileNotFoundError:
            try:
                bitmap_path = [pipeline.file_name(
                    date, 'mdi.lostarp_96m', general_path=home_path)]
                magnetogram_path = [pipeline.file_name(
                    date, 'mdi.fd_m_96m_lev182', general_path=home_path)]
            except FileNotFoundError:
                dts = np.loadtxt(
                    date_fname, dtype=np.datetime64)
                np.savetxt(date_fname, np.delete(
                    dts, i - failed_counter), fmt='%s')
                failed_counter += 1
                continue

        energy, x, y = computing.mp_energy(bitmap_path, magnetogram_path, density=density,
                                           onlyactive=True, mode='fineZ', 
                                           follow_flux=False)
        
        
        print(f'{date} - {energy}')
        energys = np.append(energys, energy)
        #continue
        if alerts:
            computing.alert_bot(f'{date} - {energy}')

        xs.append(x)
        ys.append(y)
        xy = np.vstack([xs, ys]).T

        np.savetxt(
            energy_fname, energys)
        np.savetxt(
            loc_fname, xy)

    print(f'{failed_counter} downloads failed')
    fig, ax = plots.config(xlabel='date', ylabel='enegry, erg')
    dts = np.loadtxt(
        date_fname, dtype=np.datetime64)

    ax.plot(dts, energys, 'o'
            )
    plt.show()
