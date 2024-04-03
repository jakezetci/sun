
import computing
import pandas as pd
import numpy as np
import cpp_module as cpp
import pipeline
import matplotlib.pyplot as plt
import plots
import os
import time
N = 1  # количество временных точек
day = 14  # начальная дата
year = 2011
month = '02'
noaa_ar = 11158  # номер активной области
frequency = '16h'  # частота расчётов h - hour min - minute; данные приходят раз в 12 минут
start_time = '06:00:00'
home_path = os.path.expanduser("~") + '\\sunpy\\data'

if __name__ == "__main__":
    time_spent = []
    density_array = np.hstack((np.linspace(1, 10, num=10-1), np.linspace(10.125, 20.125, 9)))[1:]

    dates = pd.date_range(start=f"{year}-{month}-{day} {start_time}",
                      freq=frequency, periods=2)
    dates = np.datetime_as_string(dates, unit='s')

    date = dates[0]
    print(date)
    for density in density_array:
        bitmap_path = [pipeline.file_name(
            date, 'hmi.mharp_720s', general_path=home_path)]
        magnetogram_path = [pipeline.file_name(
            date, 'hmi.m_720s', general_path=home_path)]
        tic = time.perf_counter()

        energy, x, y = computing.mp_energy(bitmap_path, magnetogram_path, density=density,
                                               onlyactive=True, mode='fineZ',
                                               follow_flux=False)
        toc = time.perf_counter()
        time_spent.append((toc-tic))
        computing.alert_bot(f'time spent for density {density} = {(toc-tic)}')
        np.savetxt('times_of_density2.txt', time_spent)
