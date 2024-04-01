
import computing
import pandas as pd
import numpy as np
import cpp_module as cpp
import pipeline
import matplotlib.pyplot as plt
import plots
import os


N = 12  # количество временных точек
day = 13  # начальная дата
year = 2011
month = '02'
noaa_ar = 11158  # номер активной области
frequency = '8h'  # частота расчётов h - hour min - minute; данные приходят раз в 12 минут
start_time = '00:00:00'
home_path = os.path.expanduser("~") + '\\sunpy\\data'

density_array = np.linspace(10, 1, num=20-1)

if __name__ == "__main__":
    # __magnetogram_path, __bitmap_path = pipeline.download_map_and_harp(
    #        dates[0], dates[-1], NOAA_AR=noaa_ar)
    # быстрее выходит скачивать сразу много файлов
    for density in density_array:
        dates = pd.date_range(start=f"{year}-{month}-{day} {start_time}",
                      freq=frequency, periods=N).values
        np.savetxt(
            f'dates_{noaa_ar}_{day}_dens{density}_test.txt', dates, fmt='%s')
        try:
            energys = np.loadtxt(
                f'energys_{noaa_ar}_density={density}, day={day}_test.txt')
        except:
            energys = []
        print(energys)
        try:
            start = len(energys)
        except TypeError:
            start = 1
        dates = np.datetime_as_string(dates, unit='s')
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
                dts = np.loadtxt(
                    f'dates_{noaa_ar}_{day}_dens{density}_test.txt', dtype=np.datetime64)
                np.savetxt(f'dates_{noaa_ar}_{day}_dens{density}_test.txt', np.delete(
                    dts, i - failed_counter), fmt='%s')
                failed_counter += 1
                continue
            
            energy, x, y = computing.mp_energy(bitmap_path, magnetogram_path, density=density,
                                               onlyactive=True, mode='fineZ',
                                               follow_flux=False)
            print(f'{date} - {energy}')
            energys = np.append(energys, energy)
            #continue
            computing.alert_bot(f'{date} - {energy}, dens={density}')

            np.savetxt(
                f'energys_{noaa_ar}_density={density}, day={day}_test.txt', energys)

        print(f'{failed_counter} downloads failed')
        fig, ax = plots.config(xlabel='date', ylabel='flux, erg')
        dts = np.loadtxt(
            f'dates_{noaa_ar}_{day}_dens{density}_test.txt', dtype=np.datetime64)
        ax.plot(dts, energys, '-',
                label=f'density={density}')
    
    plt.show()