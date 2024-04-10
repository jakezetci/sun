"""
computes energy for multiple active regions at once
"""

import pathlib
import computing
import pandas as pd
import numpy as np
import cpp_module as cpp
import pipeline
import matplotlib.pyplot as plt
import plots
import os
from tqdm import tqdm


alerts = False

density = 5  # chosen from tests :) although i would prefer 3
path = "D:\\jsocdata"  # path for saving data, None for default
date_start, date_end = "2010-05-02", "2024-04-01"
time_start, time_end = "00:00:00", "22:00:00"
frequency = "2h"


if __name__ == "__main__":
    __magnetogram_path, __bitmap_path, noaa_ars = pipeline.download_map_and_harp(
        date_start, date_end, path=path, returnAR=True, limitlon=60.0
    )

    dates = pd.date_range(
        start=f"{date_start} {time_start}", freq=frequency, end=f"{date_end} {time_end}"
    ).values

    dates_string = np.datetime_as_string(dates, unit="s")

    N_total = len(noaa_ars) * len(dates)

    with tqdm(total=N_total, maxinterval=0.1, colour="green") as big_pbar:

        for noaa, harp in noaa_ars.items():

            dir = f"data/{noaa}"
            pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
            energy_fname = f"{dir}/energy_{noaa}_dens{density}_const.txt"
            date_fname = f"{dir}/dates_{noaa}_dens{density}_const.txt"
            np.savetxt(date_fname, dates, fmt="%s")
            loc_fname = f"{dir}/locs_{noaa}_dens{density}_const.txt"

            try:
                energys = np.loadtxt(energy_fname)
            except:
                energys = []
            try:
                start = len(energys)
            except TypeError:
                start = 1
            failed_counter = 0
            xs, ys = [], []
            for i, date in enumerate(dates_string[start:]):
                try:
                    bitmap_path = [
                        pipeline.file_name(
                            date, "hmi.mharp_720s." + str(harp), general_path=path
                        )
                    ]
                    magnetogram_path = [
                        pipeline.file_name(date, "hmi.m_720s", general_path=path)
                    ]
                except FileNotFoundError:
                    try:
                        bitmap_path = [
                            pipeline.file_name(
                                date, "mdi.lostarp_96m." + str(harp), general_path=path
                            )
                        ]
                        magnetogram_path = [
                            pipeline.file_name(
                                date, "mdi.fd_m_96m_lev182", general_path=path
                            )
                        ]
                    except FileNotFoundError:
                        dts = np.loadtxt(date_fname, dtype=np.datetime64)
                        np.savetxt(
                            date_fname, np.delete(dts, i - failed_counter), fmt="%s"
                        )
                        failed_counter += 1
                        big_pbar.update(1)
                        continue
                energy, x, y = computing.mp_energy(
                    bitmap_path,
                    magnetogram_path,
                    density,
                    onlyactive=False,
                    mode="fineZ",
                    follow_flux=False,
                )
                energys = np.append(energys, energy)
                if alerts:
                    computing.alert_bot(f"{date} - {energy}")

                xs.append(x)
                ys.append(y)
                xy = np.vstack([xs, ys]).T

                np.savetxt(energy_fname, energys)
                np.savetxt(loc_fname, xy)
                big_pbar.update(1)
