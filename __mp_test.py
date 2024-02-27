import computing
import pipeline
import numpy as np
import time

time_1 = "2011-02-15 00:01:00"
NOAA_AR = 11158

#magnetogram, bitmap = computing.download_map_and_harp(
#    time, time, NOAA_AR=NOAA_AR)


if __name__ == '__main__':
    magnetogram = [r"C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits"]
    bitmap = [r"C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits"]
    tic = time.perf_counter()
    print(computing.mp_energy(bitmap, magnetogram, density=3, onlyactive=True, threads=10))
    toc = time.perf_counter()
    print(f'mp finished in {toc-tic:.2} s')
    print(computing.single_bitmap_energy(bitmap, magnetogram, density=3, onlyactive=True))
    tic = time.perf_counter()

    print(f'mp finished in {tic-toc:.2} s')