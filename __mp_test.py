import computing
import pipeline
import numpy as np
import time
import sunpy
import matplotlib.pyplot as plt
from astropy.io import fits


time_1 = "2011-02-13T02:00:00.000000000"
NOAA_AR = 11158

#magnetogram, bitmap = computing.download_map_and_harp(
#    time_1, time_1, NOAA_AR=NOAA_AR)


if __name__ == '__main__':
    magnetogram = [r"c:\Users\cosbo\sunpy\data\hmi.m_720s.20110213_020000_TAI.1.magnetogram.fits "]

    bitmap = [r"c:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110213_020000_TAI.bitmap.fits"]
    tic = time.perf_counter()
    #print(computing.mp_energy(bitmap, magnetogram, density=3, onlyactive=True, threads=10))
    toc = time.perf_counter()

    print(f'mp finished in {toc-tic:.2} s')
    print(computing.mp_energy(bitmap, magnetogram, density=2, onlyactive=True, threads=10, mode='opti'))
    #print(computing.single_bitmap_energy(bitmap, magnetogram, density=3, onlyactive=True))
    tic = time.perf_counter()

    print(f'cpp finished in {tic-toc:.2} s')