# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np
import cpp_module
import lib
import computing
import time
import pipeline
a = 696300 * 1000
r1 = np.asarray([1.0, 2.0, 3.0])*a
r2 = np.asarray([2.0, 4.0, 5.0])*a
debug = True


magnetogram_path = r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits'
bitmap_path = r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits'

values, points, areas = pipeline.bitmaps_to_points(False, downloaded=True, magnetogram=[magnetogram_path],
                                                   bitmaps=[bitmap_path])

time1 = time.perf_counter()
energy1 = computing.single_bitmap_energy([bitmap_path], [magnetogram_path], density=10,
                              timestamp=1000, onlyactive=True)
time2 = time.perf_counter()
energy2 = computing.fast_bitmap_energy([bitmap_path], [magnetogram_path], density=10,
                              timestamp=1000, onlyactive=True)
time3 = time.perf_counter()

time_c = time3-time2
time_def = time2-time1
print(f'выигрыш C++ в {time_def/time_c:.5} раза')
print(f'C++ - {time_c:.3} secs')
print(f'def - {time_def:.3} secs')
print(energy1)
print(energy2)
"""
time1 = time.perf_counter()
for i in range(100):
    cpp_module.green(*r1, *r2)
time2 = time.perf_counter()

for i in range(100):
    green_true = lib.GreenBl(r1, r2)
time3 = time.perf_counter()


for i in range(100):
    green_o = lib.Green_optimized(r1, r2)
time4 = time.perf_counter()


time_def = time3-time2
time_c = time2-time1
time_o = time4 - time3
print(f'C++ - {time_c:.3} secs')
print(f'green_def - {time_def:.3} secs')
print(f'green_optimized - {time_o:.3} secs')
print(f'выигрыш optimized в {time_def/time_o:.3} раза')
print(f'выигрыш C++ в {time_def/time_c:.3} раза')
"""

