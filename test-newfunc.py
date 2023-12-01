# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:51:01 2023

@author: cosbo
"""

from computing import single_bitmap_energy

bitmap_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits']
magnetogram_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits']

energy = single_bitmap_energy(bitmap_path, magnetogram_path, density=30,
                              timestamp=5, onlyactive=True)
print(energy)
