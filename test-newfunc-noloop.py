# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 16:06:03 2023

@author: cosbo
"""

import computing
import pandas as pd
import numpy as np

bitmap_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits']
magnetogram_path = [
    r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits']

energy = computing.single_bitmap_energy(bitmap_path, magnetogram_path, density=5,
                                        timestamp=10000, onlyactive=True)
