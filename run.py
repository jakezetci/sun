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

alerts = True

density = 7 #chosen from tests :) although i would prefer 3

date_start, date_end = '2011-02-11', '2011-02-11'

__magnetogram_path, __bitmap_path, noaa_ars = pipeline.download_map_and_harp(
        date_start, date_end, returnAR=True)

print(noaa_ars)

for noaa, harp in noaa_ars.items():
    pathlib.Path().mkdir(parents=True, exist_ok=True)
