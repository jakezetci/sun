# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:53:14 2023

@author: cosbo
"""
import os.path
import matplotlib.pyplot as plt
from astropy.io import fits
import sunpy.visualization.colormaps as cm
import drms
import numpy as np
import sunpy
from sunpy.net import Fido, attrs as a
import astropy.units as u
import astropy.time
import textwrap
import sunpy.map.sources
from coordinates import xyz2ll
from lib import Grid


series = 'hmi.M_720s'
# Create DRMS JSON client, use debug=True to see the query URLs
client = drms.Client(email='rrzhdanov@edu.hse.ru')

# Query series info
series_info = client.info(series)

# Print keyword info

print(f"Listing keywords for {series_info.name}:\n")
for keyword in sorted(series_info.keywords.index):
    keyword_info = series_info.keywords.loc[keyword]
    print(keyword)
    print(f"  type ....... {keyword_info.type} ")
    print(f"  recscope ... {keyword_info.recscope} ")
    print(f"  defval ..... {keyword_info.defval} ")
    print(f"  units ...... {keyword_info.units} ")

print(client.pkeys(series))
print(series_info.segments.index.values)
print(f"  note ....... {keyword_info.note} ")

series_bitmap = "hmi.M_720s"
time = "2023-08-08 00:12:00"
res_M = Fido.search(
    a.Time(time, time),
    a.jsoc.Series(series_bitmap),
    a.jsoc.Notify("rrzhdanov@edu.hse.ru"),
)

downloaded_bitmaps = Fido.fetch(res_M).data

MAP = fits.open(downloaded_bitmaps[0])
data, hdr = MAP[-1].data, MAP[-1].header

print(np.shape(data))
kwords = list(hdr.keys())
for kword in sorted(kwords):
    print(f'{kword} : {hdr[kword]}')
