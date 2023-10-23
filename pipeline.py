# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:17:33 2023

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

"""
res2 = Fido.search(a.Time('2023-08-08T00:00:00', '2023-08-08T00:20:00'),
                   a.jsoc.Series(series),
                   a.jsoc.Notify('rrzhdanov@edu.hse.ru'))

downloaded_files = Fido.fetch(res2)
print(downloaded_files)

maps = []
for file in downloaded_files:
    maps.append(fits.open(file))
    print(maps[-1].info())
data = maps[0][1].data
"""

Map = fits.open(
    'C:\\Users\\cosbo\\sunpy\\data\\hmi.m_720s.20230808_000000_TAI.3.magnetogram.fits')
data = Map[1].data
plt.imshow(data, cmap='inferno')
hdr = Map[1].header
kwords = list(hdr.keys())
CR_kwords = list(filter(lambda x: x.startswith("CR"), kwords))
for word in CR_kwords:
    print(f'{word} = {hdr[word]}')
