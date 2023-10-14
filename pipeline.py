# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:17:33 2023

@author: cosbo
"""

import aiapy
import sunpy
from sunpy.net import Fido, attrs as a
import astropy.units as u
import astropy.time
import textwrap

import drms

res = Fido.search(a.Time('2023.09.21_22:40:00_TAI', '2023.09.21_23:12:00_TAI'),
                           a.jsoc.Series('hmi.Mrmap_latlon_720s'))

print(res.show('TELESCOP', 'INSTRUME', 'T_OBS') )
# Create DRMS JSON client, use debug=True to see the query URLs
client = drms.Client()

# Query series info
series_info = client.info("hmi.Mrmap_latlon_720s")

# Print keyword info
print(f"Listing keywords for {series_info.name}:\n")
for keyword in sorted(series_info.keywords.index):
    keyword_info = series_info.keywords.loc[keyword]
    print(keyword)
    print(f"  type ....... {keyword_info.type} ")
    print(f"  recscope ... {keyword_info.recscope} ")
    print(f"  defval ..... {keyword_info.defval} ")
    print(f"  units ...... {keyword_info.units} ")
    print(f"  note ....... {keyword_info.note} ")

