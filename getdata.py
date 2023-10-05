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

response = Fido.search(a.Time('2023.09.21_22:40:00_TAI', '2023.09.21_23:12:00_TAI'),
                       a.jsoc.Series('hmi.Mrmap_latlon_720s'))