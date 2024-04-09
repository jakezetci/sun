from astropy.io import fits
import drms
import numpy as np
import sunpy
from sunpy.net import Fido, attrs as a
import astropy.units as u
import sunpy.map.sources
import urllib
import os
from collections.abc import Iterable
import matplotlib.pyplot as plt


from coordinates import xyz2ll, xyR2xyz, Coordinates
from lib import Grid, B_comp

import time


if __name__ == "__main__":
    time_s = '2010.10.27_20:48:00_TAI'
    series_tarp = 'mdi.lostarp_96m'
    res = Fido.search(a.Time(time_s, time_s),
                      a.jsoc.Series(series_tarp),
                      a.jsoc.Keyword('NOAA_AR') == 11117,
                      a.jsoc.Segment("bitmap"),
                      a.jsoc.Notify("rrzhdanov@edu.hse.ru"))

    down = Fido.fetch(res).data
    MAP_mask = fits.open(down[0])
    series_fd = 'mdi.fd_M_96m_lev182'
    res_2 = Fido.search(a.Time(time_s, time_s),
                        a.jsoc.Series(series_fd),

                        a.jsoc.Notify("rrzhdanov@edu.hse.ru"))

    down2 = Fido.fetch(res_2).data
    MAP_full = fits.open(down2[0])
    dataMap, hdrMap = MAP_full[-1].data, MAP_full[-1].header
    dataMap[np.where(dataMap == 2)] = 1
    plt.imshow(dataMap)
    plt.figure()
    data_tarp, hdrtarp = MAP_mask[-1].data, MAP_mask[-1].header
    plt.imshow(data_tarp)
    print(hdrMap)
