import pfsspy
import pfsspy.utils

import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
import pipeline
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import frame_transform_graph
import sunpy.coordinates as spc
import astropy.coordinates as coord
from astropy.io import fits
import computing


date = '2011-02-11'

# use time for gong maps
time_hms = 'hh:mm:ss'


class ASineHGC(spc.HeliographicCarrington):
    """
    This frame is a reverse transformation from the HGC frame with CUNIT2 sin(deg) to the real HGC frame.
    It is defined backwards because the WCS object for the HGC sin(deg) frame reads in as a proper HGC frame
    due to the values of CTYPE. This frame therefore converts from what WCS thinks is HGC to ASineHGC but is
    in reality converting from ASineHGC to HGC.
    """


@frame_transform_graph.transform(coord.FunctionTransform, spc.HeliographicCarrington, ASineHGC)
def hgc_to_sinehgc(hgc_coord, sinehgc_frame):
    lat = hgc_coord.lat
    lon = hgc_coord.lon

    lat_out = u.Quantity(np.arcsin(lat.value), u.rad)

    return sinehgc_frame.realize_frame(coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


@frame_transform_graph.transform(coord.FunctionTransform, ASineHGC, spc.HeliographicCarrington)
def sinehgc_to_hgc(sinehgc_coord, hgc_frame):
    lat = sinehgc_coord.lat
    lon = sinehgc_coord.lon

    lat_out = u.Quantity(np.sin(lat.value), u.deg)

    return hgc_frame.realize_frame(coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


def hmi_output(date, nrho, rss, resample_size=[360, 180],
               plot_enable=False):
    # get the files in 24h period
    date = np.datetime64(date)
    date2 = date + np.timedelta64(1, 'D')
    time = a.Time(date, date2)
    series = a.jsoc.Series('hmi.mldailysynframe_720s')
    result = Fido.search(time, series,
                         a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    files = Fido.fetch(result)

    # repair files - from sine lat to degrees

    hdus = fits.open(files[0])
    data = hdus[0].data
    header = dict(hdus[0].header)
    header['CUNIT2'] = 'deg'
    header['CDELT2'] = 180 / np.pi * header['CDELT2']
    del header['HGLN_OBS']

    # use crot data to repair NaNs
    crot = a.jsoc.PrimeKey('CAR_ROT',
                           header['CAR_ROT'])
    result_repair = Fido.search(crot, a.jsoc.Series('hmi.synoptic_mr_polfil_720s'),
                                a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    files_repair = Fido.fetch(result_repair)

    map_repair = sunpy.map.Map(files_repair[0])
    crot_data = map_repair.data
    if plot_enable:
        plt.figure()
        ax2 = plt.subplot(projection=map_repair)
        im2 = map_repair.plot(axes=ax2)

    data[np.isnan(data)] = crot_data[np.isnan(data)]

    # load map to sunpy
    hmi_map = sunpy.map.Map((data, header))

    if plot_enable:
        plt.figure()
        hmi_map.plot_settings['cmap'] = 'hmimag'
        hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
        ax = plt.subplot(projection=hmi_map)
        im = hmi_map.plot(axes=ax)
    hmi_map = hmi_map.resample(resample_size * u.pix)

    # perform pfss
    pfss_in = pfsspy.Input(hmi_map, nrho, rss)
    pfss_out = pfsspy.pfss(pfss_in)
    return pfss_out


def gong_output(gong_map, nrho, rss):
    nrho = 35
    rss = 2.5
    pfss_in = pfsspy.Input(gong_map, nrho, rss)
    return pfsspy.pfss(pfss_in)


def get_energy_pffs(time, pfss_output, NOAA_AR, density):
    values, points, areas, hdrs, (cX, cY) = computing.bitmaps_to_points(TIME=time,
                                                                        onlyactive=True,
                                                                        returnhdr=True,
                                                                        plot=False)
    grid = computing.create_3Dgrid(hdrs[0], density, cX, cY, mode='fineZ')
    # turn points to heliocentric cartesian astropy coors
    return 0


if __name__ == "__main__":
    pfss_out = hmi_output(date, 35, 2.5)
