# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:17:33 2023

@author: cosbo
"""
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



def arcsecs_to_radian(arcsecs):
    return np.radians(arcsecs / 3600)

def file_name(time, series, general_path=r"C:/Users/cosbo/sunpy/data"):
    time_f = time.replace(':', '').replace('-', '').replace('T', '_')

    
    files = os.listdir(general_path)

    for file in files:
        if series in file and time_f in file:
            return general_path + '/' + file
    
    raise FileNotFoundError

def bitmap_pixel_to_map(index, reference1, reference2):
    correction = np.full_like(index, [reference2, reference1])
    return index+correction



def bitmaps_to_points(
    TIME, onlyactive=True, downloaded=False, magnetogram=None, bitmaps=None,
    returnhdr=False, plot=False, instrument=None,
):
    """
    takes data from bitmaps and returns them as a 2d array of coordinates and B-values

    Parameters
    ----------
    TIME : str or np.datetime64
        time parsed to Fido attributes, the data is downloaded of this timestamp
    onlyactive : bool, optional
        if True, only accounts for bitmap points == 34. if False, also takes bitmap==33
    downloaded : bool, optional
        if True, ignores TIME and takes paths in magnetogram, bitmaps instead
    magnetogram,  bitmaps : str, optional
        filepath for magnetograms and bitmaps.
    returnhdr : bool, optional
        if True, adds headers and center coors for

    Returns
    -------
    values, points, areas

    """
    

    def area_simple(xindex, yindex):

        # rotating to the first quadrant
        if xindex >= 0 and yindex < 0:
            xindex, yindex = -yindex, xindex
        elif xindex < 0 and yindex < 0:
            xindex, yindex = -xindex, -yindex
        elif xindex < 0 and yindex >= 0:
            xindex, yindex = yindex, -xindex
        elif xindex >= 0 and yindex >= 0:
            xindex, yindex = xindex, yindex

        tic = time.perf_counter()
        tan_dif = (np.arctan2(yindex + 1, xindex + 1) -
                   np.arctan2(yindex, xindex))
        sqr1 = np.sqrt(1 - d_r_ratio * (xindex**2 + yindex**2))
        sqr2 = np.sqrt(1 - d_r_ratio * ((xindex + 1) ** 2 + (yindex + 1) ** 2))
        toc = time.perf_counter()
        area = np.abs(tan_dif * (sqr1 - sqr2) * r_sun**2)
        # print(f'{toc-tic:.4} sec, {area}')
        return area

    if downloaded is False:
        magnetogram, bitmaps = download_map_and_harp(TIME, TIME, instrument)
    if isinstance(magnetogram, Iterable):
        magnetogram = np.asarray(magnetogram)
        MAP = fits.open(magnetogram[0])
    else:
        MAP = fits.open(magnetogram)

    dataMap, hdrMap = MAP[1].data, MAP[1].header
    #HMIMAP = sunpy.map.sources.HMIMap(dataMap, hdrMap)
    dOBS = hdrMap["DSUN_OBS"]
    pxsizeX, pxsizeY = hdrMap["CDELT1"], hdrMap["CDELT2"]

    # plt.imshow(dataMap)
    d_pixel_arcs = np.mean([pxsizeX, pxsizeY])

    d_pixel = arcsecs_to_radian(d_pixel_arcs) * dOBS

    r_sun = hdrMap["RSUN_REF"]
    d_r_ratio = (d_pixel / r_sun)
    d_r_ratio = d_r_ratio**2

    centerX, centerY = hdrMap["CRPIX1"], hdrMap["CRPIX2"]

    points = []
    values = []
    areas = []
    if returnhdr:
        headers = []

    if not isinstance(bitmaps, Iterable):
        bitmaps = [bitmaps]
    for bitmap_path in bitmaps:
        bitmap = fits.open(bitmap_path)

        databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
        if returnhdr:
            headers.append(hdrbitmap)
        ref1, ref2 = int(hdrbitmap["CRPIX1"]), int(hdrbitmap["CRPIX2"])
        # plt.plot(ref1, ref2, 'o', ms=12)
        active_indeces = np.argwhere(databitmap == 34)
        if onlyactive is False:
            quiet_indeces = np.argwhere(databitmap == 33)
            active_indeces = np.vstack([active_indeces, quiet_indeces])
        active_onmap = bitmap_pixel_to_map(active_indeces, ref1, ref2)
        try:
            mapsizeX, mapsizeY = int(hdrbitmap['CRSIZE1']), int(hdrbitmap['CRSIZE2'])
        except KeyError:
            mapsizeY, mapsizeX = np.shape(databitmap)
            hdrbitmap['CRSIZE1'] = mapsizeX
            hdrbitmap['CRSIZE2'] = mapsizeY
        for xindex, yindex in active_indeces:
            # plt.plot(yindex, xindex, 'o', ms=4, color='pink')
            tic = time.perf_counter()
            B = dataMap[xindex+ref2, yindex+ref1]
            x_corr, y_corr = -(yindex+ref1 - centerX), -(xindex+ref2 - centerY)

            x, y = d_pixel * x_corr, d_pixel * y_corr

            points.append(xyR2xyz(x, y, r_sun))
            toc = time.perf_counter()

            values.append(B)
            areas.append(area_simple(x_corr, y_corr))
    if plot is not False:
        
        ref3, ref4 = ref1 + mapsizeX, ref2 + mapsizeY
        B_s = dataMap[ref2:ref4, ref1:ref3]
        im = plot.imshow(B_s, origin='lower',
                    extent=np.array([-ref1+centerX, -ref3+centerX, -ref2+centerY, -ref4+centerY])*d_pixel,
                    aspect='equal', cmap='inferno')
        cbar = plt.colorbar(im, ax=plot, location='right')
        cbar.ax.set_ylabel("B_z, Гаусс",
                           size='medium')
        cbar.ax.tick_params(labelsize='small', direction='in')
        plot.set_ylabel('Y, м')
        plot.set_title('Магнитограмма активной области')
       

    if returnhdr:
        return np.array(values), np.array(points), np.array(areas), headers, (centerX, centerY)
    else:
        return np.array(values), np.array(points), np.array(areas)


def bitmaps_to_points_slow(
    TIME, onlyactive=True, downloaded=False, magnetogram=None, bitmaps=None,
    returnhdr=False
):
    # same function but it uses skycoord.pixel_to_world
    def area_simple(xindex, yindex):
        # delta is consistent in dumb method so i keep it
        xindex = xindex - centerX
        yindex = yindex - centerY
        tic = time.perf_counter()
        tan_dif = np.arctan2(yindex + 1, xindex + 1) - \
            np.arctan2(yindex, xindex)
        sqr1 = np.sqrt(1 - d_r_ratio * (xindex**2 + yindex**2))
        sqr2 = np.sqrt(1 - d_r_ratio * ((xindex + 1) ** 2 + (yindex + 1) ** 2))
        toc = time.perf_counter()
        area = np.abs(tan_dif * (sqr1 - sqr2) * r_sun**2)
        # print(f'{toc-tic:.4} sec, {area}')
        return area

    if downloaded is False:
        magnetogram, bitmaps = download_map_and_harp(TIME, TIME)
    magnetogram = np.asarray(magnetogram)
    MAP = fits.open(magnetogram[0])
    dataMap, hdrMap = MAP[1].data, MAP[1].header
    HMIMAP = sunpy.map.sources.HMIMap(dataMap, hdrMap)
    dOBS = hdrMap["DSUN_OBS"]
    pxsizeX, pxsizeY = hdrMap["CDELT1"], hdrMap["CDELT2"]

    d_pixel_arcs = np.mean([pxsizeX, pxsizeY])

    d_pixel = arcsecs_to_radian(d_pixel_arcs) * dOBS

    r_sun = hdrMap["RSUN_REF"]
    d_r_ratio = (d_pixel / r_sun)
    d_r_ratio = d_r_ratio**2

    centerX, centerY = hdrMap["CRPIX1"], hdrMap["CRPIX2"]
    points = []
    values = []
    areas = []
    if returnhdr:
        headers = []
    for i, bitmap_path in enumerate(bitmaps):
        bitmap = fits.open(bitmap_path)

        databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
        if returnhdr:
            headers.append(hdrbitmap)
        ref1, ref2 = hdrbitmap["CRPIX1"], hdrbitmap["CRPIX2"]
        active_indeces = np.argwhere(databitmap == 34)
        if onlyactive is False:
            quiet_indeces = np.argwhere(databitmap == 33)
            active_indeces = np.vstack([active_indeces, quiet_indeces])
        correction = np.full_like(active_indeces, [ref1, ref2])
        active_onmap = active_indeces + correction
        tic = time.perf_counter()
        for j, (xindex, yindex) in enumerate(active_onmap):
            B = dataMap[xindex, yindex]
            skycoord = HMIMAP.pixel_to_world(
                u.Quantity(xindex, unit="pix"), u.Quantity(yindex, unit="pix")
            )

            x = skycoord.heliocentric.x.value
            y = skycoord.heliocentric.y.value
            z = skycoord.heliocentric.z.value

            values.append(B)
            areas.append(area_simple(xindex, yindex))
            if j % 100 == 0:
                toc = time.perf_counter()
                print(
                    f"values {j-100}-{j} done in {toc - tic:0.2f} seconds")
                tic = time.perf_counter()
    if returnhdr:
        return (np.array(values), np.array(points), np.array(areas), headers, (centerX, centerY), HMIMAP)
    else:
        return np.array(values), np.array(points), np.array(areas)


def download_map_and_harp(timestart, timeend, instrument=None,
                          returnAR=False, path=None, **HARPkeywords):
    if instrument is None:
        timeHMI = np.datetime64('2010-05-01T00:00')
        if np.datetime64(timestart) > timeHMI:
            instrument = 'HMI'
        else:
            instrument = 'MDI'            
    
    if instrument == 'HMI':
        series_M = "hmi.M_720s"
        series_bitmap = "hmi.Mharp_720s"
    else:
        series_M = 'mdi.fd_M_96m_lev182'
        series_bitmap = 'mdi.lostarp_96m'
    res_M = Fido.search(
        a.Time(timestart, timeend),
        a.jsoc.Series(series_M),
        a.jsoc.Notify("rrzhdanov@edu.hse.ru"),
    )
    result = None
    while result is None:
        try:
            downloaded_magnetogram = Fido.fetch(res_M, path=path).data
            result = True
        except urllib.error.URLError:
            print('error in downloading')
            pass


    HARP_args = [a.Time(timestart, timeend),
                 a.jsoc.Series(series_bitmap),
                 a.jsoc.Segment("bitmap"),
                 a.jsoc.Notify("rrzhdanov@edu.hse.ru"),
                 ]

    for key, value in HARPkeywords.items():
        HARP_args.append(a.jsoc.Keyword(key) == value)

    res_bitmap = Fido.search(*HARP_args)
    result = None
    while result is None:
        try:
            downloaded_bitmaps = Fido.fetch(res_bitmap, path=path).data
            result = True
        except urllib.error.URLError:
            print('error in downloading')
            pass

    if returnAR:
        ar_nums = {}
        for name in downloaded_bitmaps:
            bitmap = fits.open(name)
            noaa = bitmap[-1].header['NOAA_AR']
            if instrument == "HMI":
                harp = bitmap[-1].header['HARPNUM']
            else:
                harp = bitmap[-1].header['TARPNUM']
            ar_nums[noaa] = harp
        return downloaded_magnetogram, downloaded_bitmaps, ar_nums
    return downloaded_magnetogram, downloaded_bitmaps


def compute_harp_MEnergy(
    timestart,
    timeend,
    onlyactive=True,
    mu=1.25e-6,
    downloaded=False,
    magnetogram=None,
    bitmaps=None,
):
    def get_area(xind: int, yind: int):
        def get_vector(ind1, ind2):
            start = time.process_time()
            MAPMAP = sunpy.map.sources.HMIMap(dataMap, hdrMap)
            hmitime = time.process_time()
            print(f"HMIMAP {hmitime-start:.2}")
            point = MAPMAP.pixel_to_world(
                u.Quantity(ind1, unit="pix"), u.Quantity(ind2, unit="pix")
            )
            pixtime = time.process_time()
            print(f"pixtime {pixtime-start:.2}")
            vector = np.array(
                [
                    point.heliocentric.x.value,
                    point.heliocentric.y.value,
                    point.heliocentric.z.value,
                ]
            )
            return vector

        vectorA = get_vector(xind - 1, yind + 1)
        vectorB = get_vector(xind + 1, yind + 1)
        vectorC = get_vector(xind + 1, yind - 1)
        vectorD = get_vector(xind - 1, yind - 1)

        sideAB = np.linalg.norm(vectorA - vectorB)
        sideBC = np.linalg.norm(vectorB - vectorC)
        sideCD = np.linalg.norm(vectorC - vectorD)
        sideDA = np.linalg.norm(vectorD - vectorA)

        diagAC = np.linalg.norm(vectorA - vectorC)
        diagBD = np.linalg.norm(vectorB - vectorD)

        area = 0.25 * np.sqrt(
            4 * diagAC**2 * diagBD**2
            - (sideAB**2 + sideCD**2 - sideBC**2 - sideDA**2) ** 2
        )
        return area

    def energy_onetime(dataMap, bitmaps):
        harps_count = len(bitmaps)
        energy_array = np.zeros(harps_count)
        for i, bitmap_path in enumerate(bitmaps):
            energy_harp = 0
            bitmap = fits.open(bitmap_path)

            databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
            size1, size2 = hdrbitmap["CRSIZE1"], hdrbitmap["CRSIZE2"]
            ref1, ref2 = hdrbitmap["CRPIX1"], hdrbitmap["CRPIX2"]
            active_indeces = np.argwhere(databitmap == 34)
            if onlyactive is False:
                quiet_indeces = np.argwhere(databitmap == 33)
                active_indeces = np.vstack([active_indeces, quiet_indeces])
            correction = np.full_like(active_indeces, [ref1, ref2])
            active_onmap = active_indeces + correction
            for xindex, yindex in active_onmap:
                B = dataMap[xindex, yindex]
                energy_piece = B**2 / (2 * mu)
                area = area_simple(xindex - centerX, yindex - centerY)
                energy_harp = +energy_piece * area
            energy_array[i] = energy_harp
            bitmap.close()
        return energy_array

    def area_simple(xindex, yindex):
        tic = time.perf_counter()
        tan_dif = np.arctan2(yindex + 1, xindex + 1) - \
            np.arctan2(yindex, xindex)
        sqr1 = np.sqrt(1 - d_r_ratio * (xindex**2 + yindex**2))
        sqr2 = np.sqrt(1 - d_r_ratio * ((xindex + 1) ** 2 + (yindex + 1) ** 2))
        toc = time.perf_counter()
        area = np.abs(tan_dif * (sqr1 - sqr2) * r_sun**2)
        # print(f'{toc-tic:.4} sec, {area}')
        return area

    if downloaded is False:
        magnetogram, bitmaps = download_map_and_harp(timestart, timeend)
    magnetogram = np.asarray(magnetogram)
    MAP = fits.open(magnetogram[0])
    dataMap, hdrMap = MAP[1].data, MAP[1].header
    pxsizeX, pxsizeY = hdrMap["CDELT1"], hdrMap["CDELT2"]
    d_pixel = np.mean([pxsizeX, pxsizeY])
    dOBS = hdrMap["DSUN_OBS"]

    r_sun = hdrMap["RSUN_REF"]
    d_r_ratio = (arcsecs_to_radian(d_pixel) * dOBS / r_sun) ** 2
    centerX, centerY = hdrMap["CRPIX1"], hdrMap["CRPIX2"]

    return energy_onetime(dataMap, bitmaps)


def fits_to_Grid(MAP):
    """
    thats a very long way to do it.......

    Args:
        MAP (TYPE): DESCRIPTION.

    Returns:
        TYPE: DESCRIPTION.

    """

    def get_area(xind: int, yind: int):
        def get_vector(ind1, ind2):
            point = MAP.pixel_to_world(
                u.Quantity(ind1, unit="pix"), u.Quantity(ind2, unit="pix")
            )
            vector = np.array(
                point.heliocentric.x.value,
                point.heliocentric.y.value,
                point.heliocentric.z.value,
            )
            return vector

        vectorA = get_vector(xind - 1, yind + 1)
        vectorB = get_vector(xind + 1, yind + 1)
        vectorC = get_vector(xind + 1, yind - 1)
        vectorD = get_vector(xind - 1, yind - 1)

        sideAB = np.linalg.norm(vectorA - vectorB)
        sideBC = np.linalg.norm(vectorB - vectorC)
        sideCD = np.linalg.norm(vectorC - vectorD)
        sideDA = np.linalg.norm(vectorD - vectorA)

        diagAC = np.linalg.norm(vectorA - vectorC)
        diagBD = np.linalg.norm(vectorB - vectorD)

        area = 0.25 * np.sqrt(
            4 * diagAC**2 * diagBD**2
            - (sideAB**2 + sideCD**2 - sideBC**2 - sideDA**2) ** 2
        )

        return area

    data, hdr = MAP[1].data, MAP[1].header
    MAP = sunpy.map.sources.HMIMap(data, hdr)
    nonNanindeces = np.argwhere(~np.isnan(data))
    gridN = nonNanindeces.shape[0]
    sunlats, sunlons, values, areas, R = np.zeros((5, gridN))
    for i, (xindex, yindex) in enumerate(nonNanindeces):
        skycoord = MAP.pixel_to_world(
            u.Quantity(xindex, unit="pix"), u.Quantity(yindex, unit="pix")
        )
        x = skycoord.heliocentric.x.value
        y = skycoord.heliocentric.y.value
        z = skycoord.heliocentric.z.value
        if np.isnan(x) or np.isnan(y) or np.isnan(z):
            pass
        else:
            R[i], (sunlats[i], sunlons[i]) = xyz2ll(x, y, z)
            values[i] = data[xindex, yindex]
            areas[i] = get_area(xindex, yindex)
    returnGrid = Grid(R, sunlats, sunlons, area=areas)

    returnGrid.set_allvalues(values, vector=False)
    return returnGrid


"""

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


maps = []
for file in downloaded_files:
    maps.append(fits.open(file))
    print(maps[-1].info())
data = maps[0][1].data

CR_kwords = list(filter(lambda x: x.startswith("CR"), kwords))
for word in CR_kwords:
    print(f'{word} = {hdr[word]}')
    res2 = Fido.search(a.Time('2023-08-08T00:00:00', '2023-08-08T00:20:00'),
                       a.jsoc.Series(series),
                       a.jsoc.Notify('rrzhdanov@edu.hse.ru'))

downloaded_files = Fido.fetch(res2)

print(downloaded_files)

Map = fits.open(
    r"C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9862.20230808_001200_TAI.bitmap.fits")
data = Map[-1].data
# plt.imshow(data, cmap='inferno')
hdr = Map[-1].header
kwords = list(hdr.keys())
print(hdr['CRPIX1'], hdr['CRPIX2'])
series_info = client.info(series)
CR_kwords = list(filter(lambda x: x.startswith("CR"), kwords))
for keyword in CR_kwords:
    keyword_info = series_info.keywords.loc[keyword]
    print(keyword)
    print(f"  type ....... {keyword_info.type} ")
    print(f"  recscope ... {keyword_info.recscope} ")
    print(f"  defval ..... {keyword_info.defval} ")
    print(f"  units ...... {keyword_info.units} ")

for word in CR_kwords:
    print(f'{word} = {hdr[word]}')
hmimap = sunpy.map.sources.HMIMap(data, hdr)


plt.imshow(data)
energys = np.zeros_like(dates, dtype=np.float64)
for i, date in enumerate(dates):
    energys[i] = np.sum(
        compute_harp_MEnergy(date, date, onlyactive=True, downloaded=False)
    )
np.savetxt("energys3.txt", energys)


energys = np.loadtxt("energys3.txt")
fig, ax = config(logscaley=True)
ax.plot(energys, "o")
"""
if __name__ == "__main__":

    series = "hmi.M_720s"
    # Create DRMS JSON client, use debug=True to see the query URLs
    client = drms.Client(email="rrzhdanov@edu.hse.ru")
    magnetogram = [
        "C:/Users/cosbo/sunpy/data/hmi.m_720s.20230808_001200_TAI.3.magnetogram.fits"
    ]
    bitmaps = [
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9915.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9914.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9913.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9911.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9909.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9907.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9905.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9904.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9902.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9900.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9899.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9893.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9892.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9890.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9882.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9879.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9875.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9872.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9864.20230808_001200_TAI.bitmap.fits",
        "C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9862.20230808_001200_TAI.bitmap.fits",
    ]
    # en_list = compute_harp_MEnergy('2023-08-08T00:14:00',
    #                               '2023-08-08T00:14:00', onlyactive=False,
    #                               downloaded=True, magnetogram=magnetogram,
    #                               bitmaps=bitmaps)
    dates = ["2023-08-08T00:12:00"]
    B, points = bitmaps_to_points(
        dates, downloaded=True, magnetogram=magnetogram, bitmaps=bitmaps
    )
    r = Coordinates(700000 * 1000, 60, 30, latlon=True)
    print(B_comp(r, B, points))
    # dates = pd.date_range(start="2023-08-08 00:12:00",
    #                      freq="60T", periods=8).values
