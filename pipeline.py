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
import sunpy.map.sources
from coordinates import xyz2ll
from lib import Grid


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
            point = MAP.pixel_to_world(u.Quantity(ind1, unit='pix'),
                                       u.Quantity(ind2, unit='pix'))
            vector = np.array(point.heliocentric.x.value,
                              point.heliocentric.y.value, point.heliocentric.z.value)
            return vector

        vectorA = get_vector(xind-1, yind+1)
        vectorB = get_vector(xind+1, yind+1)
        vectorC = get_vector(xind+1, yind-1)
        vectorD = get_vector(xind-1, yind-1)

        sideAB = np.linalg.norm(vectorA-vectorB)
        sideBC = np.linalg.norm(vectorB-vectorC)
        sideCD = np.linalg.norm(vectorC-vectorD)
        sideDA = np.linalg.norm(vectorD-vectorA)

        diagAC = np.linalg.norm(vectorA-vectorC)
        diagBD = np.linalg.norm(vectorB-vectorD)

        area = 0.25 * np.sqrt(4*diagAC**2 * diagBD**2 -
                              (sideAB**2 + sideCD**2 - sideBC**2 - sideDA ** 2) ** 2)

        return area

    data, hdr = MAP[1].data, MAP[1].header
    MAP = sunpy.map.sources.HMIMap(data, hdr)
    nonNanindeces = np.argwhere(~np.isnan(data))
    gridN = nonNanindeces.shape[0]
    sunlats, sunlons, values, areas, R = np.zeros((5, gridN))
    for i, (xindex, yindex) in enumerate(nonNanindeces):
        skycoord = MAP.pixel_to_world(u.Quantity(
            xindex, unit='pix'), u.Quantity(yindex, unit='pix'))
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


def download_map_and_harp(timestart, timeend):
    series_M = 'hmi.M_720s'
    series_bitmap = 'hmi.Mharp_720s'
    res_bitmap = Fido.search(a.Time(timestart, timeend),
                             a.jsoc.Series(series_M),
                             a.jsoc.Notify('rrzhdanov@edu.hse.ru'))

    downloaded_magnetogram = Fido.fetch(res_bitmap).data

    res_M = Fido.search(a.Time(timestart, timeend),
                        a.jsoc.Series(series_bitmap),
                        a.jsoc.Segment('bitmap'),
                        a.jsoc.Notify('rrzhdanov@edu.hse.ru'))
    downloaded_bitmaps = Fido.fetch(res_M).data
    return downloaded_magnetogram, downloaded_bitmaps


def compute_harp_MEnergy(timestart, timeend, onlyactive=True, mu=1.25e-6):
    def energy_onetime(dataMap, bitmaps):
        harps_count = len(bitmaps)
        energy_array = np.zeros(harps_count)
        for i, bitmap_path in enumerate(bitmaps):
            energy_harp = 0
            bitmap = fits.open(bitmap_path)

            databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
            size1, size2 = hdrbitmap['CRSIZE1'], hdrbitmap['CRSIZE2']
            ref1, ref2 = hdrbitmap['CRPIX1'], hdrbitmap['CRPIX2']
            active_indeces = np.argwhere(databitmap == 34)
            if onlyactive is False:
                quiet_indeces = np.argwhere(databitmap == 33)
                active_indeces = np.vstack([active_indeces, quiet_indeces])
            correction = np.full_like(active_indeces, [ref1, ref2])
            active_onmap = active_indeces + correction
            for (xindex, yindex) in active_onmap:
                B = dataMap[xindex, yindex]
                energy_piece = B**2 / (2*mu)
                energy_harp = + energy_piece
            energy_array[i] = energy_harp
            bitmap.close()
        return energy_array

    magnetogram, bitmaps = download_map_and_harp(timestart, timeend)
    MAP = fits.open(magnetogram[0])
    dataMap, hdrMap = MAP[1].data, MAP[1].header

    return energy_onetime(dataMap, bitmaps)


series = 'hmi.Mharp_720s'
# Create DRMS JSON client, use debug=True to see the query URLs
client = drms.Client(email='rrzhdanov@edu.hse.ru')


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
#plt.imshow(data, cmap='inferno')
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
"""

en_list = compute_harp_MEnergy('2023-08-08T00:14:00',
                               '2023-08-08T00:14:00', onlyactive=True)
# Map = fits.open( 'C:/Users/cosbo/sunpy/data/hmi.mharp_720s.9862.20230808_002400_TAI.bitmap.fits')
