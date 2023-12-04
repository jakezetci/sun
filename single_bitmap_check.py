# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 18:39:38 2023

@author: cosbo
"""


from astropy.io import fits
import drms
import numpy as np
import sunpy
from sunpy.net import Fido, attrs as a
import astropy.units as u
import sunpy.map.sources
import matplotlib.pyplot as plt
import pipeline
import coordinates
import lib
import plots
import computing
import pickle
import math


def area_simple(xindex, yindex, d_r_ratio=2.692592614123163e-07, r_sun=696000*1e3):
    xindex = xindex - centerX
    yindex = yindex - centerY
    tan_dif = (np.arctan2(yindex + 1, xindex + 1) -
               np.arctan2(yindex, xindex))
    sqr1 = np.sqrt(1 - d_r_ratio * (xindex**2 + yindex**2))
    sqr2 = np.sqrt(1 - d_r_ratio * ((xindex + 1) ** 2 + (yindex + 1) ** 2))
    area = np.abs(tan_dif * (sqr1 - sqr2) * r_sun**2)
    return area


"""
NOAA_AR= 11158
"""

time = "2011-02-15 00:01:00"

magnetogram_path = r'C:\Users\cosbo\sunpy\data\hmi.m_720s.20110215_000000_TAI.1.magnetogram.fits'
bitmap_path = r'C:\Users\cosbo\sunpy\data\hmi.mharp_720s.377.20110215_000000_TAI.bitmap.fits'
bitmap = fits.open(bitmap_path)
databitmap, hdrbitmap = bitmap[-1].data, bitmap[-1].header
"""
print(hdrbitmap['HARPNUM'])
print(hdrbitmap['NOAA_AR'])
print(databitmap.shape)
plt.imshow(databitmap, origin='lower')
active = np.shape(np.where(databitmap == 34))[1]
nonactive = np.shape(np.where(databitmap == 33))[1]
print(f'active = {active}')
print(f'slightly active  = {nonactive}')
print(
    f'estimated time for coors = {0.016*active} seconds or {0.016*(nonactive+active)/60} minutes')
"""
MAP = fits.open(magnetogram_path)
dataMap, hdrMap = MAP[-1].data, MAP[-1].header
plt.figure()
plt.imshow(dataMap, origin='lower')
ref1, ref2 = int(hdrbitmap["CRPIX1"]), int(hdrbitmap["CRPIX2"])
mapsizeX, mapsizeY = int(hdrbitmap['CRSIZE1']), int(hdrbitmap['CRSIZE2'])
dOBS = hdrMap["DSUN_OBS"]
pxsizeX, pxsizeY = hdrMap["CDELT1"], hdrMap["CDELT2"]


d_pixel_arcs = np.mean([pxsizeX, pxsizeY])

d_pixel = pipeline.arcsecs_to_radian(d_pixel_arcs) * dOBS
#active_indeces = np.argwhere(databitmap == 34)
#active_onmap = pipeline.bitmap_pixel_to_map(active_indeces, ref1, ref2)
ref3, ref4 = ref1 + mapsizeX, ref2 + mapsizeY
data_needed = dataMap[ref2:ref4, ref1:ref3]
plt.figure()
plt.imshow(data_needed, origin='lower')
centerX, centerY = hdrMap["CRPIX1"], hdrMap["CRPIX2"]
R = 696000.0*1e3
lats, lons, values = [], [], []

for i in range(ref2, ref4):
    for j in range(ref1, ref3):
        B = dataMap[i, j]
        x, y = -d_pixel * (j-centerX), -d_pixel * (i-centerY)
        lat, lon = coordinates.xy2ll(x, y)
        lats.append(lat)
        lons.append(lon)
        values.append(B)

lats = np.array(lats)
lons = np.array(lons)

area = np.full_like(lats, area_simple(ref1, ref2))
grid = lib.Grid(696000*1e3, lats, lons, area=area)
grid.set_allvalues(values, vector=False)


xlim = [0.0, 3*1e8]
ylim = [-2.5*1e8, -9.85*1e7]
"""fig, ax = plots.plotmap(grid, n_lines=8, xlimit=xlim, ylimit=ylim,
                        title='NOAA 11158 15 Feb 2011',
                        xlabel='X, m', ylabel='Y, m')
"""
fig, ax = plots.config(ylabel='Z')

#################
active_indeces = np.argwhere(databitmap == 34)
#active_onmap = pipeline.bitmap_pixel_to_map(active_indeces, ref1, ref2)
lats_new, lons_new, values_new = [], [], []
for i, j in active_indeces:
    B = dataMap[i+ref2, j+ref1]
    x, y = -d_pixel * (i-centerX), -d_pixel * (j-centerY)
    x_corr, y_corr = -(i+ref1 - centerX), -(j+ref2 - centerY)

    x, y = d_pixel * x_corr, d_pixel * y_corr

    lat, lon = coordinates.xy2ll(x, y)
    lats_new.append(lat)
    lons_new.append(lon)
    values_new.append(B)

lats_new = np.array(lats_new)
lons_new = np.array(lons_new)
Lmap = lib.grid(696000*1e3, lats_new, lons_new, area=area)
Lmap.set_allvalues(values_new, vector=False)

computed = True


# 'maglines/REAL line v2.pkl', 'maglines/REAL line v3.pkl',
names = ['maglines/REAL line v10.pkl', 'maglines/REAL line v7.pkl',
         'maglines/REAL line v11.pkl']
#         'maglines/REAL line v5.pkl', 'maglines/REAL line v6.pkl']
colors = ['green', 'pink', 'cyan', 'yellow', 'black',
          'red', 'grey']

if computed == False:

    xyz = coordinates.xyR2xyz(1.54*1e8, -1.54*1e8, 696000*1e3)

    start_point = coordinates.Coordinates(*xyz)
    in_value = lib.B_comp_map(start_point, Lmap)

    line_comp = lib.Magneticline(start_point, in_value, step=300*1e3)
    steps = 100
    sign = math.copysign(1, in_value[2])
    line_comp = computing.comp_magneticline(line_comp, Lmap, returnobj=True,
                                            name=f'REAL line v12', maxsteps=steps,
                                            timestamp=10, alert=False, stoppoint=696000*1e3,
                                            sign=sign)
    xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T

    ax.plot(xx_comp, zz_comp, '--', lw=2)
else:
    for i, name in enumerate(names):

        with open(name, 'rb') as f_comp:
            line_comp = pickle.load(f_comp)
        xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
        ax.plot(xx_comp, zz_comp, '-', color=colors[i],
                label=f'computed line {i+1}', lw=2)


ax.legend(loc='best', fontsize='x-large')
