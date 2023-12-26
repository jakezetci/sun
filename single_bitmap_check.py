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
import cpp_module as cpp
from matplotlib.colors import LogNorm
import matplotlib


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
values, points, areas = pipeline.bitmaps_to_points(TIME=False, onlyactive=True,
                                                   downloaded=True, magnetogram=[magnetogram_path],
                                                   bitmaps=[bitmap_path])
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

fig1, ax1 = plots.config(xlabel='X, пиксели', ylabel='Y, пиксели',
                         title='Магнетограма 2011/02/15', figsize=(4, 3),
                         dpi=140)
pos = ax1.imshow(dataMap, origin='lower', cmap='magma')
cbar = fig1.colorbar(pos)
cbar.ax.set_ylabel("лучевая компонента магнитного поля, Гаусс", rotation=270,
                   size='large', labelpad=15)
cbar.ax.yaxis.set_ticks_position("left")

dOBS = hdrMap["DSUN_OBS"]
pxsizeX, pxsizeY = hdrMap["CDELT1"], hdrMap["CDELT2"]


d_pixel_arcs = np.mean([pxsizeX, pxsizeY])

d_pixel = pipeline.arcsecs_to_radian(d_pixel_arcs) * dOBS
# active_indeces = np.argwhere(databitmap == 34)
# active_onmap = pipeline.bitmap_pixel_to_map(active_indeces, ref1, ref2)

"""
fig2, ax2 = plots.config(xlabel='X, пиксели', ylabel='Y, пиксели',
                         title='Активная область NOAA_AR = 11158 2011/02/15',
                         figsize=(16, 9))
pos2 = ax2.imshow(data_needed, origin='upper', cmap='inferno')
cbar2 = fig2.colorbar(pos2)
cbar2.ax.set_ylabel("лучевая компонента магнитного поля, Гаусс", rotation=270,
                    size='large', labelpad=15)

cbar2.ax.yaxis.set_ticks_position("left")
"""
centerX, centerY = hdrMap["CRPIX1"], hdrMap["CRPIX2"]
R = 696000.0*1e3
lats, lons, values_map = [], [], []

ref1, ref2 = int(hdrbitmap["CRPIX1"]), int(hdrbitmap["CRPIX2"])
mapsizeX, mapsizeY = int(hdrbitmap['CRSIZE1']), int(hdrbitmap['CRSIZE2'])

ref3, ref4 = ref1 + mapsizeX, ref2 + mapsizeY

data_needed = dataMap[ref2:ref4, ref1:ref3]

for i in range(ref2, ref4):
    for j in range(ref1, ref3):
        B = dataMap[i, j]
        x, y = -d_pixel * (j-centerX), -d_pixel * (i-centerY)
        lat, lon = coordinates.xy2ll(x, y)
        lats.append(lat)
        lons.append(lon)
        values_map.append(B)

lats = np.array(lats)
lons = np.array(lons)

area = np.full_like(lats, area_simple(ref1, ref2))
grid = lib.Grid(696000*1e3, lats, lons, area=area)
grid.set_allvalues(values_map, vector=False)


xlim = [-1e7, 3.1e8]
ylim = [-2.5*1e8, -9.85*1e7]

fig, ax = plots.plotmap(grid, n_linesx=8, n_linesy=8, xlimit=xlim, ylimit=ylim,
                        title='NOAA 11158 15 Feb 2011',
                        xlabel='X, m', ylabel='Y, m', grid=False,
                        figsize=(7, 5.5), dpi=120,
                        )

figZ, axZ = plots.config(ylabel='Z')

#################
active_indeces = np.argwhere(databitmap == 34)
# active_onmap = pipeline.bitmap_pixel_to_map(active_indeces, ref1, ref2)


computed = True


# 'maglines/REAL line v2.pkl', 'maglines/REAL line v3.pkl',
names = ['maglines/REAL line v14.pkl', 'maglines/REAL line v19.pkl',
         'maglines/REAL line v20.pkl', 'maglines/REAL line v21.pkl', 'maglines/REAL line v22.pkl',
         'maglines/REAL line v25.pkl', 'maglines/REAL line v26.pkl',
         'maglines/REAL line v45.pkl', 'maglines/REAL line v43.pkl', 'maglines/REAL line v47.pkl']
#         'maglines/REAL line v5.pkl', 'maglines/REAL line v6.pkl'
colors = ['green', 'pink', 'cyan', 'yellow', 'black',
          'red', 'grey', 'purple']
cmap = matplotlib.colormaps['tab20']
colors = cmap(np.linspace(0, 1, num=20))


coor_set = [[1.220*1e8, -1.472*1e8],
            [0.817*1e8, -1.646*1e8], [1.883*1e8, -1.761*1e8], [0.948e8, -1.849e8]]


x_1, y_1 = 0.797e8, -1.816e8
x_2, y_2 = 1.951e8, -1.337e8

X, Y = np.linspace(x_1, x_2, 8), np.linspace(y_1, y_2, 8)

coor_set = np.vstack([X, Y]).T[1:]

nums = np.arange(41, 49)

if computed == False:
    for i, (coor, num) in enumerate(zip(coor_set, nums)):
        xyz = np.array(coordinates.xyR2xyz(*coor, 696500*1e3))

        in_value = cpp.b_comp(xyz, values, points, areas)

        line_comp = lib.Magneticline(xyz, in_value, step=10*1e3)
        steps = 100000
        sign = math.copysign(1, in_value[2])
        line_comp = computing.comp_magneticline(line_comp, values, points, areas, returnobj=True,
                                                name=f'REAL line v{num}', maxsteps=steps,
                                                timestamp=1000, alert=False, stoppoint=696000*1e3,
                                                sign=sign)
        xx_comp, yy_comp, zz_comp = np.array(line_comp.points).T

        ax.plot(xx_comp, yy_comp, '--', color=colors[i],
                label=f'computed line {num}', lw=2)
        axZ.plot(xx_comp, zz_comp, '--', color=colors[i],
                 label=f'computed line {num}', lw=2)

else:
    for i, name in enumerate(names):

        with open(name, 'rb') as f_comp:
            line_comp = pickle.load(f_comp)
        try:
            xx_comp, yy_comp, zz_comp = np.array(line_comp.points).T
        except ValueError:
            xx_comp, yy_comp, zz_comp = np.array(line_comp.pointsxyz).T
        ax.plot(xx_comp, yy_comp, '-', color='white',
                lw=2, alpha=0.8)
        axZ.plot(xx_comp, zz_comp, '--', color=colors[i],
                 label=f'computed line {i+1}', lw=2)
ax.plot([], [], color='white', lw=2, label='расчётные магнитные линии')
ax.legend(loc='best')

#ax.legend(loc='best', fontsize='x-large')
