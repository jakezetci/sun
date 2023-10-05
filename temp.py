# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:37:14 2023

@author: cosbo
"""
import numpy as np
from coordinates import coordinates
from lib import cell, find_nearest
import pandas as pd
import pickle
class grid:
    def __init__(self, r, latitudes, longitudes, hs=False, area=False,
                 uniformgrid=True):
        if latitudes.size != longitudes.size:
            raise ValueError('longitudes size does not match latitudes')
        self.num = np.size(latitudes)
        self.values = np.zeros_like(latitudes)
        self.valuesvector = np.zeros((self.num, 3))
        self.r = r
        self.lat = latitudes
        self.lon = longitudes
        self.latlon = list(map(list, zip(latitudes, longitudes)))
        if uniformgrid is False:
            self.cells = []
            self.area = False
            for lat, lon, size in zip(latitudes, longitudes, hs):
                center = coordinates(r, lat, lon, latlon=True)
                self.cells.append(cell(center, size))
        elif hs is not False:
            for i in range(self.num):
                self.area[i] = self.r**2 * 
            self.area = self.r**2 * np.abs(np.cos(np.radians(90-self.lat))
                                           * hs)
        elif area is not False:
            self.area = area

        self.coors_set = [coordinates(self.r, *ll, latlon=True)
                          for ll in self.latlon]
        self.progress = 0

    def set_value(self, value, lat, lon, vector=False, easy=True, index=False):
        if index is not False:
            if vector:
                self.valuesvector[index] = value
            else:
                self.values[index] = value
        elif easy is True:
            self.set_value_easy(value, lat, lon, vector=vector)
        else:
            coor, i = find_nearest(self.latlon, (lat, lon), R=self.r)
            if vector:
                self.valuesvector[i] = value
            else:
                self.values[i] = value

    def set_value_easy(self, value, lat, lon, vector=False):
        i = np.where(((np.array(self.latlon)
                       == [lat, lon]).all(1)))
        if vector:
            self.valuesvector[i] = value
        else:
            self.values[i] = value

    def set_allvalues(self, values, vector=False):
        if vector is True:
            self.valuesvector = values
        else:
            self.values = values

    def find_value(self, lat, lon, easy=True, index=False):
        if index is not False:
            return self.values[index]
        elif easy is True:
            return self.find_value_easy(lat, lon)
        else:
            coor, ind = find_nearest(self.latlon, (lat, lon), R=self.r)
            return self.values[ind]

    def find_value_easy(self, lat, lon):
        ind = np.where(((np.array(self.latlon)
                         == [lat, lon]).all(1)))
        return self.values[ind]

    def save_pkl(self, name=False):
        if name is False:
            name = f'B_map {self.r:.2}'
        with open(f'{name}.pkl',
                  'wb') as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def dataframe(self):
        d = {'latlon': self.latlon, 'lat': self.lat, 'lon': self.lon,
             'r': np.full_like(self.lat, self.r),
             'xyz': [cs.vector for cs in self.coors_set], 'area': self.area,
             'values': self.values, 'valuesvector': self.valuesvector}
        return pd.DataFrame(data=d)

    def save_csv(self, name=False):
        if name is False:
            name = f'B_map {self.r:.2}'
        df = self.dataframe(self)
        df.to_csv(f'{name}.csv', index=True)

    def progress1(self):
        self.progress = self.progress + 1

    def progress_brute(self):
        self.progress = int(max([np.count_nonzero(self.values),
                             np.count_nonzero(self.valuesvector)/3]))

    def add_fieldinfo(self, m, dipolepos):
        self.m = m
        self.dipolepos = dipolepos

    def change_coors(self):
        self.coors_set = [coordinates(self.r, *ll, latlon=True)
                          for ll in self.latlon]