# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:37:22 2023

@author: Home-Pc
"""
import numpy as np
import math
import pandas as pd
from dataclasses import dataclass
import time
from coordinates import coordinates
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


def B_dipole(r, B0=1000, R=696340, returnBl=True, returnxyz=False):
    """    r класса coordinates

    """
    B_r = -2 * B0 * (R/r.r)**3 * math.cos(r.theta)
    B_theta = -B0 * (R/r.r)**3 * math.sin(r.theta)

    if returnBl is True:
        return ((B_r * math.sin(r.theta) +
                B_theta * math.cos(r.theta)) * math.cos(r.phi))
    elif returnxyz is True:
        return [B_r * math.sin(r.theta) * math.sin(r.phi) +
                B_theta * math.cos(r.theta) * math.sin(r.phi),
                B_r * math.cos(r.theta) + B_theta*math.sin(r.theta),
                B_r * math.sin(r.theta) * math.cos(r.phi) +
                B_theta * math.cos(r.theta) * math.cos(r.phi)]
    else:
        return [B_r, 0, B_theta]


def dipolebetter(r, m, rdipole=[0, 0, 0], returnBl=False,
                 returnxyz=False, mu=1.25e-6):
    """
    r класса coordinates, m вектор (x, y, z) в системе садыкова
    """

    c = mu/(4 * math.pi * np.linalg.norm(r.vector)**3)
    m = np.asarray(m)
    rm = np.dot(r.vector, m)

    r_part = c * 3 * rm * r.vector
    m_part = - c * m
    L = np.asarray([0, 0, 1])
    if returnBl is True:
        return np.dot(r_part, L) + np.dot(m_part, L)
    elif returnxyz is True:
        return r_part + m_part




def distance_sphere(a, b, R):

    lat1, lon1 = a
    lat2, lon2 = b
    phi_1 = math.radians(lat1)
    phi_2 = math.radians(lat2)

    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    d = math.sin(delta_phi / 2.0) ** 2 + (math.cos(phi_1) * math.cos(phi_2)
                                          * math.sin(delta_lambda / 2.0) ** 2)
    try:
        c = 2 * math.atan2(math.sqrt(d), math.sqrt(1 - d))
    except ValueError:
        print(f'd={d}')
        print(f'a, b = {a, b}')
        raise ValueError
    km = R * c
    return km


def find_nearest(array, point, R=696340):
    array = list(map(list, array))
    dist_arr = np.asarray([np.abs(distance_sphere(
        a, point, R)) for a in array])
    idx = dist_arr.argmin()
    return array[idx], idx

# from here on now axis Z refers to the line of sight
# line of sight is meant to be lined up with the centre of the Sun, (0,0)
# so far we assume that sun isn't moving so line of sight is at (0,0)


def I1(r1, r2):
    x, y, z = r1 - r2  # разница векторов
    return z/((x**2 + y**2)*np.linalg.norm(r1-r2))


def I2(r1, r2):
    x, y, z = r1 - r2  # разница векторов
    num = z * (3 * x**2 + 3 * y**2 + 2 * z**2)
    den = 3 * (x**2 + y**2)**2 * (np.linalg.norm(r1-r2))**3
    a =num/den
    return num/den


def I3(r1, r2):
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    num = (-2*((x**2+y**2)**2) * z2)+((x**2+y**2)*(z**3 + 3*(z2**2)*z)) - 2*(z2**2)*(-z)**3
    den = 3 * (x**2 + y**2)**2 * (np.linalg.norm(r1-r2))**3
    return num/den


def lI1(r1, r2):
    x, y, z = r1 - r2
    return 1/(x**2 + y**2)


def lI2(r1, r2):
    x, y, z = r1 - r2
    return 2/(3 * (x**2 + y**2)**2)


def lI3(r1, r2):
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    return (x**2 + y**2 + 2*(z2**2))/(3 * (x**2 + y**2)**2)


def IlI1(r1, r2):
    num = -1
    R = np.asarray(r1 - r2)
    r = np.linalg.norm(R)
    L = np.array([0, 0, 1])
    den = (r + np.dot(L, R))*r
    return num / den


def IlI2(r1, r2):
    R = np.asarray(r1 - r2)
    r = np.linalg.norm(R)
    L = np.array([0, 0, 1])
    d1 = -2/((3*r*(r-np.dot(L,R)) * (r+np.dot(L,R))**2))
    d2 = np.dot(L, R) / (3*(r**2 - np.dot(L, R)**2)*r**3)
    return d1+d2


def IlI3(r1, r2):
    R = np.asarray(r1 - r2)
    r = np.linalg.norm(R)
    L = np.array([0, 0, 1])
    d1 = -(np.dot(L,R)**2 + np.dot(L,R)*r + r**2)/(3*(r+np.dot(L,R)) * r**3)
    d2 = -2*np.dot(L,r2)/(3*r**3)
    d3n = -(np.dot(L,r2)**2 * (2*r**3 - 3*np.dot(L,R)*r**2 + np.dot(L,R)**3))
    d3d = 3*r**3 * (r**2 - np.dot(L,R)**2)**2
    d3 = d3n/d3d
    return d1+d2+d3


def GreenBl(r1, r2, a=696340, method=False):
    r1, r2 = r1.vector, r2.vector
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    if method is True:
        ili1 = IlI1(r1, r2)
        ili2 = IlI2(r1, r2)
        ili3 = IlI3(r1, r2)
        G1 = 2 * x1 * ili1 - 3*x*((x1**2 + y1**2 - a**2)*ili2+ili3)
        G2 = 2 * y1 * ili1 - 3*y*((x1**2 + y1**2 - a**2)*ili2+ili3)
        G3 = (np.linalg.norm(r1)**2 - a**2) / (np.linalg.norm(r1-r2)**3)
    else:
        i1, i2, i3 = I1(r1, r2), I2(r1, r2), I3(r1, r2)
        li1, li2, li3 = lI1(r1, r2), lI2(r1, r2), lI3(r1, r2)
        G1 = 2*x1*(i1-li1) - 3*x*((x1**2 + y1**2 - a**2)*(i2-li2)+(i3-li3))
        G2 = 2*y1*(i1-li1) - 3*y*((x1**2 + y1**2 - a**2)*(i2-li2)+(i3-li3))
        G3 = (np.linalg.norm(r1)**2 - a**2) / (np.linalg.norm(r1-r2)**3)
    return np.array([G1, G2, G3])/(4*np.pi*a)


class cell:
    def __init__(self, center, hs, sizeType='deg'):
        self.center = center
        self.leftborder, self.rightborder = center.lon - hs, center.lon + hs
        # we pretend that cells are small enough to be considered plane
        self.area = center.r**2 * 4 * hs**2
        self.value = 0

    def set_value(self, value):
        self.value = value

class B_class:
    def __init__(self, B, progress):
        self.B = B
        self.progress = progress

    def prog(self):
        self.progress = self.progress + 1

    def B(self, B_new):
        self.B = B_new

@dataclass
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
            self.area = self.r**2 * np.abs((np.cos(np.radians(90-self.lat)+hs/2) -
                                           np.cos(np.radians(90-self.lat)-hs/2))
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

def B_comp(r, B_map, method='allmap', grid=False, debug=False, change=True):

    def B_comp_same_maps(r, B_map, heavy=True, change=True):
        """
        быстрое вычисление B в случае, когда
        сетка с числам совпадает с сеткой интегрирования
        Args:
            r (coordinates): DESCRIPTION.
            B_map (grid): DESCRIPTION.

        Returns:
            B_comp(np.array, 3d)

        """
        B = np.asarray([0.0, 0.0, 0.0], dtype=np.float32)
        R = B_map.r
        n = B_map.num
        debugging_obj = []
        if heavy:
            for i in range(n):

                lat, lon = B_map.latlon[i]
                r_2 = coordinates(R, lat, lon, latlon=True)
                B_l = B_map.find_value(lat, lon, index=i)
                add = B_l * np.asarray(GreenBl(r, r_2, method=change)) * B_map.area[i]
                B = B + add
                if debug:
                    debugging_obj.append(add)
                    if  6 < lat < 12 and 10 < lon < 14:
                        gr = GreenBl(r, r_2)
                        r1, r_2 = r.vector, r_2.vector
                        i1, i2, i3 = I1(r1, r_2), I2(r1, r_2), I3(r1, r_2)
                        li1, li2, li3 = lI1(r1, r_2), lI2(r1, r_2), lI3(r1, r_2)
                        a1, a2, a3 = i1 - li1, i2 - li2, i3 -li3
                        g = gr
                        aili1 = IlI1(r1, r_2)
                        aili2 = IlI2(r1, r_2)
                        aili3 = IlI3(r1, r_2)
                        gg = 0

        else:
            for coor, S in zip(B_map.latlon, B_map.area):
                lat, lon = coor
                r_2 = coordinates(R, lat, lon, latlon=True)
                B_l = B_map.find_value_easy(lat, lon)
                add = B_l * np.asarray(GreenBl(r, r_2)) * S
                B = B + add
        if debug:
            return B, debugging_obj
        else:
            return B

    def B_comp_diff_maps(r, grid, B_map):
        """
        r класса coordinates
        """
        B = np.asarray([0.0, 0.0, 0.0], dtype=np.float32)
        R = grid.r
        for coor, S in zip(grid.latlon, grid.area):
            lat, lon = coor
            r_2 = coordinates(R, lat, lon, latlon=True)
            B_l = B_map.find_value(lat, lon)

            add = B_l * np.asarray(GreenBl(r, r_2)) * S
            B = B + add
        return B

    if method == 'allmap':
        return B_comp_same_maps(r, B_map, change=change)
    elif method == 'diffmaps':
        return B_comp_diff_maps(r, grid, B_map)


class magneticline:
    def __init__(self, initial_point, initial_value, step):
        self.initial_point = initial_point
        self.initial_value = initial_value
        self.points = [initial_point]
        self.values = [initial_value]
        self.pointsxyz = [initial_point.vector]
        self.step = step
        self.progress = 0


    def add_value(self, func, *args):
        vec = np.asarray(self.values[-1])
        new_point = vec * self.step / np.linalg.norm(vec) + np.asarray(
            self.points[-1].vector)

        self.values.append(func(coordinates(*new_point), returnxyz=True, *args))
        self.points.append(coordinates(*new_point))
        self.pointsxyz.append(new_point)
        self.progress = self.progress + 1

    def add_value_comp(self, B_map):
        vec = np.asarray(self.values[-1])
        new_point = vec * self.step / np.linalg.norm(vec) + np.asarray(
            self.points[-1].vector)
        
        val = B_comp(coordinates(*new_point), B_map)
        self.values.append(val)
        self.points.append(coordinates(*new_point))
        self.pointsxyz.append(new_point)
        self.progress = self.progress + 1


    def line_by_length(self, func, length, *args):
        for i in range(length//self.step):
            self.add_value(self, func, *args)

    def line_by_stoppoint(self, func, stoppoint, maxsteps=1e3, *args):
        for i in range(maxsteps):
            self.add_value(self, func, *args)
            if self.points[-1].r < stoppoint:
                break

    def save_pkl(self, name=False):
        if name is False:
            name = f'{self.initial_point.r:.3}-{self.initial_value:.3}'
        with open(f'{name}.pkl',
                  'wb') as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def dataframe(self):
        d = {'points': self.points, 'values': self.values}
        return pd.DataFrame(data=d)

    def save_csv(self, name=False):
        df = self.dataframe(self)
        df.to_csv(f'{name}.csv')

    def brute_progress(self):
        self.progress = len(self.points)


def load_grid(df):
    r = df.r[0]
    Grid = grid(r, df.lat, df.lon, area=df.area)
    Grid.set_allvalues(df.values)
    Grid.set_allvalues(df.valuesvector, vector=True)
    return Grid


if __name__ == "__main__":
    pass
