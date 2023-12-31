# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:37:22 2023

@author: Home-Pc
"""
import numpy as np
import math
import pandas as pd
from dataclasses import dataclass
import cpp_module as cpp
import time
try:

    from coordinates import Coordinates
    from field import dipolebetter
except ModuleNotFoundError:

    from sun.field import dipolebetter
    from sun.coordinates import Coordinates

import pickle


def distance_sphere(a, b, R):
    """

    Parameters
    ----------
    a
        [lat, lon]
    b
        [lat, lon]
    R
        radius of sphere

    Returns
    -------
        distance

    Raises
    ------
    ValueError
    """
    lat1, lon1 = a
    lat2, lon2 = b
    phi_1 = math.radians(lat1)
    phi_2 = math.radians(lat2)

    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    d = math.sin(delta_phi / 2.0) ** 2 + (
        math.cos(phi_1) * math.cos(phi_2) * math.sin(delta_lambda / 2.0) ** 2
    )
    try:
        c = 2 * math.atan2(math.sqrt(d), math.sqrt(1 - d))
    except ValueError:
        print(f"d={d}")
        print(f"a, b = {a, b}")
        raise ValueError
    m = R * c
    return m


def find_nearest(array, point, R=696340 * 1e3):
    array = list(map(list, array))
    dist_arr = np.asarray(
        [np.abs(distance_sphere(a, point, R)) for a in array])
    idx = dist_arr.argmin()
    return array[idx], idx


def GreenBl(r1, r2, a=696340 * 1000, vector_method=False, debug=False):
    """
    Computes Green`s function based on Sadykov-Zimovets work
    This is a readable version with formuale being consistent with the article
    Parameters
    ----------
    r1 : Coordinates or np.array
        point 1 (where the field is computed)
    r2 : Coordinates or np.array
        point 2 (point on a grid)
    a : float, optional
        sun radius. The default is 696340 * 1000.
    vector_method : Bool, optional
        Either using vector method of not,
            details in the work. The default is False.

    Returns
    -------
    GreenBl : float
    """

    def I1(r1, r2):
        x, y, z = r1 - r2
        return z / ((x**2 + y**2) * np.linalg.norm(r1 - r2))

    def I2(r1, r2):
        x, y, z = r1 - r2
        num = z * (3 * x**2 + 3 * y**2 + 2 * z**2)
        den = 3 * (x**2 + y**2) ** 2 * (np.linalg.norm(r1 - r2)) ** 3
        return num / den

    def I3(r1, r2):
        x, y, z = r1 - r2
        x1, y1, z1 = r1
        x2, y2, z2 = r2
        num = (
            (-2 * ((x**2 + y**2) ** 2) * z2)
            + ((x**2 + y**2) * (z**3 + 3 * (z2**2) * z))
            - 2 * (z2**2) * (-z) ** 3
        )
        den = 3 * (x**2 + y**2) ** 2 * (np.linalg.norm(r1 - r2)) ** 3
        return num / den

    def lI1(r1, r2):
        x, y, z = r1 - r2
        return 1 / (x**2 + y**2)

    def lI2(r1, r2):
        x, y, z = r1 - r2
        return 2 / (3 * (x**2 + y**2) ** 2)

    def lI3(r1, r2):
        x, y, z = r1 - r2
        x1, y1, z1 = r1
        x2, y2, z2 = r2
        return (x**2 + y**2 + 2 * (z2**2)) / (3 * (x**2 + y**2) ** 2)

    def IlI1(r1, r2):
        num = -1
        R = np.asarray(r1 - r2)
        r = np.linalg.norm(R)
        L = np.array([0, 0, 1])
        den = (r + np.dot(L, R)) * r
        return num / den

    def IlI2(r1, r2):
        R = np.asarray(r1 - r2)
        r = np.linalg.norm(R)
        L = np.array([0, 0, 1])
        d1 = -2 / ((3 * r * (r - np.dot(L, R)) * (r + np.dot(L, R)) ** 2))
        d2 = np.dot(L, R) / (3 * (r**2 - np.dot(L, R) ** 2) * r**3)
        return d1 + d2

    def IlI3(r1, r2):
        R = np.asarray(r1 - r2)
        r = np.linalg.norm(R)
        L = np.array([0, 0, 1])
        d1 = -(np.dot(L, R) ** 2 + np.dot(L, R) * r + r**2) / (
            3 * (r + np.dot(L, R)) * r**3
        )
        d2 = -2 * np.dot(L, r2) / (3 * r**3)
        d3n = -(
            np.dot(L, r2) ** 2
            * (2 * r**3 - 3 * np.dot(L, R) * r**2 + np.dot(L, R) ** 3)
        )
        d3d = 3 * r**3 * (r**2 - np.dot(L, R) ** 2) ** 2
        d3 = d3n / d3d
        return d1 + d2 + d3

    if type(r1) == Coordinates:
        r1 = r1.vector
    if type(r2) == Coordinates:
        r2 = r2.vector
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    if vector_method is True:
        ili1 = IlI1(r1, r2)
        ili2 = IlI2(r1, r2)
        ili3 = IlI3(r1, r2)
        G1 = 2 * x1 * ili1 - 3 * x * ((x1**2 + y1**2 - a**2) * ili2 + ili3)
        G2 = 2 * y1 * ili1 - 3 * y * ((x1**2 + y1**2 - a**2) * ili2 + ili3)
        G3 = (np.linalg.norm(r1) ** 2 - a**2) / (np.linalg.norm(r1 - r2) ** 3)
    else:
        i1, i2, i3 = I1(r1, r2), I2(r1, r2), I3(r1, r2)
        li1, li2, li3 = lI1(r1, r2), lI2(r1, r2), lI3(r1, r2)
        ili1 = i1 - li1
        ili2 = i2 - li2
        ili3 = i3 - li3
        G1 = 2 * x1 * (i1 - li1) - 3 * x * (
            (x1**2 + y1**2 - a**2) * (i2 - li2) + (i3 - li3)
        )
        G2 = 2 * y1 * (i1 - li1) - 3 * y * (
            (x1**2 + y1**2 - a**2) * (i2 - li2) + (i3 - li3)
        )
        G3 = (np.linalg.norm(r1) ** 2 - a**2) / (np.linalg.norm(r1 - r2) ** 3)

    return np.array([G1, G2, G3]) / (4 * np.pi * a)


def Green_optimized(r1, r2, a=696300*1e3):

    def I1():
        return z / (x_sqr_y_sqr * r_diff_norm)

    def I2():
        num = z * (3 * x_sqr_y_sqr + 2 * z*z)
        den = 3 * x_sqr_y_sqr**2 * r_diff_norm**3
        return num / den

    def I3():
        num = (
            (-2 * (x_sqr_y_sqr**2) * z2)
            + (x_sqr_y_sqr * (z**3 + 3 * (z2*z2) * z))
            - 2 * (z2*z2) * -(z**3)
        )
        den = 3 * x_sqr_y_sqr**2 * r_diff_norm**3
        return num / den

    def lI1():
        return 1/x_sqr_y_sqr

    def lI2():
        return 2/(3*x_sqr_y_sqr**2)

    def lI3():
        return (x_sqr_y_sqr + 2*z2**2)/(3*x_sqr_y_sqr**2)
    G = np.zeros(3)
    r_diff = r1-r2
    r_diff_norm = np.linalg.norm(r_diff)
    x_sqr_y_sqr = np.linalg.norm(r_diff[:2])**2
    z = r_diff[2]
    z2 = r2[2]
    ili1 = I1() - lI1()
    ili2 = I2() - lI2()
    ili3 = I3() - lI3()
    xya = np.linalg.norm(r1[:2])**2 - a**2
    G[0] = 2 * r1[0] * ili1 - 3 * r_diff[0] * (xya * ili2 + ili3)
    G[1] = 2 * r1[1] * ili1 - 3 * r_diff[1] * (xya * ili2 + ili3)
    G[2] = (np.linalg.norm(r1) ** 2 - a**2) / (r_diff_norm ** 3)

    return G / (4 * np.pi * a)


def B_comp(r, values: np.array, points: np.array, areas: np.array):
    """
    computes
    !!! outdated function that uses python for computations!!!
    !!! use cpp_modiule.b_comp for 2000x better perfomance !!!
    Parameters
    ----------
    r : Coordinates or np.array
        point, where the field is computed
    values : np.array
    points : np.array
        points used to integrate the field

    Returns
    -------
    B : np.array
        magnetic field density in vector form

    """
    B = np.asarray([0.0, 0.0, 0.0], dtype=np.float64)
    for value, point, area in zip(values, points, areas):
        add = value * np.asarray(GreenBl(r, point)) * area
        B = B + add
    return B


def B_comp_map(r, B_map, method="allmap", grid=False, debug=False, change=False):
    def B_comp_map_same_maps(r, B_map, heavy=True, change=True):
        """
        быстрое вычисление B в случае, когда
        сетка с числам совпадает с сеткой интегрирования
        Args:
            r (Coordinates): DESCRIPTION.
            B_map (Grid): DESCRIPTION.

        Returns:
            B_comp_map(np.array, 3d)

        """
        B = np.asarray([0.0, 0.0, 0.0], dtype=np.float32)
        R = B_map.r
        n = B_map.num
        debugging_obj = []
        if heavy:
            for i in range(n):
                lat, lon = B_map.latlon[i]
                r_2 = Coordinates(R, lat, lon, latlon=True)
                B_l = B_map.find_value(lat, lon, index=i)
                add = B_l * np.asarray(GreenBl(r, r_2)) * B_map.area[i]
                B = B + add
                if debug:
                    debugging_obj.append(add)

        else:
            for coor, S in zip(B_map.latlon, B_map.area):
                lat, lon = coor
                r_2 = Coordinates(R, lat, lon, latlon=True)
                B_l = B_map.find_value_easy(lat, lon)
                add = B_l * np.asarray(GreenBl(r, r_2)) * S
                B = B + add
        if debug:
            return B, debugging_obj
        else:
            return B

    def B_comp_map_diff_maps(r, grid, B_map):
        """
        r класса Coordinates
        """
        B = np.asarray([0.0, 0.0, 0.0], dtype=np.float32)
        R = grid.r
        for coor, S in zip(grid.latlon, grid.area):
            lat, lon = coor
            r_2 = Coordinates(R, lat, lon, latlon=True)
            B_l = B_map.find_value(lat, lon)

            add = B_l * np.asarray(GreenBl(r, r_2)) * S
            B = B + add
        return B

    if method == "allmap":
        return B_comp_map_same_maps(r, B_map, change=change)
    elif method == "diffmaps":
        return B_comp_map_diff_maps(r, grid, B_map)


class cell:
    def __init__(self, center, hs, sizeType="deg"):
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
class Grid:
    """
    a class keeps a bunch of values with their Coordinates,
    works best with same-radius grid and equally spread latitudes-longitudes
    """

    def __init__(
        self, r, latitudes, longitudes, hs=False, area=False, uniformgrid=True
    ):
        if latitudes.size != longitudes.size:
            raise ValueError("longitudes size does not match latitudes")
        self.num = np.size(latitudes)
        self.values = np.zeros_like(latitudes)
        self.valuesvector = np.zeros((self.num, 3))
        if hasattr(r, "__len__"):
            self.r_array = r
            self.r = np.mean(r)
        else:
            self.r_array = np.full_like(self.values, r)
            self.r = r
        self.r = r
        self.lat = latitudes
        self.lon = longitudes
        self.latlon = list(map(list, zip(latitudes, longitudes)))
        if uniformgrid is False:
            self.cells = []
            self.area = False
            for lat, lon, size in zip(latitudes, longitudes, hs):
                center = Coordinates(r, lat, lon, latlon=True)
                self.cells.append(cell(center, size))
        elif hs is not False:
            h = hs / 2
            self.area = self.r**2 * np.abs(
                (
                    np.cos(np.radians(90 - self.lat) + h)
                    - np.cos(np.radians(90 - self.lat) - h)
                )
                * hs
            )
        elif area is not False:
            self.area = area

        self.coors_set = [Coordinates(self.r, *ll, latlon=True)
                          for ll in self.latlon]
        self.xyz = [cs.vector for cs in self.coors_set]
        self.progress = 0

    def set_value(
        self, value: float, lat=None, lon=None, vector=False, easy=True, index=False
    ):
        if index is not False:
            if vector:
                self.valuesvector[index] = value
            else:
                self.values[index] = value
        elif lat is None:
            raise ValueError
        elif easy is True:
            self.set_value_easy(value, lat, lon, vector=vector)
        else:
            coor, i = find_nearest(self.latlon, (lat, lon), R=self.r)
            if vector:
                self.valuesvector[i] = value
            else:
                self.values[i] = value

    def set_value_easy(self, value, lat, lon, vector=False):
        i = np.where(((np.array(self.latlon) == [lat, lon]).all(1)))
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
        ind = np.where(((np.array(self.latlon) == [lat, lon]).all(1)))
        return self.values[ind]

    def save_pkl(self, name=False, empty=None, vector=None):
        """
        saves the Grid to a .pkl file

        Parameters
        ----------
        name : str, optional
            name of the file, by default False
        """

        if empty is None:
            if self.progress == 0:
                empty = True
            else:
                empty = False

        if empty is True:
            folder = 'emptymaps'
        else:
            if vector is True:
                folder = 'vectormaps'
            elif vector is False:
                folder = 'Lmaps'
            elif np.count_nonzero(self.values) == 0:
                folder = 'vectormaps'
            else:
                folder = 'Lmaps'

        if name is False:
            name = f'B_map {self.r:.2}'
        with open(f'{folder}/{name}.pkl',
                  'wb') as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def dataframe(self):
        """
        Creates a pandas.DataFrame of the Grid

        Returns
        -------
            pandas.DataFrame
        """
        d = {
            "latlon": self.latlon,
            "lat": self.lat,
            "lon": self.lon,
            "r": np.full_like(self.lat, self.r),
            "xyz": [cs.vector for cs in self.coors_set],
            "area": self.area,
            "values": self.values,
            "valuesvector": self.valuesvector,
        }
        return pd.DataFrame(data=d)

    def save_csv(self, name=False):
        """
        saves the Grid to a .csv file

        Parameters
        ----------
        name, optional
            name of the file, by default False
        """
        if name is False:
            name = f"B_map {self.r:.2}"
        df = self.dataframe(self)
        df.to_csv(f"csv/{name}.csv", index=True)

    def progress1(self):
        """
        manually adds 1 to the progress counter
        """
        self.progress = self.progress + 1

    def progress_brute(self):
        """
        overrides whatever progress is saved based on the amount of values stored
        """
        self.progress = int(
            max(
                [np.count_nonzero(self.values), np.count_nonzero(
                    self.valuesvector) / 3]
            )
        )

    def add_fieldinfo(self, m, dipolepos):
        """
        adds information about the field parameters

        Parameters
        ----------
        m
            dipole moment
        dipolepos
            position of a dipole
        """
        self.m = m
        self.dipolepos = dipolepos

    def change_coors(self):
        """
        manually change coordinates (made for compatability reasons)
        """
        self.coors_set = [Coordinates(self.r, *ll, latlon=True)
                          for ll in self.latlon]


class grid(Grid):
    # a class that exists only for backward compatibility reasons
    pass


class Grid3D():
    def __init__(self, xs, ys, zs):
        self.values = np.zeros_like(xs)
        self.x = xs
        self.y = ys
        self.z = zs
        self.xyz = np.asarray([xs, ys, zs]).T
        self.r = np.linalg.norm(self.xyz, axis=1)
        self.num = np.shape(self.values)[0]

    def set_value(self, value, coor=None, ind=None):
        if coor is None:
            ind = np.where(((np.array(self.xyz) == coor).all(1)))
        self.values[ind] = value


def load_grid(df):
    r = df.r[0]
    grid = Grid(r, df.lat, df.lon, area=df.area)
    grid.set_allvalues(df.values)
    grid.set_allvalues(df.valuesvector, vector=True)
    return grid


def create_grid(latlim, lonlim, N, r=696340 * 1000, name=False):
    # N - количество ячеек на градус
    n = int((latlim[1] - latlim[0]) * N)
    m = int((lonlim[1] - lonlim[0]) * N)
    lats = np.linspace(latlim[0], latlim[1], num=n)
    lons = np.linspace(lonlim[0], lonlim[1], num=m)
    latitudes = np.repeat(lats, m)
    longitudes = np.tile(lons, n)
    L = np.asarray([latitudes, longitudes])
    latitudes, longitudes = L
    hs = (np.pi / 180) * (lonlim[1] - lonlim[0]) / m
    B_mapempty = Grid(r, latitudes, longitudes, hs)
    if name is not False:
        B_mapempty.save_pkl(name, empty=True)
    return B_mapempty


def create_3Dgrid_upon_2Dgrid(grid: Grid, N, spherical=True):
    if spherical == False:
        grid_x, grid_y, grid_z = grid.xyz
        grid_x_s, grid_y_s, grid_z_s = np.sort(
            grid_x), np.sort(grid_y), np.sort(grid_z)
        xs_unique = np.linspace(grid_x_s[0], grid_x_s[1], num=N)
        ys_unique = np.linspace(grid_y_s[0], grid_y_s[1], num=N)
        zs_unique = np.linspace(grid_z_s[0], grid_z_s[1], num=N)
        xs, ys, zs = np.meshgrid(xs_unique, ys_unique, zs_unique)
        xs, ys, zs = xs.flatten(), ys.flatten(), zs.flatten()
        basic_volume = ((xs_unique[1]-xs_unique[0]) *
                        (ys_unique[1]-ys_unique[0])*(zs_unique[1]-zs_unique[0]))

    return Grid3D(xs, ys, zs), basic_volume * 1e6


def create_3Dgrid_sph(latlim, lonlim, r_lim, N):
    rs_unique = np.linspace(r_lim[0], r_lim[1], num=N)
    lat_unique = np.linspace(latlim[0], latlim[1], num=N)
    lon_unique = np.linspace(lonlim[0], lonlim[1], num=N)
    xs = []
    ys = []
    zs = []
    for r in rs_unique:
        x = r * np.sin(lon_unique * np.pi / 180.) * \
            np.sin((90 - lat_unique) * np.pi / 180.)
        z = r * np.cos(lon_unique * np.pi / 180.) * \
            np.sin((90 - lat_unique) * np.pi / 180.)
        y = r * np.cos((90 - lat_unique) * np.pi / 180.)
        xs.append(x), ys.append(y), zs.append(z)
    basic_volume = (r_lim[1]**3 * np.cos(np.radians(90-latlim[1])) -
                    r_lim[0]**3*np.sin(np.radians(90-latlim[0]))) / (N**3)
    return Grid3D(xs, ys, zs), basic_volume * 1e6 * np.radians(lonlim[1]-lonlim[0])


class Magneticline:
    def __init__(self, initial_point, initial_value, step):
        if type(initial_point) == Coordinates:
            initial_point = initial_point.vector
        self.initial_point = initial_point
        self.initial_value = initial_value
        self.points = [initial_point]
        self.values = [initial_value]
        self.step = step
        self.progress = 0

    def add_value(self, m, dipolepos=[0, 0, 0], stoppoint=None):
        vec = np.asarray(self.values[-1])
        new_point = (vec * self.step / np.linalg.norm(vec)) + np.asarray(
            self.points[-1]
        )
        self.values.append(dipolebetter(
            new_point, m, dipolepos, returnxyz=True))
        self.points.append(new_point)
        self.progress = self.progress + 1

    def add_value_comp(self, values, points, areas, stoppoint=None, sign=+1):
        vec = np.asarray(self.values[-1]) * sign
        new_point = vec * self.step / np.linalg.norm(vec) + np.asarray(
            self.points[-1]
        )
        val = cpp.b_comp(new_point, values, points, areas)
        self.values.append(val)
        self.points.append(new_point)
        self.progress = self.progress + 1

    def line_by_length(self, func, length, *args):
        for i in range(length // self.step):
            self.add_value(self, func, *args)

    def line_by_stoppoint(self, func, stoppoint, maxsteps=1e3, *args):
        for i in range(maxsteps):
            self.add_value(self, func, *args)
            if self.points[-1].r < stoppoint:
                break

    def save_pkl(self, name=False):
        if name is False:
            name = f"{self.initial_point.r:.2}-{self.initial_value:.2}"
        with open(f"maglines/{name}.pkl", "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def dataframe(self):
        d = {"points": self.points, "values": self.values}
        return pd.DataFrame(data=d)

    def save_csv(self, name=False):
        df = self.dataframe(self)
        df.to_csv(f"{name}.csv")

    def brute_progress(self):
        self.progress = len(self.points)

    def change_step(self, step):
        self.step = step


if __name__ == "__main__":
    pass
