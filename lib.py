# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:37:22 2023

@author: Home-Pc
"""
import numpy as np
import scipy.integrate as integrate
import math

from coordinates import coordinates

def ll2pt(lat, lon):
    """
    Latitude, longitude -> spherical phi, theta

    """
    if lon >= 0:
        return math.radians(lon), math.radians(90-lat)
    else:
        return math.radians(360+lon), math.radians(90-lat)


def pt2ll(phi, theta):
    """
    Spherical phi, theta -> latitude, longitude

    """
    if phi <= math.pi:
        return math.degrees(math.pi/2 - theta), math.degrees(phi)
    else:
        return math.degrees(math.pi/2 - theta), math.degrees(phi - 2*math.pi)


def xyz2ll(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z/r)
    phi = math.copysign(1, y)*math.acos(x/math.sqrt(x**2 + y**2))
    return r, pt2ll(phi, theta)


def ll2xyz(lat, lon, r=696340):
    phi, theta = ll2pt(lat, lon)
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)
    return x, y, z


def pt2xyz(phi, theta, r=696340):
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)
    return x, y, z

def xyz2pt(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z/r)
    phi = math.copysign(1, y)*math.acos(x/math.sqrt(x**2 + y**2))
    return r, phi, theta

def distance_sphere(a, b, R):

    lat1, lon1 = a
    lat2, lon2 = b
    phi_1 = math.radians(lat1)
    phi_2 = math.radians(lat2)

    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    d = math.sin(delta_phi / 2.0) ** 2 + math.cos(phi_1) * math.cos(phi_2) * math.sin(delta_lambda / 2.0) ** 2
    c = 2 * math.atan2(math.sqrt(d), math.sqrt(1 - d))
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
    x, y, z = r1 - r2
    return z/((x**2 + y**2)*np.linalg.norm(r1-r2))


def I2(r1, r2):
    x, y, z = r1 - r2
    num = z * (3 * x**2 + 3 * y**2 + 2 * z**2)
    den = 3 * (x**2 + y**2)**2 * (np.linalg.norm(r1-r2))**3
    return num/den


def I3(r1, r2):
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    num = -2 * (x**2+y**2)**2 * z2 + (x**2+y**2)*(z**3 + 3 * z2**2 * z) - 2 * z2*z**3
    den = 3 * (x**2 + y**2)**2 * (np.linalg.norm(r1-r2))**3
    return num/den


def lI1(r1, r2):
    x, y, z = r1 - r2
    return 1/(x**2 + y**2)


def lI2(r1, r2):
    x, y, z = r1 - r2
    return 1/(3 * (x**2 + y**2)**2)


def lI3(r1, r2):
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    return (x**2 + y**2 + 2 * z2**2)/(3 * (x**2 + y**2)**2)


def GreenBl(r1, r2, a=696340):
    r1, r2 = np.asarray([r1.x, r1.y, r1.z]), np.asarray([r2.x, r2.y, r2.z])
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    i1, i2, i3 = I1(r1, r2), I2(r1, r2), I3(r1, r2)
    li1, li2, li3 = lI1(r1, r2), lI2(r1, r2), lI3(r1, r2)
    G1 = 2*x1*(i1-li1) - 3 * x * ((x1**2 + y**2 - a**2)*(i2-li2)+(i3+li3))
    G2 = 2*y1*(i1-li1) - 3 * y * ((x1**2 + y**2 - a**2)*(i2-li2)+(i3+li3))
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


class grid:
    def __init__(self, r, latitudes, longitudes, hs,
                 uniformgrid=True):
        if latitudes.size != longitudes.size:
            raise ValueError('longitudes size does not match latitudes')
        self.num = np.size(latitudes)
        self.values = np.zeros_like(latitudes)
        if uniformgrid is False:
            self.cells = []
            self.area = False
            for lat, lon, size in zip(latitudes, longitudes, hs):
                center = coordinates(r, lat, lon, latlon=True)
                self.cells.append(cell(center, size))
        self.lat = latitudes
        self.lon = longitudes
        self.latlon = list(map(list, zip(latitudes, longitudes)))

        self.area = np.full_like(latitudes, ((hs * r * 2)**2))
        self.r = r


    def set_value(self, value, lat, lon):
        aa = np.asarray(self.latlon)
        bb = np.asarray(aa == [lat,lon])
        cc = list(bb)
        dd = [c[0] * c[1] for c in cc]
        i = np.argwhere(dd)[0]
        self.values[i] = value

    def find_value(self, lat, lon):
        tt = list(map(list, self.latlon))
        coor, ind = find_nearest(self.latlon, (lat, lon), R=self.r)
        return self.values[ind]

def B_comp(r, grid, B_map):
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


if __name__ == "__main__":
    pass
