# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:07:08 2023

@author: cosbo
"""

import math
import numpy as np


def ll2pt(lat, lon):
    """
    Latitude, longitude -> spherical phi, theta

    """
    if lon >= 0:
        return math.radians(lon), math.radians(90 - lat)
    else:
        return math.radians(360 + lon), math.radians(90 - lat)


def pt2ll(phi, theta):
    """
    Spherical phi, theta -> latitude, longitude

    """
    if phi <= math.pi:
        return math.degrees(math.pi / 2 - theta), math.degrees(phi)
    else:
        return math.degrees(math.pi / 2 - theta), math.degrees(phi - 2 * math.pi)


def xyz2ll(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / r)
    phi = math.copysign(1, y) * math.acos(x / math.sqrt(x**2 + y**2))
    return r, pt2ll(phi, theta)


def pt2xyz(phi, theta, r):
    x = r * math.sin(theta) * math.sin(phi)
    y = r * math.cos(theta)
    z = r * math.sin(theta) * math.cos(phi)
    return x, y, z


def ll2xyz(lat, lon, r):
    phi, theta = ll2pt(lat, lon)
    return pt2xyz(phi, theta, r)


def xyz2pt(x, y, z):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / r)
    phi = math.copysign(1, y) * math.acos(x / math.sqrt(x**2 + y**2))
    return r, phi, theta


def xyz2truexyz(x, y, z):
    return -x, z, y


def xy2ll(x, y, r=696340):
    theta = math.acos(y / r)
    phi = math.asin(x / (r * math.sin(theta)))
    return pt2ll(phi, theta)


def xyR2xyz(x, y, R):
    """
    ну это гениальная функция
    """
    z = R**2 - x**2 - y**2
    return x, y, z


class Coordinates:
    def __init__(self, r1, r2, r3, spherical=False, latlon=False):
        if spherical is True:
            self.r, self.phi, self.theta = r1, r2, r3
            self.x, self.y, self.z = pt2xyz(self.phi, self.theta, self.r)
            self.lat, self.lon = pt2ll(self.phi, self.theta)
        elif latlon is True:
            self.r, self.lat, self.lon = r1, r2, r3
            self.phi, self.theta = ll2pt(self.lat, self.lon)
            self.x, self.y, self.z = pt2xyz(self.phi, self.theta, self.r)
        else:
            self.x, self.y, self.z = r1, r2, r3
            self.r, self.phi, self.theta = xyz2pt(self.x, self.y, self.z)
            self.lat, self.lon = pt2ll(self.phi, self.theta)
        self.xt, self.yt, self.zt = xyz2truexyz(self.x, self.y, self.z)
        self.vector = np.asarray([self.x, self.y, self.z])

    def project(self):
        if math.cos(self.theta) == 1:
            r = 0
        else:
            r = math.sin(self.theta) / (1 - math.cos(self.theta))
        return np.array([r * math.cos(self.phi), r * math.sin(self.phi)])


if __name__ == "__main__":
    pass
