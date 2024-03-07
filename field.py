# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:58:55 2023

@author: cosbo
"""

from coordinates import Coordinates

import math
import numpy as np


def B_dipole(r, B0=1000, R=696340 * 1e3, returnBl=True, returnxyz=False):
    """
    r класса Coordinates
    """
    B_r = -2 * B0 * (R / r.r) ** 3 * math.cos(r.theta)
    B_theta = -B0 * (R / r.r) ** 3 * math.sin(r.theta)

    if returnBl is True:
        return (B_r * math.sin(r.theta) + B_theta * math.cos(r.theta)) * math.cos(r.phi)
    elif returnxyz is True:
        return [
            B_r * math.sin(r.theta) * math.sin(r.phi)
            + B_theta * math.cos(r.theta) * math.sin(r.phi),
            B_r * math.cos(r.theta) + B_theta * math.sin(r.theta),
            B_r * math.sin(r.theta) * math.cos(r.phi)
            + B_theta * math.cos(r.theta) * math.cos(r.phi),
        ]
    else:
        return [B_r, 0, B_theta]


def dipolebetter(
    r,
    dipolemoment,
    rdipole=np.array([0, 0, 0]),
    returnBl=False,
    returnxyz=False,
    mu=1.25e-6,
):
    """
    r класса Coordinates, m вектор (x, y, z) в системе садыкова
    """
    if type(r) == Coordinates:
        r = r.vector
    r = Coordinates(*(r - rdipole))
    c = mu / (4 * math.pi * r.r**5)
    m = np.asarray(dipolemoment)
    rm = np.inner(r.vector, m)

    r_part = c * 3 * rm * r.vector
    m_part = -c * m * r.r**2
    L = np.asarray([0, 0, 1])
    B = r_part + m_part
    if returnBl is True:
        return np.inner(r_part, L) + np.inner(m_part, L)
    elif returnxyz is True:
        return r_part + m_part
    else:
        raise ValueError


def twomonopoles(r, q1, q2, rq1, rq2, mu=1.25e-6, returnxyz=False, returnBl=False):
    def force(r, q):
        r = np.asarray(r)
        c = mu / (4 * math.pi * np.linalg.norm(r) ** 2)
        return r * c

    rq1, rq2 = np.asarray(rq1), np.asarray(rq2)
    force1 = force(r.vector - rq1, q1)
    force2 = force(r.vector - rq2, q2)
    L = np.asarray([0, 0, 1])
    force_sum = force1 + force2
    if returnBl is True:
        return np.dot(force_sum, L)
    elif returnxyz is True:
        return force_sum
    else:
        raise ValueError


if __name__ == "__main__":
    pass
