# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:58:55 2023

@author: cosbo
"""

from coordinates import coordinates
import math
import numpy as np


def B_dipole(r, B0=1000, R=696340, returnBl=True, returnxyz=False):
    """
    r класса coordinates
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


def dipolebetter(r, m, rdipole=np.array([0, 0, 0]), returnBl=False,
                 returnxyz=False, mu=1.25e-6):
    """
    r класса coordinates, m вектор (x, y, z) в системе садыкова
    """
    r = coordinates(*(r.vector-rdipole))
    c = mu/(4 * math.pi * r.r**3)
    m = np.asarray(m)
    rm = np.dot(r.vector, m)

    r_part = c * 3 * rm * r.vector
    m_part = - c * m
    L = np.asarray([0, 0, 1])
    if returnBl is True:
        return np.dot(r_part, L) + np.dot(m_part, L)
    elif returnxyz is True:
        return r_part + m_part
    else:
        raise ValueError


if __name__ == "__main__":
    pass
