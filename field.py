# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:58:55 2023

@author: cosbo
"""

from coordinates import coordinates
import math


def B_dipole(r, B0=1000, R=696340, returnBl=True):
    """
    r класса coordinates
    """
    B_r = -2 * B0 * (R/r.r)**3 * math.cos(r.theta)
    B_theta = -B0 * (R/r.r)**3 * math.sin(r.theta)

    if returnBl is True:
        return ((B_r * math.sin(r.theta) +
                B_theta * math.cos(r.theta)) * math.cos(r.phi))
    else:
        return [B_r, 0, B_theta]


if __name__ == "__main__":
    pass