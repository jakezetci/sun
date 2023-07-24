# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 22:37:22 2023

@author: Home-Pc
"""
import numpy as np
import scipy.integrate as integrate
from numpy import sqrt
import sunpy.map


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
    return (x**2 + y**2 +2*z2**2)/(3*(x**2 + y**2)**2)

def GreenBl(r1, r2, a=696340):
    x, y, z = r1 - r2
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    i1, i2, i3 = I1(r1, r2), I2(r1, r2), I3(r1, r2)
    li1, li2, li3 = lI1(r1, r2), lI2(r1, r2), (r1, r2)
    G1 = 2*x1*(i1-li1) - 3 * x * ((x1**2 + y**2 - a**2)*(i2-li2)+(i3+li3))
    G2 = 2*y1*(i1-li1) - 3 * y * ((x1**2 + y**2 - a**2)*(i2-li2)+(i3+li3))
    G3 = (np.linalg.norm(r1)**2 - a**2) / (np.linalg.norm(r1-r2)**3)
    return np.array([G1, G2, G3])/(4*np.pi*a)

class grid:
    
    def __init__(self, matrix):
        self.num = np.size[0]

def B_comp(r, grid):
    