# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:52:18 2023

@author: cosbo
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from coordinates import coordinates
from lib import B_comp_map, Grid
from field import B_dipole
import numpy as np
import math
from plots import sphere, disk
import pickle
from plots import config


if __name__ == "__main__":
    """
    with open('maglinemodel.pkl', 'rb') as fff:
        line1 = pickle.load(fff)
    with open('maglinecomp2.pkl', 'rb') as fff:
        line2 = pickle.load(fff)
    fig2, ax3 = plt.subplots(1, 1)
    ax3.set_title("расчёты...", fontsize=20, fontname='HSE Slab')
    x1, y1 = np.array(line1.pointsxyz)[:,0], np.array(line1.pointsxyz)[:,1]
    points2 = np.array(line2.pointsxyz)[:,:2]


    x2, y2 = np.array(line2.pointsxyz)[:,0], np.array(line2.pointsxyz)[:,1]
    ax3.plot(x1, y1,'o', lw=0.5)
    ax3.plot(x2, y2, 'o', lw=0.5)
    """
    pass