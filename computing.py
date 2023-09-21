# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 21:22:29 2023

@author: cosbo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import mplstereonet
from coordinates import coordinates
from lib import B_comp, grid, dipolebetter
from plots import sphere
import os

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


def model_grid(B_map, dipole, dipolepos, vector=False, name=False,
               returnobj=False, checkpoints=721):
    """
    Model dipole grid computation with checkpoints

    Args:
        func (TYPE): DESCRIPTION.
        latlon (TYPE): DESCRIPTION.
        r (TYPE): DESCRIPTION.
        dipole (TYPE): DESCRIPTION.
        dipolepos (TYPE): DESCRIPTION.
        vector (TYPE, optional): DESCRIPTION. Defaults to False.
        name (TYPE, optional): DESCRIPTION. Defaults to False.
        returnobj (TYPE, optional): DESCRIPTION. Defaults to False.

    Returns:
        None.

    """

    def crtname():
        nme = f'{np.random.randint(1,999)}'
        if (os.path.exists(f'checkpoint2 {name}.pkl') or
                os.path.exists(f'checkpoint1 {name}.pkl') or
                os.path.exists(f'{name}.pkl')):
            return crtname()
        else:
            return nme

    if name is False:
        name = crtname()

    checkpointsnum = int(checkpoints)
    B_map.save_pkl(name=f'checkpoint2 {name}')
    B_map.save_pkl(name=f'checkpoint1 {name}')
    r = B_map.r
    for i in range(B_map.progress, B_map.num):
        lat, lon = B_map.latlon[i]
        r1 = coordinates(r, lat, lon, latlon=True)
        B_map.set_value(dipolebetter(r1, dipole, dipolepos,
                                     returnxyz=vector, returnBl=not vector),
                        lat, lon, index=i, vector=vector)
        B_map.progress1()
        if i % checkpointsnum == 0:
            os.remove('checkpoint1 {name}.pkl')
            os.replace('checkpoint2 {name}.pkl', 'checkpoint1 {name}.pkl')
            B_map.save_pkl(name='checkpoint2')
    B_map.save_pkl(name=name)
    if returnobj:
        return B_map


if __name__ == "__main__":
    pass
