# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np
import pickle
from coordinates import ll2pt

phi, theta = ll2pt(60,30)

alpha = np.arctan(np.cos(theta)/(np.sin(theta)*np.sin(phi)))
alphadeg = np.rad2deg(alpha)