# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np
import pickle

b = 'any'
with open(f'Lmaps/{b}.pkl', 'wb') as f:
    pickle.dump(b, f, pickle.HIGHEST_PROTOCOL)