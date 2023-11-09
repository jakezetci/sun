# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np
import pickle
from coordinates import ll2pt

x = np.arange(6).reshape(2, 3)
a = np.argwhere(x > 1)
add = np.full_like(a, [1, 2])
b = a + add
print(np.vstack([a, add]))

array = [1, 2]
array2 = [1, 2, 3]

array.append(array2)
print(array)
