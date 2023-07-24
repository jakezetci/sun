# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 23:34:15 2023

@author: Home-Pc
"""

import numpy as np

r = [3, 6, 4]
x,y,z = r
print(np.sqrt(x**2 + y**2 + z**2))
print(np.linalg.norm(r))