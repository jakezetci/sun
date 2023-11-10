# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np
import pickle
from coordinates import ll2pt, coordinates
from lib import Grid
import pandas as pd
from plots import config

dates = pd.date_range(start="2023-08-08 00:12:00",
                      freq="6H", periods=4).values
ts = pd.to_datetime(str(dates[0]))
d = ts.strftime('%Y.%m.%d %I:%M%p')
print(d)
print(type(dates[0]))

array = np.arange(1, 5)

fig, ax = config()
ax.plot(dates, array)
