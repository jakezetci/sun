# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 17:28:12 2023

@author: cosbo
"""

import plots
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

density = 4.7
day = 12
dates = pd.read_csv(
    f'dates_11158_{day}_dens{density}_fixed.txt', header=None).values.T[0]


energys = np.loadtxt('energys_11158_density=4.7, day=12.txt')/45


fig, ax = plots.config(
    ylabel='Энергия [эрг]', xlabel='', title='Энергия в активной области NOAA_AR=11158',
    figsize=(9, 6), dpi=240)

labels = ['12/02/2011', '13/02/2011', '14/02/2011', '15/02/2011']
x = [*dates[::6], '2011-02-15T00:00']
ax.plot(dates, energys, 'o', color='#0E2D69')
ax.plot(dates, energys, '-', color='#0E2D69')
plt.xticks(x, labels)
