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

with open('callid.txt') as f:
    lines = f.read()
lines = int(lines)
