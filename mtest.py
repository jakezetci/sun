# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 00:26:45 2023

@author: cosbo
"""

import numpy as np


def a(y, u):
    print(y, u)

def b(x,y,z):
    print(x,y,z)

def what(func, *args):
    func(*args)



what(b, 8, 9, 10)