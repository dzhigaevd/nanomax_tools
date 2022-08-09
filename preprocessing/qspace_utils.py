#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 22:00:25 2020

@author: dzhigd
"""

# Select the rocking curve

from scipy.interpolate import griddata as gd

def get_qspace_coordinates(xrd_data):
