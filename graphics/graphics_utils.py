#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 23:43:44 2020

@author: dzhigd
"""
from scipy import ndimage, interpolate
import numpy as np
import matplotlib.pyplot as plt

def show_q_space_projections(qx,qy,qz,data):    
#     Constants. They are needed for correct labeling of axes
    h                       = 4.1357e-15;                                  # Plank's constant
    c                       = 2.99792458e8;                                # Speed of light in vacuum
    energy = 11000
    pitch = 55e-6
    radius = 0.4
    
    wavelength = h*c/energy

    dq = (2*np.pi*pitch/(radius*wavelength))
    
    qxv = np.arange(np.min(qx),np.max(qx),dq)
    qyv = np.arange(np.min(qy),np.max(qy),dq)
    qzv = np.arange(np.min(qz),np.max(qz),dq)
    
    Qx,Qy,Qz = np.meshgrid(qxv, qyv, qzv)
    
    # find interpolant
    q = np.array((qz,qx,qy)).transpose()
    F = interpolate.LinearNDInterpolator(q, np.concatenate(np.concatenate(data)), fill_value=np.nan, rescale=True)
    
    # interpolate, this is the output of aligned diffraction data
    qspace_interpolated = np.nan_to_num(F(Qz.transpose(),Qx.transpose(),Qy.transpose()))
    
#    qspace_interpolated[qspace_interpolated == np.nan] = 0
    
    # sum data to plot total xrd map    
    plt.figure(num=1)
    plt.imshow(np.sum(qspace_interpolated,axis=0))
    
    plt.figure(num=2)
    plt.imshow(np.log10(np.sum(qspace_interpolated,axis=1)))
    
    plt.figure(num=3)
    plt.imshow(np.sum(qspace_interpolated,axis=2))