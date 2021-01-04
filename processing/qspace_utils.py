#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 23:51:17 2020

@author: dzhigd
"""

import numpy as np


# Q-space #############################################################
def load_q_coordinates(path):    
    q = np.load(path)
    qx = q['qx']
    qy = q['qy']
    qz = q['qz']
    return qx,qy,qz

def get_q_coordinates(xrd_roi,data_size,motor_positions,beamline): 
    # Constants. They are needed for correct labeling of axes
    h = 4.1357e-15                                                            # Plank's constant
    c = 2.99792458e8                                                          # Speed of light in vacuum
    
    wavelength = h*c/motor_positions['energy']
    
    k = 2*np.pi/wavelength; # wave vector
    
    dq = k*2*np.arctan(beamline['detector_pitch']/(2*nanomax['radius']))      # q-space pitch at the detector plane
    
    hd,vd = np.meshgrid(np.linspace(-data_size[1]/2,(data_size[1]/2-1),data_size[1]),np.linspace(-data_size[0]/2,(data_size[0]/2-1),data_size[0]))
    
    hd = (hd+(data_size[1]/2-beamline['direct_beam'][1]))*beamline['detector_pitch'];
    vd = (vd+(data_size[0]/2-beamline['direct_beam'][0]))*beamline['detector_pitch'];
    zd = np.ones(size(vd))*beamline['radius']
    
    hd = hd[xrd_roi(2):xrd_roi[3],xrd_roi[0]:xrd_roi[1],:];
    vd = vd[scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:];
    zd = zd[scan.roi(2):scan.roi(2)+scan.roi(4),scan.roi(1):scan.roi(1)+scan.roi(3),:];

    d = [np.concatenate(hd),np.concatenate(vd),np.concatenate(zd)]
    d = np.transpose(d)

    r = np.squeeze(np.sqrt(np.sum(d**2,0)));

    hq = k*(d[1,:]/r);
    vq = k*(d[2,:]/r);
    zq = k*(1-d[3,:]/r);

    q = [hq,vq,zq];
    
    return q