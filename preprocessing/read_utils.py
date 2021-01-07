#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:40:33 2020

@author: dzhigd
"""
import os
import hdf5plugin
import h5py
import numpy as np

def read_data_meta(path):
    h5file = h5py.File(path,'r')
    command = str(h5file['entry']['description'][()])[3:-2] # Reading only useful symbols
    motor_positions = {
            # Detector positions
            "delta":    h5file['entry']['snapshot']['delta'][()],
            "gamma":    h5file['entry']['snapshot']['gamma'][()],           
            "gonphi":   h5file['entry']['snapshot']['gonphi'][()],
            "gontheta": h5file['entry']['snapshot']['gontheta'][()],
            "radius":   h5file['entry']['snapshot']['radius'][()],
            "energy":   h5file['entry']['snapshot']['energy'][()]
            }
    
    try:
        scan_position_x = h5file['entry']['measurement']['pseudo']['x'][()]
    except:
        print("-- Current dataset has no x lateral scanning, continue with single position")
        scan_position_x = []
    try:
        scan_position_y = h5file['entry']['measurement']['pseudo']['y'][()]
    except:
        print("-- Current dataset has no y lateral scanning, continue with single position")
        scan_position_y = []
    try:
        scan_position_z = h5file['entry']['measurement']['pseudo']['z'][()]
    except:
        print("-- Current dataset has no z lateral scanning, continue with single position")
        scan_position_z = []
    
    try:
        incoming_intensity = h5file['entry']['measurement']['alba2']['1'][()]
    except:
        print("-- Normalization data is not found, consider manual normalization, continue without it")
        incoming_intensity = []
    
    try:
        rocking_motor = "gonphi"
        rocking_angles = h5file['entry']['measurement'][rocking_motor][()]
        print("-- Rocking motor is %s --"%rocking_motor)
    except:
        try:
            rocking_motor = "gontheta"
            rocking_angles = h5file['entry']['measurement'][rocking_motor][()]
            print("-- Rocking motor is %s"%rocking_motor)
        except:
            print("-- No rocking motor positions, pass or specify it separately!")
            rocking_angles = []
            rocking_motor = []
            pass

    return command, motor_positions, rocking_motor, rocking_angles, scan_position_x, scan_position_y, scan_position_z, incoming_intensity

def read_data_merlin(data_path,roi):
    h5file = h5py.File(data_path, 'r')
    data = h5file['entry']['measurement']['Merlin']['data'][:,roi[0]:roi[1],roi[2]:roi[3]]       
    return data

def read_data_xspress3(data_path,roi):
    h5file = h5py.File(data_path, 'r')
    entry_keys = list(h5file['/'])
    ii = 0
    nPoints = len(h5file[entry_keys[0]]['measurement']['xspress3']['data'][()])
    data = np.zeros((nPoints, len(entry_keys)))
    for key in entry_keys:
        data_temp = h5file[key]['measurement']['xspress3']['data'][:,3,roi[0]:roi[1]]
        data[:,ii] = np.sum(data_temp,1)       
        ii = ii+1
    return data

def read_mask(data_path,roi):
    data = np.load(data_path)
    mask = data['mask'][roi[0]:roi[1],roi[2]:roi[3]]
    return mask
