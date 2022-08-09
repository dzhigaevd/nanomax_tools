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

def get_scan_parameters(command):
    """
    Parse the command from the beamline and return number of points in the scan, orientation of the scan, and range in real space
    """
    coordinate_to_direction = {
    'sx': 'h',
    'sy': 'v',
    'sz': 'b'}
    # Parse scan command assuming it is the same from start_scan_number to end_scan_number: should be correct
    if command[0] == 'npointflyscan':
        fast_axis = command[1]
        fast_axis_range = [float(command[2]),float(command[3])]
        fast_axis_points = int(command[4])+1

        slow_axis = command[5]
        slow_axis_range = [float(command[6]),float(command[7])]
        slow_axis_points = int(command[8])+1
    else:
        print('Unknown command!')
    
    scan_orientation = coordinate_to_direction[fast_axis]+coordinate_to_direction[slow_axis]
    
    return fast_axis_points, slow_axis_points, fast_axis_range, slow_axis_range, scan_orientation

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
        
    h5file.close()
    return command, motor_positions, rocking_motor, rocking_angles, scan_position_x, scan_position_y, scan_position_z, incoming_intensity

def read_data_merlin(data_path,roi_xrd=None,point_list=None):
    h5file = h5py.File(data_path, 'r')
    if roi_xrd!=[]:
        if point_list!=[]:
            data = h5file['entry']['measurement']['merlin']['frames'][point_list,roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3]]       
        else:
            data = h5file['entry']['measurement']['merlin']['frames'][:,roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3]]       
    else:
        data = h5file['entry']['measurement']['merlin']['frames'][()]
    
    h5file.close()
    return data

def read_data_xspress3(data_path,roi=None):
    module = 3 # hardcoded for 4.04.2022
    h5file = h5py.File(data_path, 'r+')
    if roi:
        data = h5file['entry']['measurement']['xspress3']['frames'][:,module,roi[0]:roi[1]]
        data = np.sum(data,1)
    else:
        data = h5file['entry']['measurement']['xspress3']['frames'][:,module,:]
    print('Use XRF module #'+str(module))
    h5file.close()
    
    return data

def read_mask(data_path,roi):
    data = np.load(data_path)
    mask = data['mask'][roi[0]:roi[1],roi[2]:roi[3]]
    return mask
