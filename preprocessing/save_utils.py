#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 22:29:50 2020

@author: dzhigd
"""
import os
import scipy.io as scio
import numpy as np
from tifffile import imsave

def save_data_xrd(save_path, save_type, \
	data_xrd, xrd_roi, motor_positions, \
	rocking_motor,	rocking_angles, \
	scan_position_x, scan_position_y, scan_position_z):
	
	if save_type == "mat":
	    # Permute arrays, to have the rocking direction as a last index
	    data_xrd = np.transpose(np.single(data_xrd),(1,2,0))	    
	    data_dic = {"data": data_xrd,\
	                "command": command, \
	                "motor_positions":motor_positions,\
	                "rocking_motor":rocking_motor, \
	                "rocking_angles":rocking_angles}
	    scio.savemat(save_path, {'scan': data_dic}, do_compression=True, oned_as="column")
	    print(('++ Saved XRD and meta data to %s')%(save_path))

	elif save_type == "npz": 
		np.savez_compressed(save_path, data_xrd = data_xrd, \
    							xrd_roi = xrd_roi, \
    							motor_positions = motor_positions, \
    							rocking_motor = rocking_motor,\
    							rocking_angles = rocking_angles, \
    							scan_position_x = scan_position_x, \
    							scan_position_y = scan_position_y, \
    							scan_position_z = scan_position_z)
		print(('++ Saved XRD and meta data to %s')%(save_path))
	else:
		error('-- Unrecognized save type! Do nothing ...')

def save_data_xrf(save_path, data_xrf):
	data_xrf = np.float32(data_xrf/np.max(data_xrf))
	imsave(save_path, data_xrf)
	print(('++ Saved XRF data to %s')%(save_path))

def convert_to_single_npz(path, scan_numbers):
	# Load all 

def save_bin_data_nanomax(save_path,data):
    fid = open(save_path,'w+b')
    fid.write(data)
    fid.close()

