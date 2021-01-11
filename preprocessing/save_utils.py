#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 22:29:50 2020

@author: dzhigd
"""
import os
import scipy.io as scio
import numpy as np

def create_save_path(processing_pre_path, sample_name, scan_number, save_type, ii):
	try:
        processing_path = os.path.join(processing_pre_path,("%s_%d_%d")%(sample_name,scan_number[0],scan_number[-1]))
    except:
        processing_path = os.path.join(processing_pre_path,("%s_%d")%(sample_name,scan_number[ii]))    

    try:
        os.mkdir(processing_path)
        print("++ Processing folder created")
    except:
        print("-- Processing folder already exists")

	if save_type == "mat":
		# Save mat compressed array for matlab                 
    	save_path = os.path.join(processing_path, ("scan_%06d_merlin.mat")%(int(scan_number[ii])))    #The path to save the results
 	elif save_type == "npz":
 	   	# Save mat compressed array for matlab                 
    	save_path = os.path.join(processing_path, ("scan_%06d_merlin.npz")%(int(scan_number[ii])))    #The path to save the results                        
    
	return save_path

def save_dara_xrd_mat(processing_path, data_xrd, xrd_roi, motor_positions, rocking_motor,\
     					rocking_angles, scan_position_x, scan_position_y, scan_position_z,ii):
	
    # Permute arrays, to have the rocking direction as a last index
    data_xrd = np.transpose(np.single(data_xrd),(1,2,0))
    
    data_dic = {"data": data_xrd,\
                "command": command, "motor_positions":motor_positions,\
                "rocking_motor":rocking_motor, "rocking_angles":rocking_angles}

    scio.savemat(save_path, {'scan': data_dic}, do_compression=True,oned_as="column") 
    print(('++ Saved scan %d to %s')%(int(scan_number[ii]),save_path))

def save_data_xrd_npz(processing_path, data_xrd, xrd_roi, motor_positions, rocking_motor,\
     					rocking_angles, scan_position_x, scan_position_y, scan_position_z,ii):	
    
    
    np.savez_compressed(save_path, data_xrd = data_xrd, xrd_roi = xrd_roi, motor_positions = motor_positions, rocking_motor = rocking_motor,\
     rocking_angles = rocking_angles, scan_position_x = scan_position_x, scan_position_y = scan_position_y, scan_position_z = scan_position_z)
    print(('++ Saved XRD scan %d to %s')%(int(scan_number[ii]),save_path))

def save_data_xrf(data_xrf,processing_path):
	data_xrf = np.float32(data_xrf/np.max(data_xrf))
    save_path = os.path.join(processing_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii])))
    imsave(save_path, data_xrf)
    print(('++ Saved XRF scan %d to %s')%(int(scan_number[ii]),save_path))

def convert_to_single_npz(path, scan_numbers):


def save_bin_data_nanomax(save_path,data):
    fid = open(save_path,'w+b')
    fid.write(data)
    fid.close()

