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

def save_data_xrd(save_path, save_type, save_dictionary):
	if save_type == "mat":	        	  
	    scio.savemat(save_path, {'scan': save_dictionary}, do_compression=True, oned_as="column")
	    print(('++ Saved XRD and meta data to %s')%(save_path))
	elif save_type == "npz": 
  		np.savez_compressed(save_path, save_dictionary)
  		print(('++ Saved XRD and meta data to %s')%(save_path))
	else:
		error('-- Unrecognized save type! Do nothing ...')

def save_data_xrf(save_path, data_xrf):
	data_xrf = np.float32(data_xrf/np.max(data_xrf))
	imsave(save_path, data_xrf)
	print(('++ Saved XRF data to %s')%(save_path))

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
    elif save_type == "tif":
        # Save in the form of tiff file
        save_path = os.path.join(processing_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii])))
    return save_path