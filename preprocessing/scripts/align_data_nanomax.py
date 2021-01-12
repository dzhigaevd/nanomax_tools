#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dmitry dzhigaev and tomas stankevic

As an input a corrected and converted dataset is required
This script is adapted for 5D datasets of scanning XRD
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
from tifffile import imsave
from tifffile import imread
import tifffile as tf
from scipy import ndimage, interpolate
from skimage.transform import resize
import matplotlib.mlab as ml
from scipy.optimize import minimize
import sys
sys.path.append('/home/dzhigd/Software/')
import nanomax_tools.preprocessing.align_utils as au
import nanomax_tools.preprocessing.read_utils as ru
import nanomax_tools.preprocessing.save_utils as su

def create_save_path(save_pre_path, scan_number, save_type, ii):  
    if save_type == "mat":
        save_path = os.path.join(save_pre_path, ("scan_%06d_merlin.mat")%(int(scan_number[ii])))    #The path to save the results        
    elif save_type == "npz":
        # Save mat compressed array for matlab                 
        save_path = os.path.join(save_pre_path, ("scan_%06d_merlin.npz")%(int(scan_number[ii])))    #The path to save the results
    elif save_type == "tif":
        # Save in the form of tiff file
        save_path    = os.path.join(save_pre_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii])))
    return save_path

# INPUTS ######################################################################
year        = "2020"                                                           #The year for the experiemnt
beamtimeID  = "2020101408"                                                     #The beamtimeID
sample_name = r"0002_sample_P246_AB"                                           #The name for the p10 newfile
scan_number = np.linspace(109,128,20)                                          #The scan numbers, can be the list
rocking_motor = "gontheta" # "gonphi"

align_scans      = True
align_data_xrd   = True
save_data_interpolate = True
save_type        = ['mat']

xrd_map_crop = 40 # Start saving the map from this position
scan_n_z = 101
scan_n_x = 91

scan_pixel_size = [0.08,0.08] # um from the motor positions in the data
skip_pixel_number = 1
###############################################################################

#determing the paths"
processing_pre_path = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/process/%s_%d_%d"%(beamtimeID,sample_name,scan_number[0],scan_number[-1])
reference_path      = r"/home/dzhigd/work/projects/Qdevs_2020_NanoMAX/data/reference_%d_%d.tif"%(scan_number[0],scan_number[-1])
save_pre_path       = os.path.join(processing_pre_path,'data_aligned')

try:
    os.mkdir(save_pre_path)
    print("++ Save folder created")
except:
    print("-- Save folder already exists")

orig_pixel_size = [0.01,0.01] # um pixel size of synthetic image
offset_top_left = [-2.7,0]

with tf.TiffFile(reference_path) as tif:
            reference_image = tif.asarray()

for ii in range(0,len(scan_number)):
    print(('\033[1;33;40m #### Processing scan %d')%(int(scan_number[ii])))

    #determing the paths"
    data_npz_path = os.path.join(processing_pre_path, "scan_%06d_merlin.npz"%(int(scan_number[ii])))
    data_xrf_path = os.path.join(processing_pre_path, "scan_%06d_xspress3.tif")%(int(scan_number[ii]))   

    # Load the data ###########################################################
    data = np.load(data_npz_path,allow_pickle=True)    
    dictionary = data['arr_0'].item()

    data_xrd = dictionary['data_xrd']
    xrd_roi  = dictionary['xrd_roi']
    scan_position_x = dictionary['scan_position_x']
    scan_position_y = dictionary['scan_position_y']
    scan_position_z = dictionary['scan_position_z']
    motor_positions = dictionary['motor_positions']
    rocking_angles = dictionary['rocking_angles']

    scan_positions = np.array([scan_position_x,scan_position_z]) # IMPORTANT TO ADJUST ACCORDING TO THE SCAN MOTORS!

    with tf.TiffFile(data_xrf_path) as tif:
            data_xrf = tif.asarray()
    ###########################################################################

    if align_scans == True:                    
        # Do for each angle:
        if ii == 0:
            reference_image,mask_scan = au.generate_ref_image(data_xrf,reference_image,scan_pixel_size,orig_pixel_size,offset_top_left)
        
        # find misalignments        
        fluo_aligned, X_shift, Y_shift = au.align_image(data_xrf, reference_image, scan_positions, scan_pixel_size)
        
        if align_data_xrd == True:
            if ii == 0:
                fig_xrf, axs_xrf = plt.subplots(int(np.ceil(len(scan_number)/5)),5, figsize=(15, 6), facecolor='w', edgecolor='k')
                fig_xrf.subplots_adjust(hspace = .2, wspace=.001)
                fig_xrf.suptitle('Total XRF maps')
                axs_xrf = axs_xrf.ravel()

                fig_xrd, axs_xrd = plt.subplots(int(np.ceil(len(scan_number)/5)),5, figsize=(15, 6), facecolor='w', edgecolor='k')
                fig_xrd.subplots_adjust(hspace = .2, wspace=.001)
                fig_xrd.suptitle('Total XRD maps')
                axs_xrd = axs_xrd.ravel()

            # interpolate the misalignments from regular grid onto the motor positions so they can be just added
            X_shift = au.grid_to_image(X_shift,scan_positions,X_shift)
            
            # add misalignments to positions
            positions_new = scan_positions + np.array([Y_shift.transpose().ravel(),X_shift.transpose().ravel()])
            
            # find interpolant
            F = interpolate.LinearNDInterpolator(scan_positions.transpose(), data_xrd, fill_value=np.nan, rescale=True)
            
            # interpolate, this is the output of aligned diffraction data
            xrd_interp = F(positions_new.transpose())
            
            # plot total xrf map            
            axs_xrf[ii].imshow(fluo_aligned)
            axs_xrf[ii].set_title('Scan %d'%int(scan_number[ii]))
            plt.show(block=False)

            # sum data to plot total xrd map
            xrd_int = np.sum(np.sum(xrd_interp,axis=2),axis=1)
            xrd_int.shape
            xrd_image = np.reshape(xrd_int,data_xrf.transpose().shape).transpose()
            
            axs_xrd[ii].imshow(xrd_image)
            axs_xrd[ii].set_title('Scan %d'%int(scan_number[ii]))
            plt.show(block=False)

        # Save the interpolated xrd data ######################################
        # Save the original xrd data separately for each angle #################### 
        if save_data_interpolate == True:                  
            for value in save_type:
                save_path  = create_save_path(save_pre_path, scan_number, value, ii)                  
                if value == 'mat' or value == 'npz':
                    xrd_interp = np.reshape(xrd_interp,[scan_n_z,scan_n_x,xrd_roi[1]-xrd_roi[0],xrd_roi[3]-xrd_roi[2]])   
                    save_dictionary = {"data_xrd": xrd_interp[:,xrd_map_crop:None,:,:], \
                                        "xrd_roi": xrd_roi, \
                                        "motor_positions": motor_positions, \
                                        "rocking_motor": rocking_motor, \
                                        "rocking_angles": rocking_angles, \
                                        "scan_position_x": scan_position_x[0:xrd_map_crop], \
                                        "scan_position_y": scan_position_y, \
                                        "scan_position_z": scan_position_z[0:xrd_map_crop]}  
                if value == 'tif':
                    data_xrf = np.float32(fluo_aligned/np.max(fluo_aligned))
                    save_dictionary = data_xrf
                
                su.save_data(save_path, value, save_dictionary)
    
    print(('\033[1;32;40m #### Processing scan %d done!')%(int(scan_number[ii])))

plt.show()
print("\033[1;37;40m #### Overall Processing Done!")