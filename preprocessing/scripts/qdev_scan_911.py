#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dmitry dzhigaev and tomas stankevic
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
from tifffile import imsave
import tifffile as tf
from scipy import ndimage, interpolate
from skimage.transform import resize
import matplotlib.mlab as ml
from scipy.optimize import minimize
import sys
sys.path.append('/home/dzhigd/Software/nanomax_tools')
import nanomax_tools.preprocessing.map_alignment as ma
import nanomax_tools.preprocessing.read_utils as ru
    
# INPUTS ######################################################################
year        = "2020"                                                           #The year for the experiemnt
beamtimeID  ="2020101408"                                                      #The beamtimeID
sample_name = r"0007_MQML258BC"                                           #The name for the p10 newfile
#scan_number = np.linspace(109,128,20)                                          #The scan numbers, can be the list
#scan_number = np.linspace(137,156,20)                                          #The scan numbers, can be the list
scan_number = [911] # 146
rocking_motor = "gontheta" # "gonphi"

# Flags of processing options
normalize        = True
apply_mask       = True
xrd_analysis     = True
xrf_analysis     = False
q_space_analysis = False
save_data_xrd    = True
save_data_xrf    = False
align_scans      = False
align_data_xrd   = False
save_data_xrd_interpolate = False
save_data_xrf_interpolate = False

xrf_roi     = [870,970]
#xrd_roi     = [0,220,250,400]
xrd_map_crop = 40 # Start saving the map from this position
xrd_roi     = [0,515,0,515]
scan_n_z = 101
scan_n_x = 91

scan_pixel_size = [0.08,0.08] # um from the motor positions in the data
skip_pixel_number = 1
###############################################################################

#determing the paths"
path_root       = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/%s"%(beamtimeID,sample_name)

if len(scan_number)>1:
    processing_path = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/process/scan_%d_%d"%(beamtimeID,scan_number[0],scan_number[-1])
else:
    processing_path = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/process/scan_%d"%(beamtimeID,scan_number[0])
    
mask_path       = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/mask11kev_102020.npz"%(beamtimeID)

#merlin_mask_190222_14keV.h5,mask11kev_102020.npz
try:
    os.mkdir(processing_path)
    print("--Processing folder created--")
except:
    print("--Processing folder already exists--")
    
# DATA IMPORT #################################################################
orig_pixel_size = [0.01,0.01] # um pixel size of synthetic image
rocking_angle = np.zeros((len(scan_number),1))
mask_t = ru.read_mask(mask_path, xrd_roi)

##### Start here
offset_top_left = [-2.7,0]

for ii in range(0,len(scan_number)):
    print(('--Processing scan %d--')%(int(scan_number[ii])))
    #determing the paths"
    
    data_meta_path = os.path.join(path_root, "%06d.h5"%(int(scan_number[ii])))
    data_xrd_path  = os.path.join(path_root, "scan_%06d_merlin.hdf5")%(int(scan_number[ii]))
    data_xrf_path  = os.path.join(path_root, "scan_%06d_xspress3.hdf5")%(int(scan_number[ii]))
    
    # Load the data ###########################################################
    [command, motor_positions, scan_positions_x, scan_positions_z, incoming_intensity] = ru.read_data_meta(data_meta_path) 
    scan_positions_xz = np.array([scan_positions_x,scan_positions_z])
    
    if xrd_analysis == True:
        data_xrd  = ru.read_data_merlin(data_xrd_path,xrd_roi)
        if apply_mask == True:
            if ii == 0:
                mask = np.tile(mask_t,(np.shape(data_xrd)[0],1,1))
            data_xrd = np.multiply(data_xrd,mask)
            print("--Data was masked--")
        if normalize == True:
            for jj in range(0,np.shape(data_xrd)[0]):
                data_xrd[jj,:,:] = data_xrd[jj,:,:]/incoming_intensity[jj]
            print("--XRD data was normalized--")
            
    if xrf_analysis == True:
        data_xrf  = ru.read_data_xspress3(data_xrf_path,xrf_roi)
        if normalize == True:
            for jj in range(0,np.shape(data_xrf)[0]):
                data_xrf[jj] = data_xrf[jj]/incoming_intensity[jj]
            print("--XRF data was normalized--")
                
    rocking_angle[ii] = motor_positions[rocking_motor]
    
    ###########################################################################
    
    # Load the data ###########################################################
#    figure(num=13)
#    imshow((np.sum(data_xrd,0)))
    
#    figure(num=14)
#    imshow(np.multiply(np.sum(data_xrd,0),mask))
#    plt.colorbar()
    
    # Save the original xrd data ############################################## 
    if save_data_xrd == True:
        save_path = os.path.join(processing_path, "scan_%06d_merlin_%d_%d_%d.npz")%(int(scan_number[ii]),(data_xrd.shape[0]),(data_xrd.shape[1]),(data_xrd.shape[2]))      #The path to save the results
    #    save_bin_data_nanomax(save_path, data_xrd)        
        np.savez_compressed(save_path,data_xrd = data_xrd, \
                    command = command, \
                    motor_positions = motor_positions, \
                    rocking_motor = rocking_angle, \
                    scan_positions_z = scan_positions_z)        
        print(('Saved scan %d to %s')%(int(scan_number[ii]),save_path))
    
    if save_data_xrf == True:        
        data_xrf = np.float32(data_xrf/np.max(data_xrf))   
        imsave(os.path.join(processing_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii]))), data_xrf)
        print(('Saved scan %d to %s')%(int(scan_number[ii]),processing_path))
    
    print("--Import Done!--")
    
    if align_scans == True:
        # make reference image   
        ref_path = r"/home/dzhigd/work/projects/Qdevs_2020_NanoMAX/data/reference_%d_%d.tif"%(scan_number[0],scan_number[-1])
        
        with tf.TiffFile(ref_path) as tif:
            ref_image = tif.asarray()
        # Do for each angle:
        if ii == 0:
            reference_image,mask_scan = ma.generate_ref_image(data_xrf,ref_image,scan_pixel_size,orig_pixel_size,offset_top_left)
        
        # find misalignments
        fluo_aligned, X_shift, Y_shift = ma.align_image(data_xrf, reference_image, scan_positions_xz, scan_pixel_size)
        
        if align_data_xrd == True:
            # interpolate the misalignments from regular grid onto the motor positions so they can be just added
            X_shift = ma.grid_to_image(X_shift,scan_positions_xz,X_shift)
            #Y_shift = grid_to_image(Y_shift, positions)
            
            # add misalignments to positions
            positions_new = scan_positions_xz + np.array([Y_shift.transpose().ravel(),X_shift.transpose().ravel()])
            
            # find interpolant
            F = interpolate.LinearNDInterpolator(scan_positions_xz.transpose(), data_xrd, fill_value=np.nan, rescale=True)
            
            # interpolate, this is the output of aligned diffraction data
            xrd_interp = F(positions_new.transpose())
            
            # sum data to plot total xrd map
            xrd_int = np.sum(np.sum(xrd_interp,axis=2),axis=1)
            xrd_int.shape
            xrd_image = np.reshape(xrd_int,data_xrf.transpose().shape).transpose()
            
        #    plt.figure(num=1)
        #    plt.imshow(reference_image)
        
#            plt.figure(num=1,figsize=(15,15))
#            plt.imshow(fluo_aligned)
#            plt.grid()
#            
#            plt.figure(num=3,figsize=(15,15))
#            plt.imshow(xrd_image)
#            plt.grid()        
                         
        # Save the interpolated xrd data ######################################
        if save_data_xrd_interpolate== True and align_data_xrd == True:
            # Save NPZ compressed numpy array
#            save_path = os.path.join(processing_path, "scan_%06d_merlin.npz")%(int(scan_number[ii]))      #The path to save the results
#            np.savez_compressed(save_path,data_xrd = xrd_interp[0:xrd_map_crop,:,:], \
#                                command = command, \
#                                motor_positions = motor_positions, \
#                                scan_positions_x = scan_positions_x, \
#                                scan_positions_z = scan_positions_z)
            
            # Save binary uncompressed array
    #        save_bin_data_nanomax(processing_path, xrd_interp)
            xrd_interp = np.reshape(xrd_interp,[scan_n_z,scan_n_x,220,150])   
            
            # Save mat compressed array for matlab                 
            save_path = os.path.join(processing_path, ("scan_%06d_merlin.mat")%(int(scan_number[ii])))    #The path to save the results
            data_dic = {"data":xrd_interp[:,xrd_map_crop:None,:,:], "command": command, "motor_positions":motor_positions, "scan_positions_x":scan_positions_x[0:xrd_map_crop], "scan_positions_z":scan_positions_z[0:xrd_map_crop]}
            scio.savemat(save_path, data_dic, do_compression=True,oned_as="column") 
            
        if save_data_xrf_interpolate == True:            
            fluo_aligned = np.float32(fluo_aligned/np.max(fluo_aligned))   
            imsave(os.path.join(processing_path, ("scan_%06d_xspress3_aligned.tif")%(int(scan_number[ii]))), fluo_aligned)
            print(('Saved scan %d to %s')%(int(scan_number[ii]),processing_path))
        
    print(('--Processing scan %d done!--')%(int(scan_number[ii])))
    
print("--Overall Processing Done!--")