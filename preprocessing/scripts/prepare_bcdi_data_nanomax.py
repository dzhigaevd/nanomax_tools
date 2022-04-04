#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dmitry dzhigaev
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
from numpy.core.records import fromarrays

# INPUTS ######################################################################
#IMPORTANT INFORMATION ABOUT MERLIN DETECTOR IMAGES ORIENTATION
# IN VIEWERS! 
# The images are flipped horizontally in both the LiveView and NanoMax viewer. AND
# vertically flipped in the NanoMax viewer.
# If MERLIN is mounted 180 degrees rotated, the data has to be flipped horizontally to show a view along the beam with imshow
# The indexing of MERLIN starts at bottom left pixel in as designed orientation
#THE END OF IMPORTANT INFORMATION

year        = "2020"                                                           #The year for the experiemnt
beamtimeID  ="2020062408"                                                      #The beamtimeID
sample_name = r"sample0609"                                                    #The name for the sample folder
scan_number = [442] # The scan numbers, can be the list
rocking_motor = "gonphi" # [gonphi,gontheta]

# Options: the most general steps
merlin_180_mount = True # generalize the detector mount
spatial_scan     = False
normalize        = False
apply_mask       = True
xrd_analysis     = True
xrf_analysis     = False
q_space          = False
save_data_xrd_m  = True # save for MATLAB
save_data_xrd_p  = False # save for PYTHON
save_data_xrf    = False

xrf_roi     = [870,970]
xrd_roi     = [0,515,0,515]

skip_pixel_number = 1
###############################################################################

#determing the paths"
# Path to the raw data
path_root       = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/%s"%(beamtimeID,sample_name)
# Mask of bad pixels
mask_path       = r"/media/dzhigd/My Passport1/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/mask11kev_102020.npz"%(beamtimeID)

# DATA IMPORT #################################################################
rocking_angle = np.zeros((len(scan_number),1))
mask_t = ru.read_mask(mask_path, xrd_roi)

##### Start here

for ii in range(0,len(scan_number)):
    print(('--Processing scan %d--')%(int(scan_number[ii])))
    
    #determing the paths"
    data_meta_path = os.path.join(path_root, "%06d.h5"%(int(scan_number[ii])))
    data_xrd_path  = os.path.join(path_root, "scan_%06d_merlin.hdf5")%(int(scan_number[ii]))
    data_xrf_path  = os.path.join(path_root, "scan_%06d_xspress3.hdf5")%(int(scan_number[ii]))    
    
    # Load the data ###########################################################
    [command, motor_positions, rocking_motor, rocking_angles, scan_position_x, scan_position_z, incoming_intensity] = ru.read_data_meta(data_meta_path)

    if xrd_analysis == True:
        data_xrd  = ru.read_data_merlin(data_xrd_path,xrd_roi)
        if apply_mask == True:
            if ii == 0:
                mask = np.tile(mask_t,(np.shape(data_xrd)[0],1,1))
            data_xrd = np.multiply(data_xrd,mask)
            if merlin_180_mount == True:
                data_xrd = np.flip(data_xrd,2)
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
    ###########################################################################

    # Show the data ###########################################################
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    imagery = ax1.imshow(np.log10(np.sum(data_xrd,0)), cmap='inferno')
    fig.colorbar(imagery, ax=ax1)  
    imagery = ax2.imshow(np.log10(np.sum(data_xrd,1)), cmap='inferno')
    fig.colorbar(imagery, ax=ax2)
    imagery = ax3.imshow(np.log10(np.sum(data_xrd,2)), cmap='inferno')
    fig.colorbar(imagery, ax=ax3)
    fig.suptitle("scan_%06d_merlin"%int(scan_number[ii]), fontsize=16)
    # Save the original xrd data ##############################################
    if save_data_xrd_m == True:
        processing_path = r"/home/dzhigd/work/projects/CsPbBr3_NC_BCDI_NanoMAX/data/%s_%d"%(sample_name,scan_number[0])
        try:
            os.mkdir(processing_path)
            print("--Processing folder created--")
        except:
            print("--Processing folder already exists--")
        
        # Save mat compressed array for matlab                 
        save_path = os.path.join(processing_path, ("scan_%06d_merlin.mat")%(int(scan_number[ii])))    #The path to save the results
        # Permute arrays, to have the rocking direction as a last index
        data_xrd = np.transpose(np.single(data_xrd),(1,2,0))
        
        data_dic = {"data": data_xrd,\
                    "mask": mask,\
                    "command": command, "motor_positions":motor_positions,\
                    "rocking_motor":rocking_motor, "rocking_angles":rocking_angles}
#        scan = fromarrays([data_xrd, mask, command, motor_positions, rocking_motor, rocking_angles], names=['data_xrd', 'mask', 'command', 'motor_positions', 'rocking_motor', 'rocking_angles'])
        scio.savemat(save_path, {'scan': data_dic}, do_compression=True,oned_as="column") 

        print(('Saved scan %d to %s')%(int(scan_number[ii]),save_path))
#        save_path = os.path.join(processing_path, "scan_%06d_merlin_%d_%d_%d.npz")%(int(scan_number[ii]),(data_xrd.shape[0]),(data_xrd.shape[1]),(data_xrd.shape[2]))      #The path to save the results
#        save_bin_data_nanomax(save_path, data_xrd)
#        np.savez_compressed(save_path,data_xrd)
        
    if save_data_xrf == True:
        data_xrf = np.float32(data_xrf/np.max(data_xrf))
        imsave(os.path.join(processing_path, ("scan_%06d_xspress3.tif")%(int(scan_number[ii]))), data_xrf)
        print(('Saved scan %d to %s')%(int(scan_number[ii]),processing_path))

    print(('--Processing scan %d done!--')%(int(scan_number[ii])))

print("--Overall Processing Done!--")
