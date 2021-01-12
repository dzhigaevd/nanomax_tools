#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 13:08:32 2020

@author: dmitry dzhigaev
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/home/dzhigd/Software/') # Important to check!
import nanomax_tools.preprocessing.align_utils as au
import nanomax_tools.preprocessing.read_utils as ru
import nanomax_tools.preprocessing.save_utils as su

# define our clear function 
def clear(): 
  
    # for windows 
    if os.name == 'nt': 
        _ = os.system('cls') 
  
    # for mac and linux(here, os.name is 'posix') 
    else: 
        _ = os.system('clear') 

clear()

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

# INPUTS ######################################################################
#IMPORTANT INFORMATION ABOUT MERLIN DETECTOR IMAGES ORIENTATION
# IN VIEWERS! 
# The images are flipped horizontally in both the LiveView and NanoMax viewer. AND
# vertically flipped in the NanoMax viewer.
# If MERLIN is mounted 180 degrees rotated, the data has to be flipped horizontally to show a view along the beam with imshow
# The indexing of MERLIN starts at bottom left pixel in as designed orientation
#THE END OF IMPORTANT INFORMATION

year        = "2020"                                                           #The year for the experiemnt
beamtimeID  = "2020101408"                                                     #The beamtimeID
sample_name = r"0002_sample_P246_AB"                                           #The name for the sample folder
scan_number = np.linspace(109,128,20)                                          #The scan numbers, can be the list
rocking_motor = "gontheta" # [gonphi,gontheta]

# Options: the most general steps
merlin_180_mount  = True # generalize the detector mount
scanning_xrd_scan = True
normalize         = True
apply_mask        = True
xrd_extract       = True
xrf_extract       = True
save_average_data = True 
save_data         = True

# Numeric parameters
xrf_roi     = [870,970]
xrd_roi     = [0,200,250,400]
#xrd_roi     = [0,515,0,515]
n_scan_points     = 9191
skip_pixel_number = 1
save_type         = ["npz","tif"]
###############################################################################

# Determing the paths"
path_root           = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/%s"%(beamtimeID,sample_name)
mask_path           = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/%s/raw/mask11kev_102020.npz"%(beamtimeID)
processing_pre_path = r"/media/dzhigd/My Passport/DDzhigaev/Data/MAXIV/NanoMax/2020101408/process/";

# DATA IMPORT #################################################################
rocking_angle = np.zeros((len(scan_number),1))
mask_t = ru.read_mask(mask_path, xrd_roi)

if save_average_data == True:
    data_average = np.zeros((n_scan_points,xrd_roi[1]-xrd_roi[0],xrd_roi[3]-xrd_roi[2]))

##### Start here, each scan is processed separately within the rocking curve
for ii in range(0,len(scan_number)):
    print(('\033[1;33;40m #### Processing scan %d')%(int(scan_number[ii])))

    #determing the paths"
    data_meta_path = os.path.join(path_root, "%06d.h5"%(int(scan_number[ii])))
    data_xrd_path  = os.path.join(path_root, "scan_%06d_merlin.hdf5")%(int(scan_number[ii]))
    data_xrf_path  = os.path.join(path_root, "scan_%06d_xspress3.hdf5")%(int(scan_number[ii]))    

    # Load the meta data ###########################################################
    [command, motor_positions, rocking_motor, rocking_angles, \
    scan_position_x, scan_position_y, scan_position_z, \
    incoming_intensity] = ru.read_data_meta(data_meta_path)    

    if xrd_extract == True:
        try:
            data_xrd  = ru.read_data_merlin(data_xrd_path,xrd_roi)
            print("++ XRD data is loaded")
        except:
            raise RuntimeError("-- Cannot read data merlin!")

        if apply_mask == True:
            if ii == 0:
                mask = np.tile(mask_t,(np.shape(data_xrd)[0],1,1))
            data_xrd = np.multiply(data_xrd,mask)
            if merlin_180_mount == True:
                data_xrd = np.flip(data_xrd,2)
            print("++ XRD data was masked")
        if normalize == True:
            for jj in range(0,np.shape(data_xrd)[0]):
                data_xrd[jj,:,:] = data_xrd[jj,:,:]/incoming_intensity[jj]
            print("++ XRD data was normalized")
        
    if xrf_extract == True:
        try:
            data_xrf  = ru.read_data_xspress3(data_xrf_path,xrf_roi)
            print("++ XRF data is loaded")
        except:
            raise RuntimeError("-- Cannot read XRF data!")

        if normalize == True:
            for jj in range(0,np.shape(data_xrf)[0]):
                data_xrf[jj] = data_xrf[jj]/incoming_intensity[jj]
            print("++ XRF data was normalized")
            
    if scanning_xrd_scan == True:
        rocking_motor = "gontheta"
        rocking_angles = motor_positions[rocking_motor]
    	    	    
    ###########################################################################
    # Save the original xrd data separately for each angle #################### 
    if save_data == True:                  
        for value in save_type:
            save_path = create_save_path(processing_pre_path, sample_name, scan_number, value, ii)  
            if value == 'mat':
                # Permute arrays, to have the rocking direction as a last index            
                data_xrd = np.transpose(np.single(data_xrd),(1,2,0))
            if value == 'mat' or value == 'npz':
                save_dictionary = {"data_xrd": data_xrd, \
                            "xrd_roi": xrd_roi, \
                            "motor_positions": motor_positions, \
                            "rocking_motor": rocking_motor, \
                            "rocking_angles": rocking_angles, \
                            "scan_position_x": scan_position_x, \
                            "scan_position_y": scan_position_y, \
                            "scan_position_z": scan_position_z}  
            if value == 'tif':
                data_xrf = np.float32(data_xrf/np.max(data_xrf))
                save_dictionary = data_xrf
            su.save_data(save_path, value, save_dictionary)

    if save_average_data == True:    
        # Save the data average 
        data_average = data_average+data_xrd
            
    print(('\033[1;32;40m #### Processing scan %d done!')%(int(scan_number[ii])))

if save_average_data == True:
    try:
        data_average = data_average/(ii+1)
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        imagery = ax1.imshow(np.log10(np.sum(data_average,0)), cmap='inferno')
        fig.colorbar(imagery, ax=ax1)  
        imagery = ax2.imshow(np.log10(np.sum(data_average,1)), cmap='inferno')
        fig.colorbar(imagery, ax=ax2)
        imagery = ax3.imshow(np.log10(np.sum(data_average,2)), cmap='inferno')
        fig.colorbar(imagery, ax=ax3)
        fig.suptitle("scan_%d_to_%d"%(int(scan_number[0]),int(scan_number[-1])), fontsize=16)
        plt.show()
    except:
        print("-- No average data saved!")

print("\033[1;37;40m #### Overall Processing Done!")