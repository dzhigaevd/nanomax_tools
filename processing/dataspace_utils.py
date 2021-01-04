#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:21:36 2020

@author: dzhigd
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.measurements import center_of_mass
from skimage.viewer.canvastools import RectangleTool
from skimage.viewer import ImageViewer

# Crop the data based on selected rois
def get_bragg_peak_com(data,n_peaks):
    fig,ax = plt.subplots()
#    for ii in range(n_peaks):        
    image = ax.imshow(np.log10(np.sum(data,0)),cmap='inferno')
    
    props = {'facecolor': '#000070',
         'edgecolor': 'white',
         'alpha': 0.3}
    rect_tool = RectangleTool(image, rect_props=props)
    
    plt.show()
    print("Final selection:")
    rect_tool.callback_on_enter(rect_tool.extents)
    print("The roi is selected %s"%str(rect_tool.extents))
    plt.close()
#    return 
    
def bin2npz(data_path, dimensions):
    
    data = np.fromfile(data_path,dtype='single')
    data = np.reshape(data, (dimensions[2],dimensions[0],dimensions[1]))  
    # This is just specific case for
#    fid = fopen([save_path,sprintf('data_peak_%d_%d_%d_%d.bin',ii,size(data))],'wb');
#    fwrite(fid,data,'single');
#    fclose(fid);
    
    print(data.shape)
    
    plt.figure(1)
    plt.imshow(np.abs(data[np.int(np.round(dimensions[2]/2)),:,:]),cmap = 'inferno')
    plt.show()
    
    plt.figure(2)
    plt.imshow(np.log10(data[np.int(np.round(dimensions[2]/2)),:,:]),cmap = 'inferno')
    plt.show()
    
    np.savez_compressed(data_path[0:-4]+'.npz', data = data)