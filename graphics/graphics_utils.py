#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 23:43:44 2020

@author: dzhigd
"""
from scipy import ndimage, interpolate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import plotly.express as px
import plotly.graph_objects as go
from mpl_toolkits.axes_grid1 import make_axes_locatable

class IndexTracker:
    def __init__(self, ax, data, axis):
        self.ax = ax
        self.scroll_axis = axis
        ax.set_title('use scroll wheel to navigate images')
        self.data = data                
        
        # Start slice to show
        self.ind = 0
        
        if self.scroll_axis == 0:
            self.slices, rows, cols = data.shape
            self.im = ax.imshow(self.data[self.ind, :, :])
        elif self.scroll_axis == 1:
            rows, self.slices, cols = data.shape
            self.im = ax.imshow(self.data[:, self.ind, :])
        elif self.scroll_axis == 2:
            rows, cols, self.slices = data.shape
            self.im = ax.imshow(self.data[:, :, self.ind])
                
        # plt.colorbar(self.im, self.ax)
        self.update()

    def on_scroll(self, event):
        # print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind - 1) % self.slices
        else:
            self.ind = (self.ind + 1) % self.slices
        self.update()

    def update(self):
        if self.scroll_axis == 0:
            self.im.set_data(self.data[self.ind, :, :])
        elif self.scroll_axis == 1:
            self.im.set_data(self.data[:, self.ind, :])
        elif self.scroll_axis == 2:
            self.im.set_data(self.data[:, :, self.ind])
            
        self.ax.set_title('Slice %d along axis %d ' % (self.ind, self.scroll_axis))
        
        self.im.axes.figure.canvas.draw()        

def scroll_data(data, axis=None, colormap=None):
    fig, ax = plt.subplots(1, 1)    
    
    if axis:
        tracker = IndexTracker(ax, data, axis)
    else:
        tracker = IndexTracker(ax, data, 0)        
    
    fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
    if colormap:
        plt.set_cmap(colormap)
    else:
        plt.set_cmap('turbo')
    
    plt.show()
    
    return tracker
    
def show_q_space_projections(qx,qy,qz,data):    
#     Constants. They are needed for correct labeling of axes
    h                       = 4.1357e-15;                                  # Plank's constant
    c                       = 2.99792458e8;                                # Speed of light in vacuum
    energy = 11000
    pitch = 55e-6
    radius = 0.4
    
    wavelength = h*c/energy

    dq = (2*np.pi*pitch/(radius*wavelength))
    
    qxv = np.arange(np.min(qx),np.max(qx),dq)
    qyv = np.arange(np.min(qy),np.max(qy),dq)
    qzv = np.arange(np.min(qz),np.max(qz),dq)
    
    Qx,Qy,Qz = np.meshgrid(qxv, qyv, qzv)
    
    # find interpolant
    q = np.array((qz,qx,qy)).transpose()
    F = interpolate.LinearNDInterpolator(q, np.concatenate(np.concatenate(data)), fill_value=np.nan, rescale=True)
    
    # interpolate, this is the output of aligned diffraction data
    qspace_interpolated = np.nan_to_num(F(Qz.transpose(),Qx.transpose(),Qy.transpose()))
    
#    qspace_interpolated[qspace_interpolated == np.nan] = 0
    
    # sum data to plot total xrd map    
    plt.figure(num=1)
    plt.imshow(np.sum(qspace_interpolated,axis=0))
    
    plt.figure(num=2)
    plt.imshow(np.log10(np.sum(qspace_interpolated,axis=1)))
    
    plt.figure(num=3)
    plt.imshow(np.sum(qspace_interpolated,axis=2))
    
def imagesc(*args, cmap='turbo', xlabel=None, ylabel=None, title=None, levels = 300):
    if len(args) == 1:
        plt.figure()
        plt.imshow(args[0],cmap=cmap) 
        plt.colorbar()
        plt.show()
        
        if title!=None:
            plt.title(title)        

        if xlabel!=None:
            plt.xlabel(xlabel)

        if ylabel!=None:
            plt.ylabel(ylabel)    
    
    elif len(args) == 3:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        
        hv = args[0]
        vv = args[1]
        
        display = ax.imshow(args[2],extent=(np.min(hv), np.max(hv), np.min(vv), np.max(vv)))
        
        if cmap!=None:
            display.set_cmap(cmap)
        else:
            display.set_cmap('turbo') 

        if xlabel!=None:
            ax.set_xlabel(xlabel)

        if ylabel!=None:
            ax.set_ylabel(ylabel)    

        cbar = fig.colorbar(display)

        if title!=None:
            ax.title(title)        

        ax.set_aspect("equal")

def imagesc_central_slices(*args, cmap='turbo', xlabel=None, ylabel=None, title=None, levels = 300):
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3)
    image = args[0]
    ax1.imshow(image[int(image.shape[0]/2),:,:])
    ax2.imshow(image[:,int(image.shape[1]/2),:])
    ax3.imshow(image[:,:,int(image.shape[2]/2)])
    
def imagesc_ortho_projections(*args, cmap=None, xlabel=None, ylabel=None, title=None, levels = 300):
    def set_cmap_ortho(ax1,ax2,ax3,d1,d2,d3,cmap):
        if cmap!=None:
            d1 = ax1.set_cmap(cmap)
            d2 = ax2.set_cmap(cmap)
            d3 = ax3.set_cmap(cmap)
        else:
            d1.set_cmap('turbo')
            d2.set_cmap('turbo')
            d3.set_cmap('turbo')
            
    if len(args) == 1:
        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3)
        image = args[0]
        d1 = ax1.imshow(np.sum(image,0))
        d2 = ax2.imshow(np.sum(image,1))
        d3 = ax3.imshow(np.sum(image,2))                
        set_cmap_ortho(ax1,ax2,ax3,d1,d2,d3,cmap)
        if title!=None:
            ax.title(title)    

    elif len(args) == 4:
        fig, (ax1,ax2,ax3) = plt.subplots(ncols=3)
        # fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
        bv = args[0]*1e-10
        vv = args[1]*1e-10 
        hv = args[2]*1e-10 

        d1 = ax1.imshow(np.sum(args[3],0),extent=(np.min(hv), np.max(hv), np.min(bv), np.max(bv)))
        d2 = ax2.imshow(np.sum(args[3],1),extent=(np.min(hv), np.max(hv), np.min(vv), np.max(vv)))
        d3 = ax3.imshow(np.sum(args[3],2),extent=(np.min(vv), np.max(vv), np.min(bv), np.max(bv)))                
                
#         if xlabel!=None:
#             ax.set_xlabel(xlabel)

#         if ylabel!=None:
#             ax.set_ylabel(ylabel)
    set_cmap_ortho(ax1,ax2,ax3,d1,d2,d3,cmap)
    ax1.set_aspect('auto')
    ax2.set_aspect('auto')
    ax3.set_aspect('auto')
    plt.colorbar(d1,ax=ax1)
    plt.colorbar(d2,ax=ax2)
    plt.colorbar(d3,ax=ax3)
    plt.show()
    # cb = plt.colorbar(d1, label = 'Integrated intensity')    

#         if title!=None:
#             ax.title(title)        

        # ax.set_aspect("equal")        
def imagesc_plotly(*args, cmap='turbo', xlabel=None, ylabel=None, title=None):
    # fig = px.imshow(args[0], color_continuous_scale=cmap, origin='lower')
    fig = go.Figure(data=go.Heatmap(
          x = args[0],
          y = args[1],
          z = args[2],
          type = 'heatmap',
          colorscale = cmap))
    fig.show()
    
def scatter3_plotly(x,y,z,values):
    fig = go.Figure(data=[go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=12,
            color=values,                # set color to an array/list of desired values
            colorscale='Turbo',   # choose a colorscale
            opacity=0.8
        )
    )])

    # tight layout
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.show()