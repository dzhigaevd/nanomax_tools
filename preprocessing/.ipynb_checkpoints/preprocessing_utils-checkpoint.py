#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:40:33 2020

@author: dzhigd
"""
import matplotlib.pyplot as plt
from bcdi.utils import validation as valid
import numpy as np
from bcdi.graph import graph_utils as gu
from bcdi.utils import utilities as util
import sys

sys.path.append('/home/dzhigd/Software/')
import python_tools.numpy_extension.math as nem

def get_overlay_coordinates(n_horizontal,n_vertical,shift):
    """This function calculates linear coordinates of the scan-points in the area of overlap between angular positions after calculation of the shift that has to be read from h5 dataset
    Input: 
    :retun: size of a new map (vertical, horizontal), list of scan point coordinates to be read from h5 dataset for each shift value
    """
    
    shift.insert(0, [0,0])
    shift = np.array(shift)
    n_angles = len(shift)
    
    hv = np.arange(0,n_horizontal)
    vv = np.arange(0,n_vertical)
    h,v = np.meshgrid(hv,vv)
    
    # Add shift to every coordinate of the 2D scan
    shifted_coordinates = []
    for kk in range(len(shift)):
        t = []
        for jj in range(n_vertical):
            for ii in range(n_horizontal):                         
                    t.append(str(v[jj,ii]+shift[kk,0])+','+str(h[jj,ii]+shift[kk,1]))
        shifted_coordinates.append(t)
    
    # Reference 2D map coordinates
    reference_coordinates = shifted_coordinates[0].copy()

    # Find intersection
    for ii in range(len(shifted_coordinates)):
        t_set = frozenset(shifted_coordinates[ii])
        t = [x.split(',') for x in reference_coordinates if x in t_set]
        reference_coordinates = []
        for ii in range(len(t)):
            reference_coordinates.append(t[ii][0]+','+t[ii][1])
            
    common_coordinates = reference_coordinates
    
    if common_coordinates:
        # Find postition of common_coordinates in each 2D map: coordinates
        data_coordinates = []
        for ii in range(n_angles):
            t = []
            for coordinates in common_coordinates:
                t.append(shifted_coordinates[ii].index(coordinates))
            data_coordinates.append(t)

        c_list = []
        for ii in range(len(common_coordinates)):    
            t = common_coordinates[ii].split(',')
            c_list.append([int(float(t[0])),int(float(t[1]))])    
        c_list = np.array(c_list)

        # Size of the common scan range
        n_vertical_new = np.max(c_list[:,0])-np.min(c_list[:,0])+1
        n_horizontal_new = np.max(c_list[:,1])-np.min(c_list[:,1])+1
    else:
        print('There is no overlap between all 2D maps! Stop here.')
        raise ValueError
        
    return data_coordinates, n_vertical_new, n_horizontal_new
    
def roll_2d_frame(frame, horizontal_shift, vertical_shift):
    frame_roll = frame.copy()
    frame_roll = np.roll(frame_roll, -vertical_shift, axis = 0)    # Positive y rolls up
    frame_roll = np.roll(frame_roll, -horizontal_shift, axis = 1)     # Positive x rolls right
    return frame_roll

def normalize(data,axis = None):
    """
    normalize data 
    """
    if axis:
        return (data - np.mean(data, axis)) / (np.std(data, axis))
    else:
        return (data - np.mean(data, axis=0)) / (np.std(data, axis=0))

def calculate_projection(data,axis=None):
    if axis:
        output = np.sum(data,axis)
    else:
        output = np.sum(data,0)
        print("Projected along first dimension")
    return output

def cross_corr(data1, data2):
    """Calculates the cross correlation and lags without normalization.

    The definition of the discrete cross-correlation is in:
    https://www.mathworks.com/help/matlab/ref/xcorr.html

    Args:
    y1, y2: Should have the same length.

    Returns:
    max_corr: Maximum correlation without normalization.
    lag: The lag in terms of the index.
    """
    if len(data1) != len(data2):
        raise ValueError('The lengths of the inputs should be the same.')

    data1_auto_corr = np.dot(data1, data1) / len(data1)
    data2_auto_corr = np.dot(data2, data2) / len(data1)
    corr = np.correlate(data1, data2, mode='same')
    # The unbiased sample size is N - lag.
    unbiased_sample_size = np.correlate(np.ones(len(data1)), np.ones(len(data1)), mode='same')
    corr = corr / unbiased_sample_size / np.sqrt(data1_auto_corr * data2_auto_corr)
    shift = len(data1) // 2

    max_corr = np.max(corr)
    argmax_corr = np.argmax(corr)#
    return max_corr, argmax_corr - shift


def COM_voxels_sparse(data, Qx, Qz, Qy ):
    """
    calculates center of mass 
    """
    # center of mass calculation in reciprocal space with the meshgrids
    COM_qx = np.sum(data.flatten()* Qx.flatten())/np.sum(data.flatten())
    COM_qz = np.sum(data.flatten()* Qz.flatten())/np.sum(data.flatten())
    COM_qy = np.sum(data.flatten()* Qy.flatten())/np.sum(data.flatten())

    return COM_qx, COM_qz, COM_qy

# def COM_voxels_reciproc(data, Qx, Qz, Qy ):
#     """
#     calculates center of mass 
#     """
#     # center of mass calculation in reciprocal space with the meshgrids
#     COM_qx = np.sum(data* Qx)/np.sum(data)
#     COM_qz = np.sum(data* Qz)/np.sum(data)
#     COM_qy = np.sum(data* Qy)/np.sum(data)

#     #print( 'coordinates in reciprocal space:')
#     #print( COM_qx, COM_qz, COM_qy)
#     return COM_qx, COM_qz, COM_qy


def get_q_coordinates(setup,number_of_scans,roi_xrd=None):
    """
    Calculate q-space coordinates for each pixel in the rocking curve dataset
    
    """
    # Constants. They are needed for correct labeling of axes
    H = 4.1357e-15;                                  # Plank's constant
    C = 2.99792458e8;                                # Speed of light in vacuum
    
    number_of_pixels_horizontal = setup.detector.nb_pixel_x
    number_of_pixels_vertical = setup.detector.nb_pixel_y
    # 3rd dimension: number_of_scans
    
    # Correction angles
    detector_delta_correction = 0
    detector_gamma_correction = 0

    # NanoMax convention:
    # gamma - horizontal detector
    # delta - vertical detector
    # gonphi - rotation about vertical axis
    # gontheta - rotation about horizontal axis
    radius = setup.distance
    photon_energy   = setup.energy

    gonphi_increment = setup.tilt_angle # [deg] - can be a range of angles
    gontheta = setup.sample_outofplane
    gonphi   = setup.sample_inplane
    
    delta    = setup.outofplane_angle+detector_delta_correction # [deg] these angles are corrected with the sign respecting the rotation rules
    gamma    = setup.inplane_angle+detector_gamma_correction # [deg] 
    
    if setup.detector.pixelsize_x != setup.detector.pixelsize_y:
        print('Pixels of the detector are not square! Not implemented, stop here.')
        raise ValueError
    else:
        detector_pitch  = setup.detector.pixelsize_x # [m]

    direct_beam     = np.round([251.7859,  250.4288])

    wavelength = H*C/photon_energy

    k = 2*np.pi/wavelength # wave vector

    dq = k*2*np.arctan(detector_pitch/(2*radius)) # q-space pitch at the detector plane

    hd,vd = np.meshgrid(np.arange(-number_of_pixels_horizontal/2,number_of_pixels_horizontal/2),np.arange(number_of_pixels_vertical/2,-number_of_pixels_vertical/2,-1))

    hd = (hd+number_of_pixels_horizontal/2-direct_beam[1])*detector_pitch;
    vd = (vd+number_of_pixels_vertical/2-direct_beam[0])*detector_pitch;
    zd = np.ones(vd.shape)*radius;

    # Add functionality to image display
    # gu.imagesc(hd[0,:],vd[:,0],(np.sum(data,0)));

    # Data reduction    
    # data = data[roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3],:]
    hd = hd[roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3]]
    vd = vd[roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3]]
    zd = zd[roi_xrd[0]:roi_xrd[1],roi_xrd[2]:roi_xrd[3]]
    
    number_of_pixels_horizontal = roi_xrd[3]-roi_xrd[2]
    number_of_pixels_vertical = roi_xrd[1]-roi_xrd[0]
    
    d = np.array([hd.flatten(),vd.flatten(),zd.flatten()])

    r = np.sqrt(np.sum(d**2,0))

    hq = k*(d[0,:]/r)
    vq = k*(d[1,:]/r)
    zq = k*(1-d[2,:]/r)

    q = [hq,vq,zq]

    # Sample orientation matrix. Bounds the sample crystal with the laboratory frame
    # Angles alpha beta gamma were manually adjusted so that known peaks 
    # are exactly in their places

    # X is horizontal, perp to the beam, Y is vertical

    Rh = np.array([[1, 0,              0], # detector rotation around horizontal axis 
                  [0, nem.cosd(delta), -nem.sind(delta)],
                  [0, nem.sind(delta),  nem.cosd(delta)]])

    Rv = np.array([[nem.cosd(gamma),  0,  nem.sind(gamma)], # detector rotation around vertical axis 
                  [0,               1,  0],
                  [-nem.sind(gamma), 0,  nem.cosd(gamma)]])

    Rz = np.array([[nem.cosd(0), -nem.sind(0), 0], # detector rotation around beam axis 
                  [nem.sind(0),  nem.cosd(0), 0],
                  [0,           0,          1]])

    U = Rh@Rv@Rz
    qR = (U@q) # correct so far in real space

    # Initial coordinate of ki
    ki = np.array([0,0,k])
    kf = U@ki
    Q = kf-ki

    # Lab coordinate system: accosiated with the ki
    QLab = [qR[0,:]+Q[0], qR[1,:]+Q[1], qR[2,:]+Q[2]]
    
    # Small corrections to misalignment of the sample
    # Here the rocking curve should be introduced
    # alpha
    # beta

    # Gonphi correction
    sample_alpha_correction = 0 # Qx+Qz 
    sample_beta_correction  = 0 # Qz should be positive
    sample_gamma_correction = -gonphi_increment*number_of_scans/2 # Qz
    dphi = 0.1

    q_values = []

    for ii in range(number_of_scans):
        # Rotations to bring the q vector into sample coordinate system
        Rsh = np.array([[1,         0,                                           0], # detector rotation around horizintal axis 
                        [0,         nem.cosd(gontheta+sample_alpha_correction), -nem.sind(gontheta+sample_alpha_correction)],
                        [0,         nem.sind(gontheta+sample_alpha_correction),  nem.cosd(gontheta+sample_alpha_correction)]]) 

        Rsv = np.array([[nem.cosd(-gamma/2 + gonphi_increment*(ii-1)+sample_beta_correction),  0,  nem.sind(-gamma/2 + gonphi_increment*(ii-1)+sample_beta_correction)], # detector rotation around vertical axis 
                        [0,                                                        1,  0],
                        [-nem.sind(-gamma/2 + gonphi_increment*(ii-1)+sample_beta_correction), 0,  nem.cosd(-gamma/2 + gonphi_increment*(ii-1)+sample_beta_correction)]])

        Rsz = np.array([[nem.cosd(sample_gamma_correction), -nem.sind(sample_gamma_correction), 0], 
                        [nem.sind(sample_gamma_correction),  nem.cosd(sample_gamma_correction), 0],
                        [0,                                  0,                                 1]])

        Rs = Rsh@Rsv@Rsz

        # Sample coordinate system: accosiated with the ki
        q_values.append(Rs@QLab)

    q_values = np.array(q_values)
    
    return q_values


# adapted from bcdi package by carnisj #
def grid_bcdi_labframe_backup(
    data,
    mask,
    detector,
    setup,
    align_q=True,
    reference_axis=(1, 0, 0),
    debugging=False,
    **kwargs,
):
    """
    Interpolate BCDI reciprocal space data using a linearized transformation matrix.
    The resulting (qx, qy, qz) are in the laboratory frame (qx downstrean,
    qz vertical up, qy outboard).
    :param data: the 3D data, already binned in the detector frame
    :param mask: the corresponding 3D mask
    :param detector: an instance of the class Detector
    :param setup: instance of the Class experiment_utils.Setup()
    :param align_q: boolean, if True the data will be rotated such that q is along
     reference_axis, and q values will be calculated in the pseudo crystal frame.
    :param reference_axis: 3D vector along which q will be aligned, expressed in an
     orthonormal frame x y z
    :param debugging: set to True to see plots
    :param kwargs:
     - 'cmap': str, name of the colormap
     - 'fill_value': tuple of two real numbers, fill values to use for pixels outside
       of the interpolation range. The first value is for the data, the second for the
       mask. Default is (0, 0)
     - 'logger': an optional logger
    :return:
     - the data interpolated in the laboratory frame
     - the mask interpolated in the laboratory frame
     - a tuple of three 1D vectors of q values (qx, qz, qy)
     - a numpy array of shape (3, 3): transformation matrix from the detector
       frame to the laboratory/crystal frame
    """
    # logger = kwargs.get("logger", module_logger)
    # valid.valid_ndarray(arrays=(data, mask), ndim=3)
    # check and load kwargs
    # valid.valid_kwargs(
    #     kwargs=kwargs,
    #     allowed_kwargs={"cmap", "fill_value", "logger", "reference_axis"},
    #     name="kwargs",
    # )
    cmap = kwargs.get("cmap", "turbo")
    fill_value = kwargs.get("fill_value", (0, 0))
    # valid.valid_container(
    #     fill_value,
    #     container_types=(tuple, list, np.ndarray),
    #     length=2,
    #     item_types=Real,
    #     name="fill_value",
    # )

    # check some parameters
    if setup.rocking_angle == "energy":
        raise NotImplementedError(
            "Geometric transformation not yet implemented for energy scans"
        )
    # valid.valid_item(align_q, allowed_types=bool, name="align_q")
    # valid.valid_container(
    #     reference_axis,
    #     container_types=(tuple, list, np.ndarray),
    #     length=3,
    #     item_types=Real,
    #     name="reference_axis",
    # )
    reference_axis = np.array(reference_axis)

    # grid the data
    # logger.info(
    #     "Gridding the data using the linearized matrix, "
    #     "the result will be in the laboratory frame"
    # )
    string = "linmat_reciprocal_space_"
    (interp_data, interp_mask), q_values, transfer_matrix = setup.ortho_reciprocal(
        arrays=(data, mask),
        verbose=True,
        debugging=debugging,
        fill_value=fill_value,
        align_q=align_q,
        reference_axis=reference_axis,
        scale=("log", "linear"),
        title=("data", "mask"),
    )
    qx, qz, qy = q_values

    # check for Nan
    interp_mask[np.isnan(interp_data)] = 1
    interp_data[np.isnan(interp_data)] = 0
    interp_mask[np.isnan(interp_mask)] = 1
    # set the mask as an array of integers, 0 or 1
    interp_mask[np.nonzero(interp_mask)] = 1
    interp_mask = interp_mask.astype(int)

    # apply the mask to the data
    interp_data[np.nonzero(interp_mask)] = 0

    # save plots of the gridded data
    final_binning = (
        detector.preprocessing_binning[0] * detector.binning[0],
        detector.preprocessing_binning[1] * detector.binning[1],
        detector.preprocessing_binning[2] * detector.binning[2],
    )

    numz, numy, numx = interp_data.shape
    plot_comment = (
        f"_{numz}_{numy}_{numx}_"
        f"{final_binning[0]}_{final_binning[1]}_{final_binning[2]}.png"
    )

    max_z = interp_data.sum(axis=0).max()
    # fig, _, _ = gu.contour_slices(
    #     interp_data,
    #    (qx, qz, qy),
    #    sum_frames=True,
    #    title="Regridded data",
    #    levels=np.linspace(0, np.ceil(np.log10(max_z)), 150, endpoint=True),
    #    plot_colorbar=True,
    #    scale="log",
    #    is_orthogonal=True,
    #    reciprocal_space=True,
    #    cmap=cmap,
    # )
    # fig.savefig(detector.savedir + string + "sum" + plot_comment)
    # plt.close(fig)

    # fig, _, _ = gu.contour_slices(
    #     interp_data,
    #     (qx, qz, qy),
    #     sum_frames=False,
    #     title="Regridded data",
    #     levels=np.linspace(0, np.ceil(np.log10(interp_data.max())), 150, endpoint=True),
    #     plot_colorbar=True,
    #     scale="log",
    #     is_orthogonal=True,
    #     reciprocal_space=True,
    #     cmap=cmap,
    # )
#     fig.savefig(detector.savedir + string + "central" + plot_comment)
#     plt.close(fig)

#     fig, _, _ = gu.multislices_plot(
#         interp_data,
#         sum_frames=True,
#         scale="log",
#         plot_colorbar=True,
#         vmin=0,
#         title="Regridded data",
#         is_orthogonal=True,
#         reciprocal_space=True,
#         cmap=cmap,
#     )
#     fig.savefig(detector.savedir + string + "sum_pix" + plot_comment)
#     plt.close(fig)

#     fig, _, _ = gu.multislices_plot(
#         interp_data,
#         sum_frames=False,
#         scale="log",
#         plot_colorbar=True,
#         vmin=0,
#         title="Regridded data",
#         is_orthogonal=True,
#         reciprocal_space=True,
#         cmap=cmap,
#     )
#     fig.savefig(detector.savedir + string + "central_pix" + plot_comment)
#     plt.close(fig)
    if debugging:
        gu.multislices_plot(
            interp_mask,
            sum_frames=False,
            scale="linear",
            plot_colorbar=True,
            vmin=0,
            title="Regridded mask",
            is_orthogonal=True,
            reciprocal_space=True,
            cmap=cmap,
        )

    return interp_data, interp_mask, q_values, transfer_matrix

def grid_bcdi_labframe(
    data,
    mask,
    detector,
    setup,
    align_q=True,
    reference_axis=(1, 0, 0),
    debugging=False,
    **kwargs,
):
    """
    Interpolate BCDI reciprocal space data using a linearized transformation matrix.
    The resulting (qx, qy, qz) are in the laboratory frame (qx downstrean,
    qz vertical up, qy outboard).
    :param data: the 3D data, already binned in the detector frame
    :param mask: the corresponding 3D mask
    :param detector: an instance of the class Detector
    :param setup: instance of the Class experiment_utils.Setup()
    :param align_q: boolean, if True the data will be rotated such that q is along
     reference_axis, and q values will be calculated in the pseudo crystal frame.
    :param reference_axis: 3D vector along which q will be aligned, expressed in an
     orthonormal frame x y z
    :param debugging: set to True to see plots
    :param kwargs:
     - 'cmap': str, name of the colormap
     - 'fill_value': tuple of two real numbers, fill values to use for pixels outside
       of the interpolation range. The first value is for the data, the second for the
       mask. Default is (0, 0)
     - 'logger': an optional logger
    :return:
     - the data interpolated in the laboratory frame
     - the mask interpolated in the laboratory frame
     - a tuple of three 1D vectors of q values (qx, qz, qy)
     - a numpy array of shape (3, 3): transformation matrix from the detector
       frame to the laboratory/crystal frame
    """
    # logger = kwargs.get("logger", module_logger)
    # valid.valid_ndarray(arrays=(data, mask), ndim=3)
    # check and load kwargs
    # valid.valid_kwargs(
    #     kwargs=kwargs,
    #     allowed_kwargs={"cmap", "fill_value", "logger", "reference_axis"},
    #     name="kwargs",
    # )
    cmap = kwargs.get("cmap", "turbo")
    fill_value = kwargs.get("fill_value", (0, 0))
    # valid.valid_container(
    #     fill_value,
    #     container_types=(tuple, list, np.ndarray),
    #     length=2,
    #     item_types=Real,
    #     name="fill_value",
    # )

    # check some parameters
    if setup.rocking_angle == "energy":
        raise NotImplementedError(
            "Geometric transformation not yet implemented for energy scans"
        )
    # valid.valid_item(align_q, allowed_types=bool, name="align_q")
    # valid.valid_container(
    #     reference_axis,
    #     container_types=(tuple, list, np.ndarray),
    #     length=3,
    #     item_types=Real,
    #     name="reference_axis",
    # )
    reference_axis = np.array(reference_axis)

    # grid the data
    # logger.info(
    #     "Gridding the data using the linearized matrix, "
    #     "the result will be in the laboratory frame"
    # )
    string = "linmat_reciprocal_space_"
    (interp_data, interp_mask), q_values, transfer_matrix = setup.ortho_reciprocal(
        arrays=(data, mask),
        verbose=True,
        debugging=debugging,
        fill_value=fill_value,
        align_q=align_q,
        reference_axis=reference_axis,
        scale=("log", "linear"),
        title=("data", "mask"),
    )
    qx, qz, qy = q_values

    # check for Nan
    interp_mask[np.isnan(interp_data)] = 1
    interp_data[np.isnan(interp_data)] = 0
    interp_mask[np.isnan(interp_mask)] = 1
    # set the mask as an array of integers, 0 or 1
    interp_mask[np.nonzero(interp_mask)] = 1
    interp_mask = interp_mask.astype(int)

    # apply the mask to the data
    interp_data[np.nonzero(interp_mask)] = 0

    # save plots of the gridded data
    final_binning = (
        detector.preprocessing_binning[0] * detector.binning[0],
        detector.preprocessing_binning[1] * detector.binning[1],
        detector.preprocessing_binning[2] * detector.binning[2],
    )

    numz, numy, numx = interp_data.shape
    plot_comment = (
        f"_{numz}_{numy}_{numx}_"
        f"{final_binning[0]}_{final_binning[1]}_{final_binning[2]}.png"
    )

    max_z = interp_data.sum(axis=0).max()
    # fig, _, _ = gu.contour_slices(
    #     interp_data,
    #    (qx, qz, qy),
    #    sum_frames=True,
    #    title="Regridded data",
    #    levels=np.linspace(0, np.ceil(np.log10(max_z)), 150, endpoint=True),
    #    plot_colorbar=True,
    #    scale="log",
    #    is_orthogonal=True,
    #    reciprocal_space=True,
    #    cmap=cmap,
    # )
    # fig.savefig(detector.savedir + string + "sum" + plot_comment)
    # plt.close(fig)

    # fig, _, _ = gu.contour_slices(
    #     interp_data,
    #     (qx, qz, qy),
    #     sum_frames=False,
    #     title="Regridded data",
    #     levels=np.linspace(0, np.ceil(np.log10(interp_data.max())), 150, endpoint=True),
    #     plot_colorbar=True,
    #     scale="log",
    #     is_orthogonal=True,
    #     reciprocal_space=True,
    #     cmap=cmap,
    # )
#     fig.savefig(detector.savedir + string + "central" + plot_comment)
#     plt.close(fig)

#     fig, _, _ = gu.multislices_plot(
#         interp_data,
#         sum_frames=True,
#         scale="log",
#         plot_colorbar=True,
#         vmin=0,
#         title="Regridded data",
#         is_orthogonal=True,
#         reciprocal_space=True,
#         cmap=cmap,
#     )
#     fig.savefig(detector.savedir + string + "sum_pix" + plot_comment)
#     plt.close(fig)

#     fig, _, _ = gu.multislices_plot(
#         interp_data,
#         sum_frames=False,
#         scale="log",
#         plot_colorbar=True,
#         vmin=0,
#         title="Regridded data",
#         is_orthogonal=True,
#         reciprocal_space=True,
#         cmap=cmap,
#     )
#     fig.savefig(detector.savedir + string + "central_pix" + plot_comment)
#     plt.close(fig)
    if debugging:
        gu.multislices_plot(
            interp_mask,
            sum_frames=False,
            scale="linear",
            plot_colorbar=True,
            vmin=0,
            title="Regridded mask",
            is_orthogonal=True,
            reciprocal_space=True,
            cmap=cmap,
        )

    return q_values
