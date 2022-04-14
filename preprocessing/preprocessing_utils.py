#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:40:33 2020

@author: dzhigd
"""

import numpy as np

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