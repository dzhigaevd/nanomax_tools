#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 22:29:50 2020

@author: dzhigd
"""

def save_bin_data_nanomax(save_path,data):
    fid = open(save_path,'w+b')
    fid.write(data)
    fid.close()

