# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:28:20 2020

@author: Admin
"""

#!/usr/bin/env python
import h5py 
import numpy as np 
import math
from matplotlib import pyplot as plt
import nanomax_tools.processing.dataspace_utils as du

data_path    = "/home/dzhigd/work/projects/CsPbBr3_NC_BCDI_NanoMAX/data/sample0609_442/data_peak_2_67_67_75.bin"

du.bin2npz(data_path,[67,67,75])
