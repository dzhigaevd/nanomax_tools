# -*- coding: utf-8 -*-


"""
Scan class that defines the experimental geometry.

You can think of it as the public interface for the Beamline and Diffractometer child
classes. A script would call a method from Setup, which would then retrieve the required
beamline-dependent information from the child classes.
"""
import datetime
import logging
import multiprocessing as mp
import sys
import time
import numpy as np


class Scan:
    """
    Class for defining the scan at NanoMAX.
    
    """

    def __init__(
        self,
        setup,
        scan_command=None,
        scan_orientation=None,
        scan_numbers=None
        **kwargs,
    ):
        
        # kwargs for loading and preprocessing data
        self.dirbeam_detector_angles = kwargs.get("dirbeam_detector_angles")
        self.direct_beam = kwargs.get("direct_beam")
        self.filtered_data = kwargs.get("filtered_data", False)  # boolean
        self.custom_scan = kwargs.get("custom_scan", False)  # boolean
        self.custom_images = kwargs.get("custom_images")  # list or tuple
        self.custom_monitor = kwargs.get("custom_monitor")  # list or tuple
        self.custom_motors = kwargs.get("custom_motors")  # dictionnary
        self.actuators = kwargs.get("actuators", {})  # list or tuple
        
        

   