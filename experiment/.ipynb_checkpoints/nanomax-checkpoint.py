import datetime
import logging
import multiprocessing as mp
import sys
import time
from collections.abc import Sequence
from numbers import Integral, Real
from typing import Optional, Tuple

import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata

from bcdi.experiment.beamline import create_beamline
from bcdi.experiment.detector import Detector, create_detector
from bcdi.graph import graph_utils as gu
from bcdi.utils import utilities as util
from bcdi.utils import validation as valid

module_logger = logging.getLogger(__name__)


class Setup:
    """

    """

    labframe_to_xrayutil = {
        "x+": "y+",
        "x-": "y-",
        "y+": "z+",
        "y-": "z-",
        "z+": "x+",
        "z-": "x-",
    }  # conversion table from the laboratory frame (CXI convention)
    # (z downstream, y vertical up, x outboard) to the frame of xrayutilities
    # (x downstream, y outboard, z vertical up)

    def __init__(
        self,
        detector_name="Merlin",
        beam_direction=(1, 0, 0),
        energy=None,
        distance=None,
        outofplane_angle=None,
        inplane_angle=None,
        tilt_angle=None,
        rocking_angle=None,
        grazing_angle=None,
        **kwargs,
    ):
        valid.valid_kwargs(
            kwargs=kwargs,
            allowed_kwargs={
                "direct_beam",
                "dirbeam_detector_angles",
                "filtered_data",
                "custom_scan",
                "custom_images",
                "custom_monitor",
                "custom_motors",
                "sample_inplane",
                "sample_outofplane",
                "sample_offsets",
                "offset_inplane",
                "actuators",
                "is_series",
                "template_imagefile",
                "roi",
                "binning",
                "preprocessing_binning",
                "custom_pixelsize",
                "linearity_func",
                "logger",
            },
            name="Setup.__init__",
        )
        # kwarg for logging
        self.logger = kwargs.get("logger", module_logger)

        # kwargs for loading and preprocessing data
        self.dirbeam_detector_angles = kwargs.get("dirbeam_detector_angles")
        self.direct_beam = kwargs.get("direct_beam")
        self.filtered_data = kwargs.get("filtered_data", False)  # boolean
        self.custom_scan = kwargs.get("custom_scan", False)  # boolean
        self.custom_images = kwargs.get("custom_images")  # list or tuple
        self.custom_monitor = kwargs.get("custom_monitor")  # list or tuple
        self.custom_motors = kwargs.get("custom_motors")  # dictionnary
        self.actuators = kwargs.get("actuators", {})  # list or tuple

        # kwargs for xrayutilities, delegate the test on their values to xrayutilities
        self.sample_inplane = kwargs.get("sample_inplane", (1, 0, 0))
        self.sample_outofplane = kwargs.get("sample_outofplane", (0, 0, 1))
        self.offset_inplane = kwargs.get("offset_inplane", 0)

        # kwargs for series (several frames per point) at P10
        self.is_series = kwargs.get("is_series", False)  # boolean

        # create the detector instance
        self.detector_name = detector_name
        self.detector = create_detector(name=detector_name, **kwargs)

        # create the beamline instance
        self.beamline_name = beamline_name
        self.beamline = create_beamline(name=beamline_name, **kwargs)

        # load positional arguments corresponding to instance properties
        self.beam_direction = beam_direction
        self.energy = energy
        self.distance = distance
        self.outofplane_angle = outofplane_angle
        self.inplane_angle = inplane_angle
        self.tilt_angle = tilt_angle
        self.rocking_angle = rocking_angle
        self.grazing_angle = grazing_angle

        # initialize other attributes
        self.logfile = None
        self.incident_angles = None
        