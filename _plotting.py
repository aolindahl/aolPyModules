# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:11:32 2015

@author: antlin
"""


import numpy as _np


def limits_from_centers(centers):
    """Get the bin limits correspondinc to given centers vector"""

    centers = centers.astype(float)
    limits = _np.empty(len(centers)+1)
    limits[0] = centers[0] - _np.diff(centers[0:2])/2
    limits[-1] = centers[-1] + _np.diff(centers[-2:])/2
    limits[1:-1] = centers[:-1] + _np.diff(centers)/2
    return limits


def center_histogram(data, centers):
    """Rerurn histogram vector corresponding to given centers"""

    limits = limits_from_centers(centers)
    hist, _ = _np.histogram(data, limits)
    return hist


def center_histogram_2d(x_data, y_data, x_centers, y_centers=None):
    """Rerurn histogram array corresponding to given centers"""

    x_limits = limits_from_centers(x_centers)
    if y_centers is None:
        y_limits = x_limits
    else:
        y_limits = limits_from_centers(y_centers)
    hist, _, _ = _np.histogram2d(x_data, y_data, [x_limits, y_limits])
    return hist.T
