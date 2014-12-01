# -*- coding: utf-8 -*-
"""
Created on Fri May  9 09:55:49 2014

@author: hnick
"""
from pylab import median
import pywt

def stand_mad(x):
    """Calculates median average deviation"""
    return 1.4826*median(abs(x-median(x)))
    
def wavelet_filt(data, thresh=None, W='db5', levels=6, printTh=False):    
    """ Array filtering based on hard-threshholding of wavelet coefficients"""
    arraywc = pywt.wavedec(data,W,level=levels,mode='sym')
    if thresh == None:
        thresh = 4*stand_mad(arraywc[-1])
    if printTh:
        print thresh
    arraywcfilt = map(lambda x: pywt.thresholding.hard(x,thresh), arraywc)
    return pywt.waverec(arraywcfilt, W)


