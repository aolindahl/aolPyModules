# -*- coding: utf-8 -*-
"""
Created on Thu May  1 17:27:33 2014

@author: alindahl
"""

import json
import os.path

from aolUtil import struct

def convertToStrings(input):
    if isinstance(input, dict):
        return {convertToStrings(key): convertToStrings(value) for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [convertToStrings(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

def loadConfiguration(fileName, default=None, quiet=False):
    if (not os.path.isfile(fileName)):
        if default == None:
            raise IOError(
            "The file %s could not be read and no default was given." % fileName)
            exit()
        else:
            return struct(default)
            
    f = open(fileName, 'r')
    configDict = json.load(f)
    f.close()
    configDict = convertToStrings(configDict)
    if not quiet:
        print "Configuration file read:\n\
%s" % json.dumps(configDict, sort_keys=True, indent=4)

    return struct(configDict)
