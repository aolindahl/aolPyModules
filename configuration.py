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
        return {convertToStrings(key): convertToStrings(value) 
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [convertToStrings(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


def load_configuration_dict(file_name, default=None, verbose=False):
    if (not os.path.isfile(fileName)):
        if default is None:
            raise IOError(('The file {} could not be read' +
                ' and no default was given.').format(file_name)
        
        return struct(default)
        
    with open(file_name, 'r') as fp:
        config_dict = json.load(fp)
    
    config_dict = convertToStrings(config_dict)

    if verbose:
        print 'Configuration file read:'
        print json.dumps(config_dict, sort_keys=True, indent=4)

    return config_dict


def loadConfiguration(file_name, default=None, quiet=False):
       return struct(load_config_dict(file_name, default, ~quiet))
