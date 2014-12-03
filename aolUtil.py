# -*- coding: utf-8 -*-
"""
Created on Fri May  2 12:04:54 2014

@author: alindahl
"""

class struct(object):
    def __init__(self, dataDict=None):
        if dataDict is not None:
            self.fromDict(dataDict)

    def members(self):
        print [k for k in self.__dict__]

    def toDict(self):
        return self.__dict__

    def fromDict(self, dataDict):
        for k,v in dataDict.items():
            setattr(self, k, v)

    def __repr__(self):
        s = 'aolUtil.struct( {'
        for k,v in self.__dict__.iteritems():
            s += '\n\t"{}":{}'.format(k,v)
        s += ' })'
        return s

