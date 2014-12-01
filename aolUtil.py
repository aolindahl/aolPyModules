# -*- coding: utf-8 -*-
"""
Created on Fri May  2 12:04:54 2014

@author: alindahl
"""

class struct(object):
    def __init__(self, dataDict=None):
        self._memberList = []
        if dataDict is not None:
            self.fromDict(dataDict)

    def members(self):
        print self._memberList

    def toDict(self):
        dict = {}
        for k in self._memberList:
            dict[k] = getattr(self, k)

    def fromDict(self, dataDict):
        for k,v in dataDict.items():
            setattr(self, k, v)
            self._memberList.append(k)

        self._memberList.sort()

