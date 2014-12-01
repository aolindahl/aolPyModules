# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

# -*- coding: utf-8 -*-
"""
Created on Thu May  1 11:04:42 2014

@author: alindahl
"""

import zmq

class line:
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y

class linePlot:
    def __init__(self, lines=None):
        if lines == None:
            self.lines = []
        else:
            self.lines = lines

class zmqSender:
    """\
    Class that takes care of the zmq comunication of data to the client.\
    """
    def __init__(self, socketNumber=19820809, hwm=1):
        self._context = zmq.Context()
        self._socket = self._context.socket(zmq.PUB)
        self._socket.set_hwm(hwm)
        #self._socket.setsockopt(zmq.SNDHWM, hwm)
        self._socket.bind("tcp://*:%d" % socketNumber)
        
    def __del__(self):
        self._socket.close()
        
    def sendObject(self, obj, topic='data'):
        self._socket.send(topic, zmq.SNDMORE)
        self._socket.send_pyobj(obj)
