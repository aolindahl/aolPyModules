# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:50:18 2014

@author: hartmang
"""
import numpy as np
#detectors' angles
phi = np.arange(0, 360, 22.5)
phir = np.radians(phi)
#basis vectors
Bt = np.array([ np.ones_like(phir), np.cos(2*phir), np.sin(2*phir)]).T
B = Bt / np.sqrt( (Bt**2).sum(axis=0) )

# angular distribution
# I set the random numbers to fixed parameters in order to compare it with the Matlab Script -fine!
a = np.random.normal(size=1)
b = 1+np.random.normal(size=1)
I0 = (100+200*np.random.normal(size=1))
d = np.radians(360*np.random.normal(size=1))

a = 0.817
b = 1.9
I0 = 125
d = np.radians(328)

Itest0 = I0 * (1 + b * 0.25 * (1 + 3*a*np.cos(2*(phir+d))))
Inoise = np.sqrt(Itest0) * np.random.normal(size=Itest0.shape)
# check: Inoise=[0 for i in range(0,16)] -fine!

Itest = [Itest0[i]+float(Inoise[i]) for i in range(0,16)]

c = B.T.dot(Itest)

#operating status
hassig0 = np.array([1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1], dtype=bool)

#matrix stuff
hassig0b = ~hassig0.astype(bool)

Bn2 = np.diag(hassig0b).dot(B)
Bfir = np.linalg.inv(np.eye(3) - Bn2.T.dot(Bn2))
Bsec = B.T.dot(hassig0*Itest)

cgab = Bfir.dot(Bsec)

print c
print cgab




