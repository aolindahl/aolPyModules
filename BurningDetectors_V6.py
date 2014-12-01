# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:50:18 2014

@author: hartmang
"""
import numpy as np
class projector:
    def __init__(self):
        #detectors' angles
        self._phir = np.radians(np.arange(0, 360, 22.5))

        self._norm = np.diag( [1./4., 1./np.sqrt(8.), 1./np.sqrt(8.)] )

        #basis vectors
        BNotNorm = np.array(
                [ np.ones_like(self._phir),
                    np.cos(2 * self._phir),
                    np.sin(2 * self._phir) ]).T
        self._B = BNotNorm / np.sqrt( (BNotNorm**2).sum(axis=0) )

        self._BAll = self._norm.dot( self._B.T )
        self._BGap = self._B

        self.setFitMask()

    def setFitMask(self, fitMask=np.ones(16)):
        self._fitMask = fitMask.astype(bool)

        self._BGap = np.linalg.inv( np.eye(3) -
                self._B[~self._fitMask,:].T.dot(
                    self._B[~self._fitMask,:] ) ).dot(
                            self._B[self._fitMask,:].T )
        #normalization
        self._BGap = self._norm.dot( self._BGap )

        #self._fitMask = self._fitMask.astype(int)

    def solve(self, data, beta):
        coef = self._BGap.dot(data[self._fitMask])
        return self.params(coef, beta)

    def solveAllThere(self, data, beta):
        coef = self._BAll.dot(data)
        return self.params(coef, beta)

    def params(self, coef, beta):
        I0 = 4. * coef[0] / (4. + beta)
        linear = (4. + beta) * np.sqrt(coef[1]**2 + coef[2]**2) / (3. * beta *
                coef[0])
        tilt = 0.5 * np.arctan2(coef[2], coef[1])

        if linear < 0:
            linear = -linear
        while tilt < -np.pi/2:
            tilt += pi
        while tilt > np.pi/2:
            tilt -= np.pi

        return I0, linear, tilt 

def model(I0, beta, linear, tilt, angles):
        return I0 * (1. + beta / 4. * (1. + 3. * linear * np.cos( 2. * (angles -
            tilt))))

if __name__ == '__main__':
    # angular distribution
    # I set the random numbers to fixed parameters in order to compare it with the Matlab Script -fine!
    a = np.random.random()
    a = 1
    #b=1+np.random.normal(size=1)
    b = 2.
    #A = np.random.random() * 200 + 100
    A = 100
    d = np.radians(180 * np.random.random() - 90)
    d = np.radians(22.5)
    
    #a=0.817
    #b=1.9
    #I0=125
    #d=np.radians(328)
    
    angles = np.radians(np.arange(0,360,22.5))

    Itest0 = model(I0=A, beta=b, linear=a, tilt=d, angles=angles)
    Inoise = np.sqrt(Itest0)*np.random.normal(size=16)
    Inoise = Inoise * (Inoise > 0)
    # check: Inoise=[0 for i in range(0,16)] -fine!
    
    Itest = Itest0 + Inoise
    
    
    #operating status
    mask = np.array([1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1], dtype=bool)
    
    proj = projector()
    proj.setFitMask(mask)

    c = proj.solveAllThere(Itest0, b)
    cMissing = proj.solve(Itest0, b)
    
    print [A, a, d]
    print c
    print cMissing

    import matplotlib.pyplot as plt
    plt.ion()
    
    fig = plt.figure(1); plt.clf()
    fig.add_subplot(111, polar=True)
    plt.plot(angles, Itest0)
    plt.plot(angles, Itest)
    plt.plot(angles, model(I0=cMissing[0], beta=b, linear=cMissing[1],
        tilt=cMissing[2], angles=angles))
