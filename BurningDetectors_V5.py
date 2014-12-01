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
        # B is array of basis vectors
        BNotNorm = np.array(
                [ np.ones_like(self._phir),
                    np.cos(2 * self._phir),
                    np.sin(2 * self._phir) ]).T
        self._B = BNotNorm / np.sqrt( (BNotNorm**2).sum(axis=0) )
        self._B = np.array(
                [ np.ones_like(self._phir)/16.,
                    np.cos(2 * self._phir)/8.,
                    np.sin(2 * self._phir)/8. ])

    def setFitMask(self, fitMask=np.ones(16)):
        self._fitMask = fitMask.astype(bool)
        Q = np.diag( self._fitMask.astype(float) )
        Qnot = np.diag( (~self._fitMask.astype(bool)).astype(float) )

        B = self._B
        Bt = B.transpose()
        BtB = Bt.dot(B)
        print 'BtB:', BtB
        BtBinv = np.linalg.inv(BtB)
        print 'BtBinv', BtBinv
        np.allclose(

        BQnot = B.dot(Qnot)

        C = np.eye(3) - BQnot.dot(BtBinv).dot(Bt)

        Cinv = np.linalg.inv(C)

        self._BCorrected = Cinv.dot(B).dot(Q)

        #self._BCorrected = np.linalg.inv( np.eye(3) - self._B.dot( Qnot ).dot(
        #    np.linalg.inv( self._B.T.dot(self._B) ) ).dot( self._B.T ) ).dot(
        #            self._B ).dot( Q )

        #self._fitMask = self._fitMask.astype(int)
        
    def solve(self, data, beta):
        coef = self._BCorrected.dot(data)
        print coef
        return self.params(coef, beta)

    def solveAllThere(self, data, beta):
        coef = self._B.dot(data)
        print coef
        return self.params(coef, beta)

    def params(self, coef, beta):
        I0 = 4. * coef[0] / (4. + beta)
        linear = (4. + beta) * np.sqrt(coef[1]**2 + coef[2]**2) / (3. * beta *
                coef[0])
        tilt = 0.5 * np.arctan2(coef[2], coef[1])

        return I0, linear, tilt 

def model(I0, beta, linear, tilt, angles):
        return I0 * (1. + beta / 4. * (1. + 3. * linear * np.cos( 2. * (angles - d))))

if __name__ == '__main__':
    # angular distribution
    # I set the random numbers to fixed parameters in order to compare it with the Matlab Script -fine!
    a = np.random.random()
    a = 0.5
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
