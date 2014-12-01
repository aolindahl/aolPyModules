import numpy as np
import tofData2 as tofData
from configuration import loadConfiguration as loadConfig
from aolUtil import struct
import sys
import random
import lmfit
from BurningDetectors_V6 import projector

# A bunch of methods to take care of the cookie box data


def modelFunction(params, x, y=None, eps=None):
    A = params['A'].value
    beta = params['beta'].value
    tilt = params['tilt'].value
    linear = params['linear'].value

    mod = A * ( 1 + beta * 0.25 * ( 1 + 3 * linear * np.cos( 2*(x-tilt) ) ) ) 

    if y is None:
        ret = mod
    elif eps is None:
        ret = mod-y
    else:
        ret = (mod-y)/eps

    return ret

def initialParams(yData=None):
    params = lmfit.Parameters()
    params.add('A', 10, min=0)
    params.add('beta', 2, min=-1, max=2)
    params.add('tilt', 0, min = -np.pi/2, max=np.pi/2)
    params.add('linear', 0.5, min=0)
    #params.add('tilt', np.pi, min=-2*np.pi, max=2*np.pi)
    #params.add('linear', 0.5, min=-0.5, max=1.5)

    if yData!=None:
        params['A'].value = yData[np.isfinite(yData)].mean()
        tilt = initialAngle(yData)
        #params['tilt'].value = tilt
        #params['tilt'].min = tilt - 2*np.pi
        #params['tilt'].max = tilt + 2*np.pi

    return params

_anglesDeg = np.arange(0, 360, 22.5)
_anglesRad = _anglesDeg * np.pi / 180
_sinVect = np.sin( 2*_anglesRad )
_cosVect = np.cos( 2*_anglesRad )

def initialAngle(yData):
    vec = yData.copy()
    vec[np.isnan(vec)] = 0
    A = yData.dot(_sinVect)
    B = yData.dot(_cosVect)
    phi = 0.5 * ( np.arctan(A/B)
            + np.pi * ( -1 * ((A<0) & (B<0)) + 1 * ((A>0) & (B<0)) ) )
    return phi



def sliceFromRange(data, range):
    return [slice(d.searchsorted(np.min(range)), d.searchsorted(np.max(range)))
            for d in data]


class CookieBox:
    'Class that handels the cookiebox data'
    def __init__(self, config, verbose=False):
        '''\
        Initialization method.\
        
        The configuration should be a single object or a list of 16 objects.\
        '''
        # A list of the tofData objects
        self._tofList = []
        for conf in config.tofConfigList:
            self._tofList.append(tofData.tofData(conf, quiet=not verbose))


        self._phiDeg = _anglesDeg
        self._phiRad = _anglesRad

        self._timeAmplitudesUpToDate = False

        self.proj = projector()

    def setupScales(self, env=None,
            energyAxisBinLimits=np.linspace(800, 1000, 129),
            timeScale=None, verbose=False):
        for  tof in self._tofList:
            tof.setupScales(env, energyAxisBinLimits, timeScale)

    def setBaselineSubtractionAveraging(self, weightLast):
        for tof in self._tofList:
            tof.setBaselineSubtractionAveraging(weightLast)

    def setRawData(self, evt, verbose=False, newDataFactor=1):
        for tof in self._tofList:
            tof.setRawData(evt, newDataFactor=newDataFactor)
        self._timeAmplitudesUpToDate = False

    def getTimeScales_us(self, roi=None, verbose=False):
        return [tof.getTimeScale_us(roi=roi) for tof in self._tofList]

    def getTimeAmplitudes(self, roiSlices=None, verbose=False):
        if not self._timeAmplitudesUpToDate:
            self._timeAmplitudes = [tof.getTimeAmplitude() for tof in
                    self._tofList]
        if roiSlices is None:
            return self._timeAmplitudes
        return [amp[s] for amp, s in zip(self._timeAmplitudes, roiSlices)]

    def getTimeAmplitudesFiltered(self, roi=None, verbose=False):
        return [tof.getTimeAmplitudeFiltered(roi=roi) for tof in self._tofList]

    def getEnergyScales_eV(self, roi=None, verbose=False):
        return [tof.getEnergyScale_eV(roi=roi) for tof in self._tofList]
    
    def getEnergyAmplitudes(self, roi=None, verbose=False):
        return [tof.getEnergyAmplitude(roi=roi) for tof in self._tofList]

    def getMoments(self, domain='Time', roi=None):
        return [ tof.getMoments(domain=domain, roi=roi) for tof in
                self._tofList]
    
    def getPhotonEnergy(self, energyShift=0):
        moments = np.array(self.getMoments(domain='Energy'))
        amplitudes = self.getIntensityDistribution(domain='Energy')

        return (np.average(moments, axis=0, weights=amplitudes)
                + np.array([energyShift, 0]))


    def getIntensityDistribution(self, roiSlices=[slice(None)]*16,
            domain='Time', verbose=False, detFactors=[1]*16):
        if verbose:
            print 'Geting the intensity distribution',
            if roi is None:
                print '.'
            else:
                print 'in roi {}.'.format(roi)
        intensity = []
        if domain=='Energy':
            if verbose:
                print 'Using energy domain.'
            ampFunk = self.getEnergyAmplitudes
        else:
            ampFunk = self.getTimeAmplitudesFiltered
        for amp, factor, s in zip(ampFunk(verbose=verbose), detFactors,
                roiSlices):
            if amp is None:
                intensity.append(np.nan)
            else:
                intensity.append(amp[s].sum() * factor)

        if verbose:
            print 'Returning vector of length {}.'.format(len(intensity))

        return np.array(intensity)

    def getAngles(self, kind='rad'):
        if kind=='rad':
            return self._phiRad
        else:
            return self._phiDeg

    def randomizeAmplitudes(self, verbose=False):
        '''This is only for testing. Direct manipulation of the private
        variables of the tofData class (as done here) is not recommended.'''
        params = initialParams()
        params['linear'].value = random.random()
        params['tilt'].value = random.random()*np.pi - np.pi/2
        params['A'].value = 1
        if verbose:
            for par in params.itervalues():
                print par
        factors = ( modelFunction(params, self.getAngles()) *
                np.random.normal(1, 0.05, (16,)) )
        #factors = modelFunction(params, self.getAngles())
        for factor, tof in zip(factors, self._tofList):
            tof._timeAmplitude *= factor
            #tof._timeAmplitudeFiltered *= factor
            tof._energyAmplitude *= factor

        return params
        



if __name__ == '__main__':
    import psana
    import time
    import matplotlib.pyplot as plt
    plt.ion()
    
    ds = psana.DataSource('exp=amoi0114:run=33')


    #config = loadConfig('cookieBoxDefaultConfig.json')
    #config = [config]*16
    import CookieBox3DefaultConfig as config
    reload(config)
    #config.makeTofConfigList(online=False)

    config

    if hasattr(config, 'domainToDisplay') and config.domainToDisplay=='Energy':
        domain = 'Energy'
    else:
        domain = 'Time'
    cb = CookieBox(config, verbose=False)

    #randTilt = []
    #projTilt = []
    t_temp = None
    for i, evt in enumerate(ds.events()):
        if t_temp is not None:
            print 'event processing time:', time.time()-t_temp, 's'
        t_temp = time.time()
        if i >= 1:
            break
        if i%10==0:
            plot=True
        else:
            plot=False
        if i==0:
            cb.setupScales(ds.env(), config.energyScaleBinLimits)
            
            energyScales = cb.getEnergyScales_eV()
            #xScalesRoi0 = cb.getEnergyScales_eV(roi=0)
            #xScalesRoi1 = cb.getEnergyScales_eV(roi=1)
            timeScales = cb.getTimeScales_us()
            tRoi0s = sliceFromRange(timeScales, config.timeRoi0_us)
            tRoi1s = sliceFromRange(timeScales, config.timeRoi1_us)
            timeScalesRoi0 = [t[s] for t,s in zip(timeScales, tRoi0s)]
            timeScalesRoi1 = [t[s] for t,s in zip(timeScales, tRoi1s)]
            
            fig1 = plt.figure(1); plt.clf()
            for k, x , xRoi0, xRoi1 in zip(range(16), timeScales,
                    timeScalesRoi0, timeScalesRoi1):
                plt.subplot(4,4,k+1)
                plt.title('det at {} deg'.format(cb.getAngles('deg')[k]))
                plt.plot(x, np.zeros_like(x))
                plt.plot(xRoi0, np.zeros_like(xRoi0), 'r')
                plt.plot(xRoi1, np.zeros_like(xRoi1), 'g')
                yTMin, yTMax = 0, 0

            fig2 = plt.figure(2);
            fig2.clf()
            fig2ax = fig2.add_subplot(111, polar=True)
            angles = cb.getAngles('rad')
            anglesFit = np.linspace(0, 2*np.pi, 10000)
            fig2ax.plot(
                    angles, np.ones_like(angles), 'ro',
                    angles, np.ones_like(angles), 'gs',
                    anglesFit, np.zeros_like(anglesFit), 'm-')
                    #anglesFit, np.zeros_like(anglesFit), 'm--')

            fig3 = plt.figure(3)
            fig3.clf()
            fig3ax = fig3.add_subplot(111)
            fig3ax.plot(
                    energyScales[0], np.zeros_like(energyScales[0]))
        print 'event number', i


        cb.setRawData(evt)
        #randParams = cb.randomizeAmplitudes(verbose=True)
        #randTilt.append(randParams['tilt'].value) 

        if plot:
            amplitudes0 = cb.getIntensityDistribution(roi=0, domain=domain,
                verbose=True)
            amplitudes1 = cb.getIntensityDistribution(roi=1, domain=domain,
                verbose=True)
            
            energyData = cb.getEnergyAmplitudes()
            print len(energyData)
            print energyData[0].shape
            #spectrum = np.average(energyData, axis=0)
            spectrum = energyData[15]
            print spectrum
            #energyDataRoi0 = cb.getEnergyAmplitudes(roi=0)
            #energyDataRoi1 = cb.getEnergyAmplitudes(roi=1)
            timeData = cb.getTimeAmplitudes()
            timeDataRoi0 = [t[s] for t,s in zip(timeData, tRoi0s)]
            timeDataRoi1 = [t[s] for t,s in zip(timeData, tRoi1s)]
            
            tMin = np.min(timeData)
            tMax = np.max(timeData)
            rescaleFlag = False
            if tMin < yTMin:
                yTMin = tMin
                rescaleFlag = True
            if tMax > yTMax:
                yTMax = tMax
                rescaleFlag = True
        
            for ax, y, yRoi0, yRoi1 in zip(fig1.axes, timeData, timeDataRoi0,
                    timeDataRoi1):
                ax.lines[0].set_ydata(y)
                ax.lines[1].set_ydata(yRoi0)
                ax.lines[2].set_ydata(yRoi1)
                if rescaleFlag:
                    ax.set_ybound(yTMin, yTMax)

        amplitudes = cb.getIntensityDistribution(roi=0, domain=domain,
                verbose=True)
        params = initialParams(amplitudes)
        #projTilt.append(params['tilt'].value)
        params['beta'].vary=False
        res = lmfit.minimize(modelFunction, params, args=(angles, amplitudes),
                method='leastsq')
        print res.nfev, 'function evaluations'
        print 'Fit', ('succeded' if res.success else 'failed')
        print res.message
        
        print lmfit.fit_report(params)


        moments = cb.getMoments(domain=domain, roi=0)

        if plot:
            fig2ax.lines[0].set_ydata(amplitudes0)
            fig2ax.lines[1].set_ydata(amplitudes1)
            fig2ax.lines[2].set_ydata(modelFunction(params, anglesFit))
            #fig2ax.lines[3].set_ydata(amplitudes0.mean() * modelFunction(randParams, anglesFit))

            fig2ax.relim()
            fig2ax.autoscale_view()


            fig3ax.lines[0].set_ydata(spectrum)
            fig3ax.relim()
            fig3ax.autoscale_view()

            for fig in [fig1, fig2, fig3]:
                fig.canvas.draw()


    #randTilt = np.array(randTilt)
    #projTilt = np.array(projTilt)
    raw_input('Press enter to exit...')
