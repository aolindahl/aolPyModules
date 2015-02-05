#from setupEnvironment import *
import numpy as np
from configuration import loadConfiguration as loadConfig
try:
    import psana
except ImportError:
    psana = None
    print 'Using the module "{}" without psana capabililties.'.format(__name__)
import time
import wiener
from scipy.sparse import coo_matrix
import sys


_useWavelet = True
if _useWavelet:
    try:
        from wavelet_filter import wavelet_filt as wavelet
    except:
        print 'Wavelet filtering not avaliable. pywt not found.'
        _useWavelet = False

m_e_eV = 0.510998928e6 # http://physics.nist.gov/cgi-bin/cuu/Value?me 2014-04-21
c_0_mps = 299792458 # http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in! 2014-04-21

###################################################
# Functions interacting with psana based data



def getSignalScaling(env, sourceString, channel, verbose=False):
    """Get information on how to rescale the raw signal.

    env -- psana environment.
    sourceString
        String describing the psana dataSource of the aquiris board(s).
    channel
        The channel to get data for.
    verbose
        If True the function prints debug information. Default = False.

    Returns: scaling, offset

    scaling
        Multiplicative factor (V).
    offset
        Offset (V).
    """
    # Get the configuration
    try:
        acqirisConfig = env.configStore().get(psana.Acqiris.ConfigV1,
                getSource(sourceString) )
    except:
        return None

    # Get the scaling constants for the vertical scale.
    # convenience reference
    vertScale = acqirisConfig.vert()[channel]
    # The vertical scale information is given as the full scale voltage over
    # all the 2**16 bits.
    # Here the voltage per bit is calculated
    scaling = vertScale.fullScale() * 2**-16
    # The scale also has an offset in voltage
    offset = vertScale.offset()

    return scaling, offset



def getRawSignals(evt, ch, sourceString):
    """Get raw aquiris data.

    evt
        (psana.Event) Event to extract data from.
    ch
        (Int) Aquiris channel to get the data from.
    souceString
        (String) Name of data source in psana.
    """

    try:
        # try to get the acqiris data
        acqirisData = evt.get(psana.Acqiris.DataDescV1, getSource(sourceString))
    except:
        return None

    if acqirisData == None:
        return None

    return acqirisData.data(ch).waveforms()[0]

class tofData(object):
    """Class to handle TOF data.

    Conversions to energy scale including rescaling of the amplitudes\
    """

    def __init__(self, config, quiet=True):
        """\
        Initialize the tofData class giving it a configutaion object.\
        """

        # Extract the data source
        self._source = psana.Source(config['detectorSource'])
        # Extract the acqiris channel
        self._acqCh = config['acqCh']
        # Load the callibration file pointed to in the configuration
        self._calib = loadConfig(config['calibFile'], quiet=quiet)

        # Store the configuration
        self._config = config

        # Setup the basic rescaling
        self._acqVertScaling = 1
        self._acqVertOffset = 0

        # Initialy the class does not contain any data or any scales
        self._noData = True
        self._noScales = True

        # Basic info about filtering
        self._filterTimeDomain = False
        self._timeAmplitudeFiltered = None

        self._quiet = quiet

        if 'filterMethod' in self._config.keys():
            if self._config['filterMethod'] == "wienerDeconv":
                self.filterTimeDomainData(method='wienerDeconv',
                        SNR=self._config['filterWienerSNR'],
                        response=self._config['filterWienerResponse'])
            elif self._config['filterMethod'] == 'wavelet':
                self.filterTimeDomainData(method='wavelet',
                        levels=self._config["filterWaveletLevels"])
            elif self._config['filterMethod'] == 'average':
                if not self._quiet:
                    print 'Using averaging with {} points'.format(
                            self._config['filterAverageNumPoints'])
                self.filterTimeDomainData(method='average',
                        numPoints=self._config['filterAverageNumPoints'])


        self._timeAmplitude = None
        self._timeRoiSlice = [None, None]
        self._energyRoiSlice =[None, None]

        self._bgWeight = None



    def setupScales(self, env=None,
                    energyAxisBinLimits=np.linspace(800, 1000, 129),
                    timeScale=None, retardation=0):
        """Setup the information about the acqiris channel used for the TOF.

        Reads the scale factors for the raw aquiris data and also calculates
        the conversion to the energy domain.
        """

        if not self._quiet:
            print 'Seting up the scales.'

        segment = 0 # Somehow the aquiris channels suses only segment 0

        # Try to get the acqiris configuration form the environment object
        try:
            acqirisConfig = env.configStore().get(psana.Acqiris.ConfigV1)
        except:
            acqirisConfig = None

        # If the returned configuration is None tnd here is no valide scale
        # either
        if acqirisConfig is None and timeScale is None:
            print ('WARNING: No acqiris configuration obtained,\n' +
                    ' and no time scale given explicitly. WARNING')
            self._noScales = True
            return

        ################
        # Seting up the scales in time domain

        if timeScale is None:
            # Make the time scale vector for the acqiris channel.
            # This is just for convenience
            timeScale = acqirisConfig.horiz()
            # Start time
            t0 = timeScale.delayTime()
            # Time step
            dt = timeScale.sampInterval()
            # Number of samples
            nSample = timeScale.nbrSamples()
            # Make the time scale vector from the above information and rescale it
            # to microseconds
            self._timeScale_us = np.arange(t0, dt*nSample, dt)*1e6

            # If only a slice of the time should be used
            if self._config['tSlice']:
                # Define the time slice for the configuration
                self._timeSlice = slice(
                        np.searchsorted(self._timeScale_us,
                            self._config['tMin_us']),
                        np.searchsorted(self._timeScale_us,
                            self._config['tMax_us']))
                self._timeScale_us = self._timeScale_us[self._timeSlice]
            else:
                self._timeSlice = slice(None)


            # Get the scaling constants for the vertical scale.
            # convenience reference
            vertScale = acqirisConfig.vert()[self._acqCh]
            # The vertical scale information is given as the full scale voltage over
            # all the 2**16 bits.
            # Here the voltage per bit is calculated
            self._acqVertScaling = vertScale.fullScale() * 2**-16
            # The scale also has an offset in voltage
            self._acqVertOffset = vertScale.offset()
        else:
            self._timeScale_us = timeScale.copy()
            dt = np.diff(self._timeScale_us).mean()
            t0 = self._timeScale_us[0]
            nSample = len(timeScale)
            self._timeSlice = slice(None)

        #print 'dt =', dt

        # Calculate the bin edges in the time domain
        tEdges = np.concatenate([self._timeScale_us - dt/2,
            self._timeScale_us[-1:] + dt/2])
        # and the corresponding bin edges in the energy domain.
        # The conversion is given by:
        # E = D^2 m_e 10^6 / ( 2 c_0^2 (t - t_p)^2 ) + E0,
        # where:
        # D is the flight distance in mm
        # m_e is the electon rest mass expressed in eV
        # the 10^6 factor accounts the the otherwhise missmachching prefixes
        # c_0 is the speed of light in m/s
        # E is the energy in eV
        # E0 is an energy offset in eV, should be determined in a callibration
        # fit
        # t_p is the arrival time of the prompt signal in microseconds
        eEdges = (self._calib.D_mm**2 * m_e_eV * 1e6
                / (c_0_mps**2 * 2 * (tEdges - self._calib.prompt_us
                    - self._calib.tOffset_us)**2)
                + self._calib.EOffset_eV - self._calib.EOffsetRetardation *
                retardation)



        # Number of energy bins
        nEBins = len(energyAxisBinLimits) - 1
        # Energy bin size
        dE = np.diff(energyAxisBinLimits).mean()

        # Make matrixes out of the edges vectors in energy domain
        Me = np.concatenate([eEdges.reshape(1,-1)]*nEBins)
        ME = np.concatenate([energyAxisBinLimits.reshape(-1,1)] *
                len(self._timeScale_us), axis=1)

        # Compute the start and end energies for the conversion from the time axis
        # to energy axis
        highE = ( np.minimum( Me[:,:-1], ME[1:,:] ) )
        lowE = ( np.maximum( Me[:,1:], ME[:-1,:] ) )
        # Only where the high energy is more than the low energy the conversion makes anny sense
        I = lowE < highE
        # Allocate a tempoaraty conversion matrix
        tempMat = np.zeros((nEBins, len(self._timeScale_us)))
        # Calculate the elements of the conversion matrix
        # For each time bin the measured amplitude is multiploed by the bin size
        # in order to errive at the integral of the signal. Then the it is
        # determined how much of each time bin contributes to each energy bin.
        # This is done by comparing the edge positions of the time and energy
        # bins and assigning the correct proportion of the integral in time
        # domain to integral in the energy domain. Finally the total integral is
        # divided by the energy bin size in order to return to an amplitude.
        # Summation over all time bins is performed in the matrix multiplication
        # of the conversion matrix with the time domain amplitude vector.
        tempMat[I] = dt * (highE[I] - lowE[I]) / ( Me[:,:-1] - Me[:,1:] )[I] / dE
        # The conversion matrix is highly sparse, thus make a sparse matrix to
        # spped up the calculationss
        self._timeToEnergyConversionMatrix = coo_matrix(tempMat)*1e9


        #plt.figure(111); plt.clf()
        #plt.imshow(self._timeToEnergyConversionMatrix.todense(),
        #        interpolation='none', aspect='auto',
        #        extent=(self._timeScale_us.min(), self._timeScale_us.max(),
        #            energyAxisBinLimits.min(), energyAxisBinLimits.max()),
        #        origin='lower' )
        #plt.gcf().show()




        # Create the energy scale

        # Find out which part of the time scale is after the prompt peak
        # The calibration information is used for this
        I = self._timeScale_us > self._calib.prompt_us + self._calib.tOffset_us
        # Allocate an energy scale with -1 value (unphysical).
        self._rawEnergyScale_eV = -np.ones(self._timeScale_us.shape)

        # Calculate the energy in eV that corresponds to all times after the
        # time for prompt signal. This is done as:
        # E = D**2 * m_e / (c_0**2 * (t - t_prompt - t_fitOffset)**2) *1e6
        #   - E_fitOffset
        # The 1e6 factor accounts for the unit mismatch
        self._rawEnergyScale_eV[I] = (self._calib.D_mm**2
            * m_e_eV * 1e6 / (c_0_mps**2 * 2
                * (self._timeScale_us[I] - self._calib.prompt_us
                    - self._calib.tOffset_us)**2)
            - self._calib.EOffset_eV)
        # Set the Energies that correspond to imes befor the prompt to the
        # maximum energy.
        self._rawEnergyScale_eV[~I] = np.max(self._rawEnergyScale_eV)


        # Make the energy scale from the bin limits. Calculating the center
        # of the bins.
        self._energyScale_eV = (np.diff(energyAxisBinLimits)/2
                                + energyAxisBinLimits[:-1])

        # get the energy bin size
        self._energyBinSize = self._energyScale_eV[1] - self._energyScale_eV[0]

        # Set the region of interest slices
        if not self._quiet:
            print 'Looking for ROIs'
        for domain, roiBase in zip(
                ['Time', 'Energy'],
                ['timeRoi{}_us', 'energyRoi{}_eV']):
            for iRoi in range(2):
                roi = roiBase.format(iRoi)
                if roi in self._config:
                    if not self._quiet:
                        print '\t{} found'.format(roi)
                        print self._config[roi]
                    self.setBaseRoi(
                            min = self._config[roi][0],
                            max = self._config[roi][1],
                            roi = iRoi,
                            domain = domain)


        # Make a backgeound slice
        if self._config['baselineSubtraction'] == 'early':
            self._bgSlice = slice(self._timeScale_us.searchsorted(
                        self._config['baselineEnd_us']))
        elif self._config['baselineSubtraction'] == 'roi':
            try:
                self._bgSlice = slice(
                        min( [i.stop for i in self._timeRoiSlice] ),
                        max( [i.start for i in self._timeRoiSlice] ) )
            except:
                print "Could not set the gsSlice from the roi's."
                print "Attempted slice({}, {})".format(
                        min( [i.stop for i in self._timeRoiSlice] ),
                        max( [i.start for i in self._timeRoiSlice] ) )
                self._bgSlice = slice(0,0)
        else:
                self._bgSlice = slice(0,0)

        # Check if there is actually something to calculate the background
        # from
        if len(self._timeScale_us[self._bgSlice]) < 1:
            self._config['baselineSubtraction'] = 'none'



        self._noScales = False

    def setBaseRoi(self, min, max, roi=0, domain='Time'):
        if domain=='Time':
            a = self._timeScale_us.searchsorted(min)
            b = a + self._timeScale_us[a:].searchsorted(max)
            self._timeRoiSlice[roi] = slice(a,b)
        elif domain=='Energy':
            a = self._energyScale_eV.searchsorted(min)
            b = a + self._energyScale_eV[a:].searchsorted(max)
            self._energyRoiSlice[roi] = slice(a,b)

    def setBaselineSubtractionAveraging(self, weightLast):
        self._bgWeightLast = weightLast
        self._bgWeightHistory = 1-weightLast
        self._bgHistory = None


    def setRawData(self, evt=None, timeAmplitude_V=None, newDataFactor=None):
        """\
        Set waveform data and compute the scaled data.\
        """

        if not self._quiet:
            print 'Seting raw data.'


        # If a psan event was passed as a parametere
        if evt is not None:
            if not self._quiet:
                print 'Event object given.'
            try:
                # try to get the acqiris data
                acqirisData = evt.get(psana.Acqiris.DataDescV1, self._source)
            except:
                self._noData = True
                return

            if acqirisData == None or self._noScales:
                self._noData = True
                return

            # If everything worked out calculate the rescaled amplitude
            new = -(acqirisData.data(self._acqCh).waveforms()[0][self._timeSlice]
                    * self._acqVertScaling - self._acqVertOffset)

            if (self._timeAmplitude is None or newDataFactor == 1 or
                    newDataFactor is None):
                self._timeAmplitude = new
            else:
                self._timeAmplitude = (self._timeAmplitude * (1.0-newDataFactor)
                        + new * newDataFactor)

            # If set up to do that, subtract the baseline
            if self._config['baselineSubtraction'] != 'none':
                if self._bgWeight is None:
                    self._timeAmplitude -= \
                            self._timeAmplitude[self._bgSlice].mean()
                else:
                    if self._bgHistory is None:
                        self._bgHistory =  \
                                self._timeAmplitude[self._bgSlice].mean()
                    else:
                        self._bgHistory *= self._bgWeightHistory
                        self._bgHistory += self._bgWeightLast \
                                * self._timeAmplitude[self._bgSlice].mean()

                    self._timeAmplitude -= self._bgHistory

        elif timeAmplitude_V is  None:
            if not self._quiet:
                print 'Niether event nor time scale given.'
            # If neither psana event, nor amplitude was given
            self._noData = True
            return

        else:
            # If no psana event was given but threre is an amplitude
            self._timeAmplitude = timeAmplitude_V[self._timeSlice]

        self._filterTimeDomainData()

        self.calcEnergyAmplitude()

        self._noData = False

    def calcEnergyAmplitude(self, filtered=None):
        if self._filterTimeDomain and filtered is not False:
            tAmp = self._timeAmplitudeFiltered
        else:
            tAmp = self._timeAmplitude

        # Calculate the signal amplitude in the energy domain.
        self._energyAmplitude = self._timeToEnergyConversionMatrix.dot(tAmp)
        return


    def getTimeScale_us(self, roi=None):
        if self._noScales:
            return None
        if roi!=None and self._timeRoiSlice!=None:
            return self._timeScale_us[self._timeRoiSlice[roi]]
        return self._timeScale_us

    def getTimeAmplitude(self, roi=None):
        if self._noData:
            return None
        if roi!=None and self._timeRoiSlice!=None:
            return self._timeAmplitude[self._timeRoiSlice[roi]]
        return self._timeAmplitude

    def getTimeAmplitudeFiltered(self, roi=None):
        if self._noData:
            return None
        if roi!=None: #and self._timeRoiSlice!=None:
            return self._timeAmplitudeFiltered[self._timeRoiSlice[roi]]
        return self._timeAmplitudeFiltered

    def getEnergyAmplitude(self, roi=None):
        if self._noData:
            return None
        if self._noData:
            return None
        if roi!=None and self._energyRoiSlice!=None:
            return self._energyAmplitude[self._energyRoiSlice[roi]]
        return self._energyAmplitude

    def getEnergyScale_eV(self, roi=None):
        if self._noScales:
            return None
        if roi!=None and self._energyRoiSlice!=None:
            return self._energyScale_eV[self._energyRoiSlice[roi]]
        return self._energyScale_eV

    def getRawEnergyScale(self, roi=None):
        if self._noScales:
            return None
        if roi!=None and self._timeRoiSlice!=None:
            return self._energyScale_eV[self._timeRoiSlice[roi]]
        return self._rawEnergyScale_eV

    def filterTimeDomainData(self, method='wienerDeconv', numPoints=4, levels=6,
            SNR=1, response=None):
        if method is False:
            self._filterTimeDomain = False
            return

        self._filterMethod = method
        self._filterNumPoints = numPoints
        self._filterLevels = levels
        self._filterTimeDomain = True

        if method == 'wienerDeconv':
            if type(SNR) == str:
                self._SNR = np.fromfile(SNR)
            else:
                self._SNR = SNR
            if type(response) == str:
                self._response = np.fromfile(response)
            else:
                self._response = response


    def _filterTimeDomainData(self):
        if self._filterTimeDomain is False:
            self._timeAmplitudeFiltered = self._timeAmplitude
            return

        if self._filterMethod == 'average':
            self._timeAmplitudeFiltered = \
                    np.convolve(self._timeAmplitude,
                            np.ones((self._filterNumPoints,))
                            / self._filterNumPoints, mode='same')
            return

        if self._filterMethod == 'wavelet' and _useWavelet:
            #print 'Wavelet filteing'
            self._timeAmplitudeFiltered = wavelet(self._timeAmplitude,
                    levels=self._filterLevels)
            return

        if self._filterMethod == 'wienerDeconv':
            self._timeAmplitudeFiltered = wiener.deconvolution(
                    self._timeAmplitude, self._SNR, self._response)
            return

        print '{} is not a valid filtering method when "_useWavelet" is {}.'\
                .format(self._filterMethod, _useWavelet)


    def getTraceBounds(self, threshold_V=0.02, minWidth_eV=2, EnergyOffset=0, useRel=False,
            threshold_rel=0.5):
        if self._noData:
            return [np.nan for i in range(3)]

        if useRel:
            threshold_temp = threshold_rel * \
            np.max(self._energyAmplitude[np.isfinite(self._energyAmplitude)])
            if threshold_temp < threshold_V:
                return [np.nan for i in range(3)]
            else:
                threshold_V = threshold_temp
        nPoints = np.round(minWidth_eV/self._energyBinSize)

        min = 0
        for i in range(1, self._energyAmplitude.size):
            if self._energyAmplitude[i] < threshold_V:
                min = i
                continue
            if i-min >= nPoints:
                break
        else:
            min = np.nan


        max = self._energyAmplitude.size - 1
        for i in range(self._energyAmplitude.size-1, -1, -1):
            if self._energyAmplitude[i] < threshold_V:
                max = i
                continue
            if max-i >= nPoints:
                break
        else: max = np.nan

        if min == 0 and max == self._energyAmplitude.size - 1:
            min = np.nan
            max = np.nan

        if not np.isnan(min):
                min = self._energyScale_eV[min] - EnergyOffset
                max = self._energyScale_eV[max] - EnergyOffset

        return min, max, threshold_V



    def getPulseDuration(self, lo, hi):
        if self._noData:
            return None
        amplitude = self._config.tof_maxStreaking_eV
        cutoff = self._config.tof_streakingCutoff_eV
        if hi > cutoff or lo < -cutoff:
            return None
        dur = (np.arccos(lo/amplitude) - np.arccos(hi/amplitude)) / np.pi * \
            self._config.tof_quarterCycle_fs

        return dur

    def getMoments(self, domain='Time', roi=None):
        if domain == 'Time':
            x = self.getTimeScale_us(roi=roi)
            y = self.getTimeAmplitudeFiltered(roi=roi)
        elif domain == 'Energy':
            x = self.getEnergyScale_eV(roi=roi)
            y = self.getEnergyAmplitude(roi=roi)
        else:
            print 'Error: {} is not a valid domain'.format(domain)
            return None

        if y.sum() == 0:
            return np.nan, np.nan

        center = np.average(x, weights=y)
        width = np.sqrt( np.average((x-center)**2, weights=y-y.min()) )

        return center, width



def psanaTester():
    import psana
    import matplotlib.pyplot as plt
    import scipy.signal
    from ZmqSender import zmqSender
    import time


    # Load the config file
    import cookieBoxDefaultConfig as config
    # Create the sender, but only if zmq should be used
    if config.useZmq:
        sender = zmqSender()
    else:
        plt.ion()

    # Create the tofData object
    tof = tofData(config.basicTofConfig, quiet=False)
    tof.filterTimeDomainData(method='wavelet', numPoints=4)

    EMin = config.minE_eV
    EMax = config.maxE_eV
    threshold = 0.02
    minWidth_eV = 3
    yValuesBar = np.array((0,2,1,1,2,0,1,1,2,0))
    xValuesBar = np.zeros(np.size(yValuesBar))

    # Connect to data source
    print 'Connecting to data soutrce:', config.dataSource
    ds = psana.DataSource(config.dataSource)
    print 'done'
    for num, evt in enumerate(ds.events()):
        print 'event ', num
        if num >= config.nEvents:
            break

        if num is 0:
            #Initialize the scalse
            tof.setupScales(ds.env(), np.linspace(EMin, EMax,
                config.nEnergyBins))
            # Get the x scales
            t = tof.getTimeScale_us()
            E = tof.getEnergyScale_eV()
            rawE = tof.getRawEnergyScale()

            # setup local plotting
            f1 = plt.figure(1)
            f1.clf()
            a1 = f1.add_subplot(211)
            l11, l12, = a1.plot(t,t, t,t, 'r')
            a1.autoscale_view(scalex=False)
            a2 = f1.add_subplot(212)
            l21, l22, l23= a2.plot(rawE,rawE, E,E, 'r.', yValuesBar, yValuesBar,
                    'k')
            a2.set_xlim(EMin, EMax)
            f1.show()

        # Get the y data
        tof.setRawData(evt)
        tAmpRaw = tof.getTimeAmplitude()
        tAmpF = tof.getTimeAmplitudeFiltered()
        EAmp = tof.getEnergyAmplitude()

        t_timing = time.time()
        min, max, th = tof.getTraceBounds(threshold, minWidth_eV)
        print 'thresholding time:', time.time()-t_timing
        print 'min =', min, 'max =', max, 'threshold = ', th
        print 'bar y:', yValuesBar
        xValuesBar[:3] = min
        xValuesBar[3:7] = (min+max)/2
        xValuesBar[7:] = max


        # local plotting
        l11.set_ydata(tAmpRaw)
        l12.set_ydata(tAmpF)
        a1.relim()
        a1.autoscale_view(scalex=False)

        l21.set_ydata(tAmpF)
        l22.set_ydata(EAmp)
        l23.set_ydata(yValuesBar*th)
        l23.set_xdata(xValuesBar)
        a2.relim()
        a2.autoscale_view(scalex=False)
        #a2.xlim(EMin, EMax)
        f1.canvas.draw()

        # Remote plotting
        if config.useZmq:
            packet = []
            if num is 0:
                plot1 = linePlot((line(t, tAmpRaw),
                    line(t, tAmpF)))
                plot2 = linePlot((line(E, EAmp), line(xValuesBar, yValuesBar*th)))
            else:
                plot1 = linePlot((line(y=tAmpRaw), line(y=tAmpF)))
                plot2 = linePlot((line(y=EAmp), line(x=xValuesBar, y=yValuesBar*th)))

            packet.append(plot1)
            packet.append(plot2)
            sender.sendObject(packet)

        time.sleep(1)

    raw_input('Press enter to exit...')
    if config.useZmq:
        del sender




def nonPsanaTester():
    import matplotlib.pyplot as plt
    import scipy.signal
    import h5py


    # Load the config file
    config = \
            loadConfig('/reg/neh/home/alindahl/amoc8114/configFiles/configTofDataModuleTester.json')
    plt.ion()

    # Open an hdf5 file
    fileName = ('/reg/neh/home/alindahl/amoc8114/output/keepers' +
            '/amoc8114_run109_2014-6-16_1.hdf5')
    file = h5py.File(fileName, 'r')

    # Make references to some of the data in the h5 file
    tAx = file['tof_timeScale_us']
    tAmpVec = file['tof_timeAmplitude_V']


    # Create the tofData object
    tof = tofData(config, quiet=False)
    #tof.filterTimeDomainData(False)
    tof.filterTimeDomainData(method='average')
    #tof.filterTimeDomainData(method='wienerDeconv',
    #        SNR=np.fromfile('../h5Analysis/SNRrun109.npy'),
    #        response=np.fromfile('../h5Analysis/responseFunction108.npy'))

    EMin = config.tof_minE_eV
    EMax = config.tof_maxE_eV
    threshold = 0.02
    minWidth_eV = 3
    yValuesBar = np.array((0,2,1,1,2,0,1,1,2,0))
    xValuesBar = np.zeros(np.size(yValuesBar))

    for num, tAmpI in enumerate(tAmpVec):
        if num >= config.nEvents:
            break
        print 'event ', num

        if num is 0:
            #Initialize the scalse
            tof.setupScales(None,
                    np.linspace(EMin, EMax, config.tof_nEnergyBins),
                    tAx[:])
            # Get the x scales
            t = tof.getTimeScale_us()
            E = tof.getEnergyScale_eV()
            rawE = tof.getRawEnergyScale()

            # setup local plotting
            f1 = plt.figure(1)
            f1.clf()
            a1 = f1.add_subplot(211)
            a1.plot(t,t, label='raw signal')
            a1.plot(t, t, 'r', label='filtered signal')
            l11, l12 = a1.get_lines()
            a1.autoscale_view(scalex=False)
            a1.legend()
            a2 = f1.add_subplot(212)
            a2.plot(rawE,rawE, label='raw energy data')
            a2.plot(E,E, 'r.', label='energy rescaled data')
            a2.plot(yValuesBar, yValuesBar, 'k', label='peak finding results')
            l21, l22, l23 = a2.get_lines()
            a2.set_xlim(EMin, EMax)
            a2.legend()
            f1.show()

        # Get the y data
        tof.setRawData(timeAmplitude_V=tAmpI[:])
        tAmpRaw = tof.getTimeAmplitude()
        if tAmpRaw == None:
            continue
        tAmpF = tof.getTimeAmplitudeFiltered()
        EAmp = tof.getEnergyAmplitude()


        t_timing = time.time()
        min, max, th = tof.getTraceBounds(threshold, minWidth_eV)
        print 'thresholding time:', time.time()-t_timing
        print 'min =', min, 'max =', max, 'threshold = ', th
        print 'bar y:', yValuesBar
        xValuesBar[:3] = min
        xValuesBar[3:7] = (min+max)/2
        xValuesBar[7:] = max


        # local plotting
        l11.set_ydata(tAmpRaw)
        l12.set_ydata(tAmpF)
        a1.relim()
        a1.autoscale_view(scalex=False)

        l21.set_ydata(tAmpF / tAmpF.max() * EAmp.max())
        l22.set_ydata(EAmp)
        l23.set_ydata(yValuesBar*th)
        l23.set_xdata(xValuesBar)
        a2.relim()
        a2.autoscale_view(scalex=False)
        #a2.xlim(EMin, EMax)
        f1.canvas.draw()


        time.sleep(1)





if __name__ == '__main__':
    psanaTester()
    #nonPsanaTester()
