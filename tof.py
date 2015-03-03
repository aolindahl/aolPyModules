#from setupEnvironment import *
import numpy as np
from configuration import loadConfiguration, load_configuration_dict
import time
import wiener
from scipy.sparse import coo_matrix
import sys

import simplepsana
import aolUtil

_useWavelet = True
if _useWavelet:
    try:
        from wavelet_filter import wavelet_filt as wavelet
    except:
        print 'Wavelet filtering not avaliable. pywt not found.'
        _useWavelet = False

m_e_eV = 0.510998928e6 # http://physics.nist.gov/cgi-bin/cuu/Value?me 2014-04-21
c_0_mps = 299792458 # http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=universal_in! 2014-04-21


def edges_from_centers(centers):
    step = np.diff(centers).mean()
    return np.concatenate([centers - step/s, centers[-1] + step/2])


def get_acqiris_scales(env, source_string, channel, verbose=False):
    # If a time scale is given, use it, oterhwhise try to get is from the
    # env object.
    if verbose:
        print 'Get the time scale of the acqiris using:'
        print '\tenv =', env
        print '\tsource_string =', source_string
    time_scale_us = simplepsana.get_acqiris_time_scale_us(env, source_string,
                                                          verbose = verbose)
    if time_scale_us is None:
        print ('WARNING: No acqiris configuration obtained,' +
                ' no scales aquired.')
        return None, 1, 0
    
    if verbose:
        print 'Get the vertical scaling of the acqiris channel.'
    vert_scaling_V, vert_offset_V = \
        simplepsana.get_acqiris_signal_scaling(env, source_string, channel,
                                               verbose = verbose)

    return time_scale_us, vert_scaling_V, vert_offset_V


def energy_from_time_physical(time, D_mm=None, prompt_us=None, t_offset_us=0,
                     E_offset_eV=0):
    # The conversion is given by:
    # E = D^2 m_e 10^6 / ( 2 c_0^2 (t - t_p)^2 ) + E0,
    # where:
    # D is the flight distance in mm
    # m_e is the electon rest mass expressed in eV
    # the 10^6 factor accounts the the otherwhise
    # missmachching prefixes
    # c_0 is the speed of light in m/s
    # E is the energy in eV
    # E0 is an energy offset in eV, should be
    # determined in a callibration fit
    # t_p is the arrival time of the prompt signal in microseconds
    return (D_mm**2 * m_e_eV * 1e6 /
        (c_0_mps**2 * 2 * (time - prompt_us - t_offset_us)**2) + E_offset_eV)



def get_time_to_energy_conversion(time_scale_us, energy_scale_eV, verbose=False,
                                  D_mm=None, prompt_us=None, t_offset_us=0,
                                  E_offset_eV=0):
    if verbose:
        print 'In get_energy_scale_and conversion()'
    # Get basic data about the time scale
    dt = np.diff(time_scale_us).mean()
    t0 = time_scale_us[0]
    num_time_bins = len(time_scale_us)
    # Calculate bin edges in the time domain
    time_scale_t_edges = aolUtil.limits_from_centers(time_scale_us)
    # and the corresponding bin edges in the energy domain.
    time_scale_E_edges = energy_from_time_physical(time_scale_t_edges,
                                                   D_mm=D_mm,
                                                   prompt_us=prompt_us,
                                                   t_offset_us=t_offset_us,
                                                   E_offset_eV=E_offset_eV)
    if verbose:
        print 'Time scale E edges are:', time_scale_E_edges
    # Number of energy bins
    num_energy_bins = len(energy_scale_eV)
    # Energy bin size
    dE = np.diff(energy_scale_eV).mean()
    # energy scale bin limits
    energy_scale_E_edges = aolUtil.limits_from_centers(energy_scale_eV)
    
    # Make matrixes out of the edges vectors in energy domain
    mat_time_E = np.concatenate([time_scale_E_edges.reshape(1,-1)] *
        num_energy_bins)
    mat_energy_E= np.concatenate([energy_scale_E_edges.reshape(-1,1)] *
        num_time_bins, axis=1)

    # Compute the start and end energies for the conversion from the time axis
    # to energy axis
    high_E_limit = ( np.minimum( mat_time_E[:,:-1], mat_energy_E[1:,:] ) )
    low_E_limit = ( np.maximum( mat_time_E[:,1:], mat_energy_E[:-1,:] ) )
    # Only where the high energy is more than the low energy the conversion
    # makes anny sense
    I = low_E_limit < high_E_limit
    # Allocate a tempoarary conversion matrix
    temp_conversion_mat = np.zeros((num_energy_bins, num_time_bins))
    # Calculate the elements of the conversion matrix
    # For each time bin the measured amplitude is multiplied by the bin size
    # in order to arrive at the integral of the signal. Then it is
    # determined how much of each time bin contributes to each energy bin.
    # This is done by comparing the edge positions of the time and energy
    # bins and assigning the correct proportion of the integral in time
    # domain to integral in the energy domain. Finally the total integral is
    # divided by the energy bin size in order to return to an amplitude.
    # Summation over all time bins is performed in the matrix multiplication
    # of the conversion matrix with the time domain amplitude vector.
    temp_conversion_mat[I] = (dt * (high_E_limit[I] - low_E_limit[I]) / 
            ( mat_time_E[:,:-1] - mat_time_E[:,1:] )[I] / dE)
    # The conversion matrix is highly sparse, thus make a sparse matrix to
    # spped up the calculationss
    conversion_mat = coo_matrix(temp_conversion_mat)

    # Create the energy scale
    # Find out which part of the time scale is after the prompt peak
    # The calibration information is used for this
    I = time_scale_us > prompt_us + t_offset_us
    # Allocate an energy scale with -1 value (unphysical).
    raw_energy_scale_eV = -np.ones_like(time_scale_us)
    raw_energy_scale_eV[I] = energy_from_time_physical(time_scale_us,
                                                       D_mm=D_mm,
                                                       prompt_us=prompt_us,
                                                       t_offset_us=t_offset_us,
                                                       E_offset_eV=E_offset_eV)
    # Set the Energies that correspond to times befor the prompt to the
    # maximum energy.
    raw_energy_scale_eV[~I] = np.max(self._rawEnergyScale_eV)

    return conversion_mat, raw_energy_scale_eV

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

    def __init__(self, config, verbose=False):
        """\
        Initialize the tofData class giving it a configutaion object.\
        """

        # Extract the data source
        self._source = simplepsana.get_source(config['detectorSource'])
        self._acqiris_source_string = config['detectorSource']
        # Extract the acqiris channel
        self._acqiris_channel = config['acqCh']
        # Load the callibration file pointed to in the configuration
        self._calibration= load_configuration_dict(config['calibFile'],
                                                   verbose=verbose)
        

        # Store the configuration
        self._config = config

        # Setup the basic rescaling
        self._acqVertScaling = 1
        self._acqVertOffset = 0

        # Initialy the class does not contain any data or any scales
        self._noData = True
        self._no_scales = True

        # Basic info about filtering
        self._filterTimeDomain = False
        self._timeAmplitudeFiltered = None

        self._verbose = verbose

        if 'filterMethod' in self._config.keys():
            if self._config['filterMethod'] == "wienerDeconv":
                self.filterTimeDomainData(method='wienerDeconv',
                        SNR=self._config['filterWienerSNR'],
                        response=self._config['filterWienerResponse'])
            elif self._config['filterMethod'] == 'wavelet':
                self.filterTimeDomainData(method='wavelet',
                        levels=self._config["filterWaveletLevels"])
            elif self._config['filterMethod'] == 'average':
                if self._verbose:
                    print 'Using averaging with {} points'.format(
                            self._config['filterAverageNumPoints'])
                self.filterTimeDomainData(method='average',
                        numPoints=self._config['filterAverageNumPoints'])


        self._timeAmplitude = None
        self._timeRoiSlice = [None, None]
        self._energyRoiSlice =[None, None]

        self._bgWeight = None



    def setupScales(self, energy_scale_eV, env=None,
                    time_scale_us=None, t_min_us=None, t_max_us=None):
        """Setup the information about the acqiris channel used for the TOF.

        Reads the scale factors for the raw aquiris data and also calculates
        the conversion to the energy domain.
        """

        if self._verbose:
            print 'Seting up the scales.'

        # Time scale
        if time_scale_us is None:
            self._time_scale_us, self._acq_vert_scaling, self._acq_vert_offset = \
                get_acqiris_scales(env, self._acqiris_source_string,
                                   self._acqiris_channel,
                                   verbose=self._verbose)
        else:
            self._time_scale_us = time_scale_us.copy()

        if self._time_scale_us is None:
            if self._verbose:
                print 'No scales found.'
            self._no_scales = True
            return

        # Time sliceing
        if self._verbose:
            print 'Seting up the time slicing.'
        if t_min_us is None:
            slice_start = None
        else:
            slice_start = self._time_scale_us.searchsorted(t_min_us)
        if t_max_us is None:
            slice_end = None
        else:
            slice_end = self._time_scale_us.searchsorted(t_max_us)
        self._time_slice = slice(slice_start, slice_end)

        # Energy scale and conversion matrix
        self._energy_scale_eV = energy_scale_eV
        self._time_to_energy_matrix, self._raw_energy_scale_eV = \
            get_time_to_energy_conversion(self._time_scale_us[self._time_slice],
                                          self._energy_scale_eV,
                                          verbose=self._verbose,
                                          **self._calibration)


        # Set the region of interest slices
        if self._verbose:
            print 'Looking for ROIs'
        for domain, roiBase in zip(
                ['Time', 'Energy'],
                ['timeRoi{}_us', 'energyRoi{}_eV']):
            for iRoi in range(2):
                roi = roiBase.format(iRoi)
                if roi in self._config:
                    if self._verbose:
                        print '\t{} found'.format(roi)
                        print self._config[roi]
                    self.setBaseRoi(
                            min = self._config[roi][0],
                            max = self._config[roi][1],
                            roi = iRoi,
                            domain = domain)


        # Make a backgeound slice
        if self._config['baselineSubtraction'] == 'early':
            self._bgSlice = slice(self._time_scale_us.searchsorted(
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
        if len(self._time_scale_us[self._bgSlice]) < 1:
            self._config['baselineSubtraction'] = 'none'



        self._no_scales = False

    def setBaseRoi(self, min, max, roi=0, domain='Time'):
        if domain=='Time':
            a = self._time_scale_us.searchsorted(min)
            b = a + self._time_scale_us[a:].searchsorted(max)
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

        if self._verbose:
            print 'Seting raw data.'


        # If a psan event was passed as a parametere
        if evt is not None:
            if self._verbose:
                print 'Event object given.'
            try:
                # try to get the acqiris data
                acqirisData = evt.get(psana.Acqiris.DataDescV1, self._source)
            except:
                self._noData = True
                return

            if acqirisData == None or self._no_scales:
                self._noData = True
                return

            # If everything worked out calculate the rescaled amplitude
            new = -(acqirisData.data(self._acqiris_channel).waveforms()[0][
                        self._time_slice]
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
            if self._verbose:
                print 'Niether event nor time scale given.'
            # If neither psana event, nor amplitude was given
            self._noData = True
            return

        else:
            # If no psana event was given but threre is an amplitude
            self._timeAmplitude = timeAmplitude_V[self._time_slice]

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
        if self._no_scales:
            return None
        if roi!=None and self._timeRoiSlice!=None:
            return self._time_scale_us[self._timeRoiSlice[roi]]
        return self._time_scale_us

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
        if self._no_scales:
            return None
        if roi!=None and self._energyRoiSlice!=None:
            return self._energyScale_eV[self._energyRoiSlice[roi]]
        return self._energyScale_eV

    def getRawEnergyScale(self, roi=None):
        if self._no_scales:
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
    import matplotlib.pyplot as plt
    import scipy.signal
    from ZmqSender import zmqSender
    import time


    # Load the config file
    import cookiebox_default_config as config
    # Create the sender, but only if zmq should be used
    if config.useZmq:
        sender = zmqSender()
    else:
        plt.ion()

    # Create the tofData object
    tof = tofData(config.basic_tof_config, verbose=True)
    tof.filterTimeDomainData(method='wavelet', numPoints=4)

    EMin = config.minE_eV
    EMax = config.maxE_eV
    threshold = 0.02
    minWidth_eV = 3
    yValuesBar = np.array((0,2,1,1,2,0,1,1,2,0))
    xValuesBar = np.zeros(np.size(yValuesBar))

    # Connect to data source
    print 'Connecting to data soutrce:', config.dataSource
    ds = simplepsana.get_data_source(config.dataSource)
    print 'done'
    for num, evt in enumerate(ds.events()):
        print 'event ', num
        if num >= config.nEvents:
            break

        if num is 0:
            #Initialize the scalse
            tof.setupScales(np.linspace(EMin, EMax, config.nEnergyBins),
                            ds.env())
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
    tAx = file['tof_time_scale_us']
    tAmpVec = file['tof_timeAmplitude_V']


    # Create the tofData object
    tof = tofData(config, verbose=True)
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
