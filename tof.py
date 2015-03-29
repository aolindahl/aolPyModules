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
                     E_offset_eV=0, verbose=False):
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
    if verbose:
        print 'In tof.energy_from_time_physical()'
    return (D_mm**2 * m_e_eV * 1e6 /
        (c_0_mps**2 * 2 * (time - prompt_us - t_offset_us)**2) + E_offset_eV)



def get_time_to_energy_conversion(time_scale_us, energy_scale_eV, verbose=False,
                                  D_mm=None, prompt_us=None, t_offset_us=0,
                                  E_offset_eV=0):
    if verbose:
        print 'In "tof.get_energy_scale_and conversion()."'
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
    #print energy_scale_eV
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
    raw_energy_scale_eV[I] = energy_from_time_physical(time_scale_us[I],
                                                       D_mm=D_mm,
                                                       prompt_us=prompt_us,
                                                       t_offset_us=t_offset_us,
                                                       E_offset_eV=E_offset_eV)
    # Set the Energies that correspond to times befor the prompt to the
    # maximum energy.
    raw_energy_scale_eV[~I] = np.max(raw_energy_scale_eV)

    return conversion_mat, raw_energy_scale_eV


def get_acqiris_data(evt, source_string, channel, scaling=1., offset=0,
                     invert=True, selection=slice(None), verbose=False):
    raw_data = simplepsana.get_acqiris_waveform(evt, source_string,channel,
                                                verbose=verbose)
    if raw_data is None:
        return None

    invert_factor = -1. if invert else 1.
    return invert_factor * (raw_data[selection] * scaling - offset)


class TofData(object):
    """Class to handle TOF data.

    Conversions to energy scale including rescaling of the amplitudes\
    """

    def __init__(self, config, verbose=False):
        """\
        Initialize the TofData class giving it a configutaion object.\
        """

        # Extract the data source
        self._source_string = config['detectorSource']
        self._source = simplepsana.get_source(config['detectorSource'])
        self._acqiris_source_string = config['detectorSource']
        # Extract the acqiris channel
        self._acqiris_channel = config['acqCh']
        # Load the callibration file pointed to in the configuration
        if verbose:
            print 'Load the calibration from file "{}".'.format(
                    config['calib_file'])
        self._calibration= load_configuration_dict(config['calib_file'],
                                                   verbose=verbose)
        

        # Store the configuration
        self._config = config

        # Setup the basic rescaling
        self._acq_vert_scaling = 1
        self._acq_vert_offset = 0

        # Initialy the class does not contain any data or any scales
        self._no_data = True
        self._no_scales = True

        # Basic info about filtering
        self._filter_time_domain = False
        self._filter_method = None
        self._time_amplitude_filtered = None

        self._verbose = verbose

        if 'filterMethod' in self._config.keys():
            if self._config['filterMethod'] == "wienerDeconv":
                self.setup_time_domain_filtering(method='wienerDeconv',
                        SNR=self._config['filterWienerSNR'],
                        response=self._config['filterWienerResponse'])
            elif self._config['filterMethod'] == 'wavelet':
                self.setup_time_domain_filtering(method='wavelet',
                        levels=self._config["filterWaveletLevels"])
            elif self._config['filterMethod'] == 'average':
                if self._verbose:
                    print 'Using averaging with {} points'.format(
                            self._config['filterAverageNumPoints'])
                self.setup_time_domain_filtering(method='average',
                        numPoints=self._config['filterAverageNumPoints'])


        self._time_amplitude = None
        self._time_roi_slice = [None, None, None, None]
        self._energy_roi_slice =[None, None, None, None]

        self._bgWeight = None



    def setup_scales(self, energy_scale_eV, env=None,
                    time_scale_us=None, retardation=0):
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
        if ('t_slice' in self._config) and (self._config['t_slice'] == True):
            if 't_min_us' in self._config:
                slice_start = self._time_scale_us.searchsorted(
                        self._config['t_min_us'])
            else:
                slice_start = None

            if 't_max_us' in self._config:
                slice_end = self._time_scale_us.searchsorted(
                        self._config['t_max_us'])
            else:
                slice_end = None
            self._time_slice = slice(slice_start, slice_end)
        else:
            self._time_slice = slice(None)
        if self._verbose:
            print 'Time slice is: {}.'.format(self._time_slice)

        # Adjust the time scale
        self._time_scale_us = self._time_scale_us[self._time_slice]

        # Energy scale and conversion matrix
        self._calibration['E_offset_eV'] += retardation
        self._energy_scale_eV = energy_scale_eV
        #self._energy_scale_eV = energy_scale_eV[
        #        energy_scale_eV.searchsorted(self._calibration['E_offset_eV']):]
        self._energy_bin_size = np.diff(energy_scale_eV).mean()
        self._time_to_energy_matrix, self._raw_energy_scale_eV = \
            get_time_to_energy_conversion(self._time_scale_us,
                                          self._energy_scale_eV,
                                          verbose=self._verbose,
                                          **self._calibration)


        # Set the region of interest slices
        if self._verbose:
            print 'Looking for ROIs'
        for domain, roiBase in zip(
                ['Time', 'Energy'],
                ['time_roi_{}_us', 'energy_roi_{}_eV']):
            for iRoi in range(4):
                roi = roiBase.format(iRoi)
                roi_bg = roiBase.format('{}_bg'.format(iRoi))
                if self._verbose:
                    print 'Looking for {}.'.format(roi)
                if roi in self._config:
                    if self._verbose:
                        print '{} found'.format(roi)
                        print self._config[roi]
                    self.set_base_roi(
                            min = self._config[roi][0],
                            max = self._config[roi][1],
                            roi = iRoi,
                            domain = domain)
                if roi_bg in self._config:
                    pass


        # Make a backgeound slice
        if self._verbose:
            print 'Make background slice.'
        if self._config['baselineSubtraction'] == 'early':
            self._bgSlice = slice(self._time_scale_us.searchsorted(
                        self._config['baselineEnd_us']))
        elif self._config['baselineSubtraction'] == 'roi':
            try:
                self._bgSlice = slice(
                        min( [i.stop for i in self._time_roi_slice] ),
                        max( [i.start for i in self._time_roi_slice] ) )
            except:
                print "Could not set the gsSlice from the roi's."
                print "Attempted slice({}, {})".format(
                        min( [i.stop for i in self._time_roi_slice] ),
                        max( [i.start for i in self._time_roi_slice] ) )
                self._bgSlice = slice(0,0)
        else:
                self._bgSlice = slice(0,0)

        # Check if there is actually something to calculate the background
        # from
        if len(self._time_scale_us[self._bgSlice]) < 1:
            self._config['baselineSubtraction'] = 'none'
            if self._verbose:
                print 'Background slice is empty.'

        # There are scales, change the flag
        self._no_scales = False

    def set_base_roi(self, min, max, roi=0, domain='Time'):
        if domain=='Time':
            a = self._time_scale_us.searchsorted(min)
            b = a + self._time_scale_us[a:].searchsorted(max)
            self._time_roi_slice[roi] = slice(a,b)
        elif domain=='Energy':
            a = self._energy_scale_eV.searchsorted(min)
            b = a + self._energy_scale_eV[a:].searchsorted(max)
            self._energy_roi_slice[roi] = slice(a,b)

    def set_baseline_subtraction_averaging(self, weight_last):
        self._bg_weight_last = weight_last
        self._bg_weight_history = 1-weight_last
        self._bg_history = None


    def set_raw_data(self, evt=None, timeAmplitude_V=None, newDataFactor=None):
        """
        Set waveform data and compute the scaled data.
        """

        if self._verbose:
            print 'In tof.TofData.set_raw_data().'


        # If a psana event was passed as a parametere
        if evt is not None:
            if self._verbose:
                print 'Event object given.'

            new = get_acqiris_data(evt, self._source_string,
                                   self._acqiris_channel,
                                   scaling=self._acq_vert_scaling,
                                   offset=self._acq_vert_offset,
                                   invert=True, selection=self._time_slice,
                                   verbose=self._verbose)
            if new is None:
                if self._verbose:
                    print 'Could not extract any valide data.'
                self._no_data = True
                return


            if (self._time_amplitude is None or newDataFactor == 1 or
                    newDataFactor is None):
                if self._verbose:
                    print 'Using the new data.'
                self._time_amplitude = new
            else:
                if self._verbose:
                    print 'Updating rolling average.'
                self._time_amplitude = (self._time_amplitude * (1.0-newDataFactor)
                        + new * newDataFactor)

            # If set up to do that, subtract the baseline
            if self._config['baselineSubtraction'] != 'none':
                if self._verbose:
                    print 'Performing baseling subtraction.'
                if self._bgWeight is None:
                    self._time_amplitude -= \
                            self._time_amplitude[self._bgSlice].mean()
                else:
                    if self._bg_history is None:
                        self._bg_history =  \
                                self._time_amplitude[self._bgSlice].mean()
                    else:
                        self._bg_history *= self._bg_weight_history
                        self._bg_history += self._bg_weight_last \
                                * self._time_amplitude[self._bgSlice].mean()

                    self._time_amplitude -= self._bg_history

        elif timeAmplitude_V is  None:
            if self._verbose:
                print 'Niether event nor time scale given.'
            # If neither psana event, nor amplitude was given
            self._no_data = True
            return

        else:
            # If no psana event was given but threre is an amplitude
            self._time_amplitude = timeAmplitude_V[self._time_slice]

        if self._verbose:
            print 'Apply time domain filter.'
        self._filter_time_domain_data()

        if self._verbose:
            print 'Calculate energy amplitudes.'
        self.calc_energy_amplitude()

        self._no_data = False

    def calc_energy_amplitude(self, filtered=None):
        if self._filter_time_domain and filtered is not False:
            tAmp = self._time_amplitude_filtered
        else:
            tAmp = self._time_amplitude

        # Calculate the signal amplitude in the energy domain.
        self._energy_amplitude = self._time_to_energy_matrix.dot(tAmp)
        #self._energy_amplitude = np.ones(len(self._energy_scale_eV))
        return


    def get_time_scale_us(self, roi=None):
        if self._no_scales:
            return None
        if roi!=None and self._time_roi_slice!=None:
            return self._time_scale_us[self._time_roi_slice[roi]]
        return self._time_scale_us

    def get_time_amplitude(self, roi=None):
        #if self._verbose:
        #    print 'In tof.TofData.get_time_amplitude.'
        #    print 'Has', 'no' if self._no_data else None, 'data.'
        #    print self._time_amplitude
        if self._no_data:
            return None
        if (roi is not None) and (self._time_roi_slice is not None):
            return self._time_amplitude[self._time_roi_slice[roi]]
        return self._time_amplitude

    def get_time_amplitude_filtered(self, roi=None):
        if self._no_data:
            return None
        if roi!=None: #and self._time_roi_slice!=None:
            return self._time_amplitude_filtered[self._time_roi_slice[roi]]
        return self._time_amplitude_filtered

    def get_energy_amplitude(self, roi=None):
        if self._no_data:
            return None
        if roi!=None and self._energy_roi_slice!=None:
            return self._energy_amplitude[self._energy_roi_slice[roi]]
        return self._energy_amplitude

    def get_energy_scale_eV(self, roi=None):
        if self._no_scales:
            return None
        if roi!=None and self._energy_roi_slice!=None:
            return self._energy_scale_eV[self._energy_roi_slice[roi]]
        return self._energy_scale_eV

    def get_raw_energy_scale(self, roi=None):
        if self._no_scales:
            return None
        if roi!=None and self._time_roi_slice!=None:
            return self._energy_scale_eV[self._time_roi_slice[roi]]
        return self._raw_energy_scale_eV

    def setup_time_domain_filtering(self, method='wienerDeconv', numPoints=4, levels=6,
            SNR=1, response=None):
        if self._verbose:
            'In tof.TofData.setup_time_domain_filtering().'
        if method is False:
            self._filter_time_domain = False
            return

        self._filter_method = method
        self._filterNumPoints = numPoints
        self._filterLevels = levels
        self._filter_time_domain = True

        if method == 'wienerDeconv':
            if type(SNR) == str:
                self._SNR = np.fromfile(SNR)
            else:
                self._SNR = SNR
            if type(response) == str:
                self._response = np.fromfile(response)
            else:
                self._response = response


    def _filter_time_domain_data(self):
        if self._verbose:
            print ('In tof.TofData._filter_time_domain_data(),' +
                    ' with method = {}.'.format(self._filter_method))
        if self._filter_time_domain is False:
            self._time_amplitude_filtered = self._time_amplitude
            return

        if self._filter_method == 'average':
            self._time_amplitude_filtered = \
                    np.convolve(self._time_amplitude,
                            np.ones((self._filterNumPoints,))
                            / self._filterNumPoints, mode='same')
            return

        if self._filter_method == 'wavelet' and _useWavelet:
            #print 'Wavelet filteing'
            self._time_amplitude_filtered = wavelet(self._time_amplitude,
                    levels=self._filterLevels)
            return

        if self._filter_method == 'wienerDeconv':
            self._time_amplitude_filtered = wiener.deconvolution(
                    self._time_amplitude, self._SNR, self._response)
            return

        print '{} is not a valid filtering method when "_useWavelet" is {}.'\
                .format(self._filter_method, _useWavelet)


    def get_trace_bounds(self, threshold_V=0.02, min_width_eV=2, energy_offset=0,
                         useRel=False, threshold_rel=0.5):
        if self._no_data:
            return [np.nan for i in range(3)]

        if useRel:
            threshold_temp = threshold_rel * \
            np.max(self._energy_amplitude[np.isfinite(self._energy_amplitude)])
            if threshold_temp < threshold_V:
                return [np.nan for i in range(3)]
            else:
                threshold_V = threshold_temp
        nPoints = np.round(min_width_eV/self._energy_bin_size)

        min = 0
        for i in range(1, self._energy_amplitude.size):
            if self._energy_amplitude[i] < threshold_V:
                min = i
                continue
            if i-min >= nPoints:
                break
        else:
            min = np.nan


        max = self._energy_amplitude.size - 1
        for i in range(self._energy_amplitude.size-1, -1, -1):
            if self._energy_amplitude[i] < threshold_V:
                max = i
                continue
            if max-i >= nPoints:
                break
        else: max = np.nan

        if min == 0 and max == self._energy_amplitude.size - 1:
            min = np.nan
            max = np.nan

        if not np.isnan(min):
                min = self._energy_scale_eV[min] - energy_offset
                max = self._energy_scale_eV[max] - energy_offset

        return min, max, threshold_V



    def get_pulse_duration(self, lo, hi):
        if self._no_data:
            return None
        amplitude = self._config.tof_maxStreaking_eV
        cutoff = self._config.tof_streakingCutoff_eV
        if hi > cutoff or lo < -cutoff:
            return None
        dur = (np.arccos(lo/amplitude) - np.arccos(hi/amplitude)) / np.pi * \
            self._config.tof_quarterCycle_fs

        return dur

    def get_moments(self, domain='Time', roi=None):
        if domain == 'Time':
            x = self.get_time_scale_us(roi=roi)
            y = self.get_time_amplitude_filtered(roi=roi)
        elif domain == 'Energy':
            x = self.get_energy_scale_eV(roi=roi)
            y = self.get_energy_amplitude(roi=roi)
        else:
            print 'Error: {} is not a valid domain'.format(domain)
            return None

        if y.sum() == 0:
            return np.nan, np.nan

        center = np.average(x, weights=y)
        width = np.sqrt( np.average((x-center)**2, weights=y-y.min()) )

        return center, width



def psanaTester(do_plot=False):
    import scipy.signal
    import time
    if do_plot:
        from ZmqSender import zmqSender
        import matplotlib.pyplot as plt


    # Load the config file
    import cookiebox_default_config as config

    # Create the sender, but only if zmq should be used
    if do_plot:
        if config.useZmq:
            sender = zmqSender()
        else:
            plt.ion()

    # Create the TofData object
    tof = TofData(config.basic_tof_config, verbose=True)
    tof.setup_time_domain_filtering(method='wavelet', numPoints=4)

    EMin = config.minE_eV
    EMax = config.maxE_eV
    threshold = 0.02
    min_width_eV = 3
    yValuesBar = np.array((0,2,1,1,2,0,1,1,2,0))
    xValuesBar = np.zeros(np.size(yValuesBar))

    # Connect to data source
    print 'Connecting to data soutrce:', config.dataSource
    ds = simplepsana.get_data_source(config.offline_source)
    print 'done'
    for num, evt in enumerate(ds.events()):
        print 'event ', num
        if num >= config.nEvents:
            break

        if num is 0:
            #Initialize the scalse
            print '\nInitialize the scales in the tof object.'
            tof.setup_scales(np.linspace(EMin, EMax, config.nEnergyBins),
                            env=ds.env())
            # Get the x scales
            print '\nGet the time scale.'
            t = tof.get_time_scale_us()
            print '\nGet the energy scale.'
            E = tof.get_energy_scale_eV()
            print '\nGet the raw energy scale.'
            rawE = tof.get_raw_energy_scale()

            # setup local plotting
            if do_plot:
                print 'Initialize plotting.'
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
        print '\nSet raw data ov event to tof object.'
        tof.set_raw_data(evt)
        print '\nGet some vectors.'
        tAmpRaw = tof.get_time_amplitude()
        tAmpF = tof.get_time_amplitude_filtered()
        EAmp = tof.get_energy_amplitude()

        t_timing = time.time()
        min, max, th = tof.get_trace_bounds(threshold, min_width_eV)
        print 'thresholding time:', time.time()-t_timing
        print 'min =', min, 'max =', max, 'threshold = ', th
        print 'bar y:', yValuesBar
        xValuesBar[:3] = min
        xValuesBar[3:7] = (min+max)/2
        xValuesBar[7:] = max


        if do_plot:
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
    if config.useZmq and do_plot:
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


    # Create the TofData object
    tof = TofData(config, verbose=True)
    #tof.setup_time_domain_filtering(False)
    tof.setup_time_domain_filtering(method='average')
    #tof.setup_time_domain_filtering(method='wienerDeconv',
    #        SNR=np.fromfile('../h5Analysis/SNRrun109.npy'),
    #        response=np.fromfile('../h5Analysis/responseFunction108.npy'))

    EMin = config.tof_minE_eV
    EMax = config.tof_maxE_eV
    threshold = 0.02
    min_width_eV = 3
    yValuesBar = np.array((0,2,1,1,2,0,1,1,2,0))
    xValuesBar = np.zeros(np.size(yValuesBar))

    for num, tAmpI in enumerate(tAmpVec):
        if num >= config.nEvents:
            break
        print 'event ', num

        if num is 0:
            #Initialize the scalse
            tof.setup_scales(None,
                    np.linspace(EMin, EMax, config.tof_nEnergyBins),
                    tAx[:])
            # Get the x scales
            t = tof.get_time_scale_us()
            E = tof.get_energy_scale_eV()
            rawE = tof.get_raw_energy_scale()

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
        tof.set_raw_data(timeAmplitude_V=tAmpI[:])
        tAmpRaw = tof.get_time_amplitude()
        if tAmpRaw == None:
            continue
        tAmpF = tof.get_time_amplitude_filtered()
        EAmp = tof.get_energy_amplitude()


        t_timing = time.time()
        min, max, th = tof.get_trace_bounds(threshold, min_width_eV)
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
