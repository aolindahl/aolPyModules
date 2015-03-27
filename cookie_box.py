import numpy as np
import tof
from configuration import loadConfiguration as loadConfig
from aolUtil import struct
import sys
import random
import lmfit
from BurningDetectors_V6 import projector
import simplepsana

# A bunch of methods to take care of the cookie box data

_source_dict = {}
def get_source(source_string):
    global _source_dict
    if source_tring not in _source_dict:
        _source_dict[source_string] = psana.Source(source_string)
    return _source_dict[source_string]


proj = projector()

def model_function(params, x, y=None, eps=None):
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

def initial_params(y_data=None):
    params = lmfit.Parameters()
    params.add('A', 10, min=0)
    params.add('beta', 2, min=-1, max=2)
    params.add('tilt', 0, min = -np.pi/2, max=np.pi/2)
    params.add('linear', 0.5, min=0)
    #params.add('tilt', np.pi, min=-2*np.pi, max=2*np.pi)
    #params.add('linear', 0.5, min=-0.5, max=1.5)

    if y_data!=None:
        params['A'].value = y_data[np.isfinite(y_data)].mean()
        tilt = initial_angle(y_data)
        #params['tilt'].value = tilt
        #params['tilt'].min = tilt - 2*np.pi
        #params['tilt'].max = tilt + 2*np.pi

    return params

phi_deg = np.arange(0, 360, 22.5)
phi_rad = phi_deg * np.pi / 180
_angles_deg = np.arange(0, 360, 22.5)
_angles_rad = _angles_deg * np.pi / 180
_sin_vect = np.sin( 2*_angles_rad )
_cos_vect = np.cos( 2*_angles_rad )

def initial_angle(y_data):
    vec = y_data.copy()
    vec[np.isnan(vec)] = 0
    A = y_data.dot(_sin_vect)
    B = y_data.dot(_cos_vect)
    phi = 0.5 * ( np.arctan(A/B)
            + np.pi * ( -1 * ((A<0) & (B<0)) + 1 * ((A>0) & (B<0)) ) )
    return phi



def slice_from_range(data, range):
    data_list = isinstance(data, list)
    range_list = isinstance(range[0], list)
    if data_list and range_list:
        if len(data) != len(range):
            return None
        return [ slice( d.searchsorted(np.min(r)), d.searchsorted(np.max(r)) )
                    for d, r in zip(data, range) ]

    if data_list and (not range_list):
        return [ slice( d.searchsorted(np.min(range)),
            d.searchsorted(np.max(range)))
            for d in data ]

    if (not data_list) and range_list:
        return [ slice( data.searchsorted(np.min(r)),
            data.searchsorted(np.max(r)))
            for r in range ]

    return slice( data.searchsorted(np.min(range)),
            data.searchsorted(np.max(range)))


def get_raw_signals(evt, source_string, time_slice=slice(None), verbose=False):
    if verbose:
        print 'Trying to grab raw signals from source{}.'.format(get_source(source_string))
    try:
        # try to get the acqiris data
        acqiris_data = evt.get(psana.Acqiris.DataDescV1, get_source(source_string))
    except:
        if verbose:
            print 'Fail. Exception was thrown.'
        return None
            
    if acqiris_data is None:
        return None

    return np.array([ acqiris_data.data(ch).waveforms()[0][time_slice] for ch in
        range(16) ])

#def sumSignals(signals, slices):
#    if type(slices) != slice:
#        slices = [slices]*16
#    [

def get_signal_scaling(env, source_string, verbose=False):
    temp = np.array( [tof.get_acqiris_scales(env, source_string, ch,
        verbose=verbose) for ch in range(16)] )
    return temp[:,1], temp[:,2]


class CookieBox:
    'Class that handels the cookiebox data'
    def __init__(self, config, verbose=False):
        '''\
        Initialization method.\
        
        The configuration should be a single object or a list of 16 objects.\
        '''
        if verbose:
            print 'In cookie_box.CookieBox.__init__().'

        self._verbose = verbose

        # A list of the TofData objects
        if verbose:
            print 'Make the list of TofData objects.'
        self._tof_list = []
        for conf in config.tof_config_list:
            self._tof_list.append(tof.TofData(conf, verbose=verbose))


        self._phi_deg = _angles_deg
        self._phi_rad = _angles_rad

        self._time_amplitudes_up_to_date = False

        self.proj = projector()

    def setup_scales(self, energy_scale_eV,
            env=None, time_scale=None, retardation=0, verbose=None):
        #print 'energy scale: ', energy_scale_eV
        if verbose is None:
            verbose = self._verbose
        if verbose:
            print 'In cookie_box.CookieBox.setup_scales().'
        for  tof in self._tof_list:
            tof.setup_scales(energy_scale_eV, env, time_scale,
                             retardation=retardation)


    def setBaselineSubtractionAveraging(self, weightLast):
        for tof in self._tof_list:
            tof.setBaselineSubtractionAveraging(weightLast)

    def set_raw_data(self, evt, verbose=False, newDataFactor=1):
        for tof in self._tof_list:
            tof.set_raw_data(evt, newDataFactor=newDataFactor)
        self._time_amplitudes_up_to_date = False

    def get_time_scales_us(self, roi=None, verbose=False):
        return [tof.get_time_scale_us(roi=roi) for tof in self._tof_list]

    def get_time_amplitudes(self, roiSlices=None, verbose=False):
        if not self._time_amplitudes_up_to_date:
            self._timeAmplitudes = [tof.get_time_amplitude() for tof in
                    self._tof_list]
        if roiSlices is None:
            return self._timeAmplitudes
        return [amp[s] for amp, s in zip(self._timeAmplitudes, roiSlices)]

    def get_time_amplitudes_filtered(self, roi=None, verbose=False):
        return [tof.get_time_amplitude_filtered(roi=roi) for tof 
                in self._tof_list]

    def get_energy_scales_eV(self, roi=None, verbose=False):
        return [tof.get_energy_scale_eV(roi=roi) for tof in self._tof_list]
    
    def get_energy_amplitudes(self, roi=None, verbose=False):
        return [tof.get_energy_amplitude(roi=roi) for tof in self._tof_list]

    def get_moments(self, domain='Time', roi=None):
        return [ tof.get_moments(domain=domain, roi=roi) for tof in
                self._tof_list]

    def get_positions(self):
        moments = np.array(self.get_moments(domain='Time', roi=0))
        positions = np.array([moments[i,0] - moments[i+8,0] for i in range(8)])
        return positions

    
    def get_photon_energy(self, energyShift=0):
        moments = np.array(self.get_moments(domain='Energy'))
        amplitudes = self.get_intensity_distribution(domain='Energy')
    
        return (np.average(moments, axis=0, weights=amplitudes)
                + np.array([energyShift, 0]))


    def get_intensity_distribution(self, rois=[slice(None)]*16,
            domain='Time', verbose=None, detFactors=[1]*16):
        if verbose is None:
            verbose = self._verbose
        if verbose:
            print 'Geting the intensity distribution',
            if rois is None:
                print '.'
            else:
                print 'in roi {}.'.format(rois)
        intensity = []
        if domain=='Energy':
            if verbose:
                print 'Using energy domain.'
            ampFunk = self.get_energy_amplitudes
        else:
            ampFunk = self.get_time_amplitudes_filtered
        for amp, factor, roi in zip(ampFunk(verbose=verbose,
                                            roi=(rois if isinstance(rois, int)
                                                else None)),
                                            detFactors,
                                            rois if isinstance(rois, list) else
                                                [rois]*16):
            if amp is None:
                intensity.append(np.nan)
            else:
                intensity.append(amp[roi if isinstance(roi, slice) else
                    slice(None)].sum() * factor)

        if verbose:
            print 'Returning vector of length {}.'.format(len(intensity))

        return np.array(intensity)

    def getAngles(self, kind='rad'):
        if kind=='rad':
            return self._phi_rad
        else:
            return self._phi_deg

    def randomize_amplitudes(self, verbose=False):
        '''This is only for testing. Direct manipulation of the private
        variables of the TofData class (as done here) is not recommended.'''
        params = initial_params()
        params['linear'].value = random.random()
        params['tilt'].value = random.random()*np.pi - np.pi/2
        params['A'].value = 1
        if verbose:
            for par in params.itervalues():
                print par
        factors = ( model_function(params, self.getAngles()) *
                np.random.normal(1, 0.05, (16,)) )
        #factors = model_function(params, self.getAngles())
        for factor, tof in zip(factors, self._tof_list):
            tof._timeAmplitude *= factor
            #tof._timeAmplitude_filtered *= factor
            tof._energyAmplitude *= factor

        return params
        



if __name__ == '__main__':
    do_plot = False
    verbose = True
    import time
    if do_plot:
        import matplotlib.pyplot as plt
        plt.ion()
    
    if verbose:
        print 'Connect to data source.'
    ds = simplepsana.get_data_source('exp=amoi0114:run=33', verbose=verbose)


    if verbose:
        print 'Import the configutration.'
    import cookie_box_default_config as config
    reload(config)

    config

    if hasattr(config, 'domainToDisplay') and config.domainToDisplay=='Energy':
        domain = 'Energy'
    else:
        domain = 'Time'

    if verbose:
        print 'Create the CookieBox object.'
    cb = CookieBox(config, verbose=verbose)

    t_temp = None
    if verbose:
        print 'Go through the events.'
    for i, evt in enumerate(ds.events()):
        if t_temp is not None:
            print 'event processing time:', time.time()-t_temp, 's'
        t_temp = time.time()
        if i >= 1:
            break
        if i%10==0 and do_plot:
            plot=True
        else:
            plot=False
        if i==0:
            if verbose:
                print 'Set up the time and energy scales.'
            cb.setup_scales(config.energy_scale_eV, ds.env())
            
            if verbose:
                print 'Get the energy scales.'
            energy_scales = cb.get_energy_scales_eV()
            #x_scales_roi_0 = cb.get_energy_scales_eV(roi=0)
            #x_scales_roi_1 = cb.get_energy_scales_eV(roi=1)
            if verbose:
                print 'Get the time scales.'
            time_scales = cb.get_time_scales_us()
            if verbose:
                print 'Get the ROI slices in time domain.'
            t_roi_0s = slice_from_range(time_scales, config.time_roi_0_us)
            t_roi_1s = slice_from_range(time_scales, config.time_roi_1_us)
            if verbose:
                print 'Get the time scales corresponding to the ROIs.'
            time_scales_roi_0 = [t[s] for t,s in zip(time_scales, t_roi_0s)]
            time_scales_roi_1 = [t[s] for t,s in zip(time_scales, t_roi_1s)]
            
            angles = cb.getAngles('rad')
            if do_plot:
                fig1 = plt.figure(1); plt.clf()
                for k, x , x_roi_0, x_roi_1 in zip(range(16), time_scales,
                        time_scales_roi_0, time_scales_roi_1):
                    plt.subplot(4,4,k+1)
                    plt.title('det at {} deg'.format(cb.getAngles('deg')[k]))
                    plt.plot(x, np.zeros_like(x))
                    plt.plot(x_roi_0, np.zeros_like(x_roi_0), 'r')
                    plt.plot(x_roi_1, np.zeros_like(x_roi_1), 'g')
                    yTMin, yTMax = 0, 0
    
                fig2 = plt.figure(2);
                fig2.clf()
                fig2ax = fig2.add_subplot(111, polar=True)
                angles_fit = np.linspace(0, 2*np.pi, 10000)
                fig2ax.plot(
                        angles, np.ones_like(angles), 'ro',
                        angles, np.ones_like(angles), 'gs',
                        angles_fit, np.zeros_like(angles_fit), 'm-')
                        #angles_fit, np.zeros_like(angles_fit), 'm--')

                fig3 = plt.figure(3)
                fig3.clf()
                fig3ax = fig3.add_subplot(111)
                fig3ax.plot(energy_scales[0],
                            np.zeros_like(energy_scales[0]))
        print 'event number', i


        if verbose:
            print 'Set the data of the event to the data structure.'
        cb.set_raw_data(evt)
        #rand_params = cb.randomize_amplitudes(verbose=True)
        #rand_tilt.append(rand_params['tilt'].value) 

        if plot:
            amplitudes0 = cb.get_intensity_distribution(roi=0, domain=domain,
                verbose=True)
            amplitudes1 = cb.get_intensity_distribution(roi=1, domain=domain,
                verbose=True)
            
            energy_data = cb.get_energy_amplitudes()
            #print len(energy_data)
            #print energy_data[0].shape
            #spectrum = np.average(energy_data, axis=0)
            spectrum = energy_data[15]
            #print spectrum
            #energy_data_roi_0 = cb.get_energy_amplitudes(roi=0)
            #energy_data_roi_1 = cb.get_energy_amplitudes(roi=1)
            time_data = cb.get_time_amplitudes()
            time_data_roi_0 = [t[s] for t,s in zip(time_data, t_roi_0s)]
            time_data_roi_1 = [t[s] for t,s in zip(time_data, t_roi_1s)]
            
            tMin = np.min(time_data)
            tMax = np.max(time_data)
            rescale_flag = False
            if tMin < yTMin:
                yTMin = tMin
                rescale_flag = True
            if tMax > yTMax:
                yTMax = tMax
                rescale_flag = True
        
            for ax, y, y_roi_0, y_roi_1 in zip(fig1.axes, time_data, time_data_roi_0,
                    time_data_roi_1):
                ax.lines[0].set_y_data(y)
                ax.lines[1].set_y_data(y_roi_0)
                ax.lines[2].set_y_data(y_roi_1)
                if rescale_flag:
                    ax.set_ybound(yTMin, yTMax)

        if verbose:
            print 'Get the signal amplitudes.'
        amplitudes = cb.get_intensity_distribution(rois=0, domain=domain,
                verbose=verbose)
        params = initial_params(amplitudes)
        #proj_tilt.append(params['tilt'].value)
        params['beta'].vary=False
        res = lmfit.minimize(model_function, params, args=(angles, amplitudes),
                method='leastsq')
        print res.nfev, 'function evaluations'
        print 'Fit', ('succeded' if res.success else 'failed')
        print res.message
        
        print lmfit.fit_report(params)


        moments = cb.get_moments(domain=domain, roi=0)

        if plot:
            fig2ax.lines[0].set_y_data(amplitudes0)
            fig2ax.lines[1].set_y_data(amplitudes1)
            fig2ax.lines[2].set_y_data(model_function(params, angles_fit))
            #fig2ax.lines[3].set_y_data(amplitudes0.mean() * model_function(rand_params, angles_fit))

            fig2ax.relim()
            fig2ax.autoscale_view()


            fig3ax.lines[0].set_y_data(spectrum)
            fig3ax.relim()
            fig3ax.autoscale_view()

            for fig in [fig1, fig2, fig3]:
                fig.canvas.draw()


    #rand_tilt = np.array(rand_tilt)
    #proj_tilt = np.array(proj_tilt)
    raw_input('Press enter to exit...')
