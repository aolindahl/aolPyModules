# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 13:42:23 2015

@author: Anton O. Lindahl
"""
import sys
import numpy as np
try:
    import psana
except ImportError as exc:
    print ''.join([
        'Import of module "{}" failed due to an exception:'.format(__name__),
        '"{}"'.format(exc.message),
        'when attempting to "import psana".'])
#    raise ImportWarning('\n\t'.join([
#        'Import of module "{}" failed due to an exception:'.format(__name__),
#        '"{}"'.format(exc.message),
#        'when attempting to "import psana".'
#        ]))

def allow_corrupt_epics():
    psana.setOption('psana.allow-corrupt-epics', True)

def get_data_source(source_string, verbose=False):
    if verbose:
        print 'In "simplepsana.get_data_source()" function.'
    return psana.DataSource(source_string)


# Internal dictionary in the module to keep track of the used data sources
_sourceDict = {}
def get_source(source_string, verbose=False):
    if verbose:
        print 'In "simplepsana.get_source()" function.'
    if psana is None:
        print 'ERROR: Function "simplepsana.getSoutrce" cannot be used without psana.'
        sys.exit()
    global _sourceDict
    if source_string not in _sourceDict:
        _sourceDict[source_string] = psana.Source(source_string)
    return _sourceDict[source_string]



def get_acqiris_time_scale_us(env, source_string, verbose=False):
    """Get the time scale of the aquiris trace from an environment object.

    Returns None at failiure.
    Unit is microseconds."""
    if verbose:
        print 'In "simplepsana.get_acqiris_time_scale_us()" function.'
    if psana is None:
        print 'ERROR: Function implepsana.getTimeScales_us" cannot be used without psana.'
        sys.exit()

    # Get the configuration
    try:
        acqiris_config = env.configStore().get(psana.Acqiris.ConfigV1,
                get_source(source_string, verbose=verbose))
    except:
        raise
        return None

    # make the time scale vector for the acqiris channel.
    # This is just for convenience
    timeScale = acqiris_config.horiz()
    # Start time
    #t0 = timeScale.delayTime()
    t0 = 0
    # Time step
    dt = timeScale.sampInterval()
    # Number of samples
    nSample = timeScale.nbrSamples()
    # Make the time scale vector from the above information and rescale it
    # to microseconds
    return np.arange(t0, dt*nSample, dt)*1e6

def get_acqiris_signal_scaling(env, source_string, channel, verbose=False):
    """Get information on how to rescale the raw signal.

    env -- psana environment.
    source-string
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
    if verbose:
        print 'In "simplepsana.get_acqiris_signal_scaling()".'
    # Get the configuration
    try:
        acqiris_config = env.configStore().get(psana.Acqiris.ConfigV1,
                get_source(source_string, verbose=verbose) )
        if verbose:
            print 'Found acqiris configuration.'
    except:
        if verbose:
            print 'Did not find any configuration, returning (1, 0).'
        return 1, 0

    # Get the scaling constants for the vertical scale.
    # convenience reference
    vert_scale = acqiris_config.vert()[channel]
    # The vertical scale information is given as the full scale voltage over
    # all the 2**12 bits but they need to be shifted.
    # Here the voltage per bit is calculated
    scaling = vert_scale.fullScale() * 2**-12
    # The scale also has an offset in voltage
    offset = vert_scale.offset()

    return scaling, offset


def get_acqiris_waveform(evt, source_string, channel, segment=0, verbose=False):
    if verbose:
        print 'In simplepsana.get_acqiris_waveform()'
        print '\tGet the acqiris data structure from the event.'

    try:
        acqiris_data = evt.get(psana.Acqiris.DataDescV1,
                               get_source(source_string, verbose=verbose))
    except:
        acqiris_data = None

    if acqiris_data is None:
        if verbose:
            print '\tNo Acqiris data object in event.'
        return None

    channels = acqiris_data.data_shape()[0]
    if verbose:
        print 'Channel {} of {} (max {}) requested.'.format(channel, channels,
                channels-1)

    if channels <= channel:
        raise IndexError('Requested Acqiris channel = {} out'.format(channel) + 
                         'of bounds, max = {}.'.format(channels-1))

    if verbose:
        print '\tReturn waveform.'
    return acqiris_data.data(channel).waveforms()[segment] / 16

if __name__ == '__main__':
    import sys
    print 'Running module tests'
    verbose = True

    # Connect to a data source
    if len(sys.argv) < 2:
        source_string = 'exp=amoc8114:run=24'
    else:
        soiurceString = sys.argv[1]
    print '\nConnecting to data source "{}".'.format(source_string)
    ds = get_data_source(source_string, verbose=verbose)
    print '\tDone'

    # Set up an Acqiris based detector
    acq_source_string = 'DetInfo(AmoETOF.0:Acqiris.0)'
    acq_channel = 1
    # Get the time scale
    if verbose:
        print '\nTrying to get acqiris time scale.'
    time_scale_us = get_acqiris_time_scale_us(ds.env(), acq_source_string,
                                              verbose=True)
    if verbose:
        print '\nTrying to get acqiris vertical scaling.'
    scaling, offset = get_acqiris_signal_scaling(ds.env(), acq_source_string,
                                                acq_channel, verbose=True)
