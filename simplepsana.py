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
    raise ImportError('\n\t'.join([
        'Import of module "{}" failed due to an exception:'.format(__name__),
        '"{}"'.format(exc.message),
        'when attempting to "import psana".'
        ]))


def get_data_source(source_string):
    return psana.DataSource(source_string)


# Internal dictionary in the module to keep track of the used data sources
_sourceDict = {}
def get_source(source_string):
    if psana is None:
        print 'ERROR: Function "getSoutrce" cannot be used without psana.'
        sys.exit()
    global _sourceDict
    if source_string not in _sourceDict:
        _sourceDict[source_string] = psana.Source(source_string)
    return _sourceDict[source_string]



def get_acqiris_time_scale_us(env, source_string, verbose=False):
    """Get the time scale of the aquiris trace from an environment object.

    Returns None at failiure.
    Unit is microseconds."""

    if psana is None:
        print 'ERROR: Function "getTimeScales_us" cannot be used without psana.'
        sys.exit()

    # Get the configuration
    try:
        acqiris_config = env.configStore().get(psana.Acqiris.ConfigV1,
                get_source(source_string) )
    except:
        raise
        return None

    # make the time scale vector for the acqiris channel.
    # This is just for convenience
    timeScale = acqiris_config.horiz()
    # Start time
    t0 = timeScale.delayTime()
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
        print 'Is get_acqiris_signal_scaling().'
    # Get the configuration
    try:
        acqiris_config = env.configStore().get(psana.Acqiris.ConfigV1,
                get_source(source_string) )
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
    # all the 2**16 bits.
    # Here the voltage per bit is calculated
    scaling = vert_scale.fullScale() * 2**-16
    # The scale also has an offset in voltage
    offset = vert_scale.offset()

    return scaling, offset



if __name__ == '__main__':
    import sys
    print 'Running module tests'

    # Connect to a data source
    if len(sys.argv) < 2:
        source_string = 'exp=amoc8114:run=24'
    else:
        soiurceString = sys.argv[1]
    print 'Connecting to data source "{}".'.format(source_string)
    ds = psana.DataSource(source_string)
    print '\tDone'

    # Set up an Acqiris based detector
    acq_source_string = 'DetInfo(AmoETOF.0:Acqiris.0)'
    acq_channel = 1
    # Get the time scale
    time_scale_us = get_acqiris_time_scale_us(ds.env(), acq_source_string,
                                              verbose=True)
    scaling, offset = get_acqiris_signal_scaling(ds.env(), acq_source_string,
                                                acq_channel, verbose=True)
