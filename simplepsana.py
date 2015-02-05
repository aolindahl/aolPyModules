# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 13:42:23 2015

@author: Anton O. Lindahl
"""
import sys
try:
    import psana
except ImportError as exc:
    raise ImportError('\n\t'.join([
        'Import of module "{}" failed due to an exception:'.format(__name__),
        '"{}"'.format(exc.message),
        'when attempting to "import psana".'
        ]))


# Internal dictionary in the module to keep track of the used data sources
_sourceDict = {}
def getSource(sourceString):
    if psana is None:
        print 'ERROR: Function "getSoutrce" cannot be used without psana.'
        sys.exit()
    global _sourceDict
    if sourceString not in _sourceDict:
        _sourceDict[sourceString] = psana.Source(sourceString)
    return _sourceDict[sourceString]



def getTimeScale_us(env, sourceString, verbose=False):
    """Returnes time scale of aquiris trace from environment object

    Returns None at failiure.
    Unit is microseconds."""

    if psana is None:
        print 'ERROR: Function "getTimeScales_us" cannot be used without psana.'
        sys.exit()

    # Get the configuration
    try:
        acqirisConfig = env.configStore().get(psana.Acqiris.ConfigV1,
                getSource(sourceString) )
    except:
        return None

    # make the time scale vector for the acqiris channel.
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
    return np.arange(t0, dt*nSample, dt)*1e6


if __name__ == '__main__':
    print 'Running module tests'