import psana
import numpy as np


EBeamTypeList = (psana.Bld.BldDataEBeamV0,
        psana.Bld.BldDataEBeamV1,
        psana.Bld.BldDataEBeamV2,
        psana.Bld.BldDataEBeamV3,
        psana.Bld.BldDataEBeamV4,
        psana.Bld.BldDataEBeamV5,
        psana.Bld.BldDataEBeamV6)

_EBeamType = None
_EBeamSource = psana.Source('BldInfo(EBeam)')

feeTypeList = (psana.Bld.BldDataFEEGasDetEnergy,
        psana.Bld.BldDataFEEGasDetEnergyV1)

_feeType = None
_feeSource = psana.Source('BldInfo(FEEGasDetEnergy)')

_evt = None

def setEvent(evt, verbose=False):
    global _evt
    if verbose:
        print 'Updating the internal event.'
    _evt = evt

def _checkEvent(evt, verbose=False):
    if evt is not None:
        if verbose:
            print 'Event given.'
        setEvent(evt, verbose)
        return True
    if _evt is None:
        if verbose:
            print 'Internal event empty.'
        return False
    return True

# e-beam data
def getEBeamEnergyL3_MeV(evt=None, verbose=False):
    if not _checkEvent(evt, verbose):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamL3Energy()

def getEBeamEnergyBC2_MeV(evt=None, verbose=False):
    if not _checkEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    # The BC2 energy is calculated using the dispersion. The vispersion
    # value for BC2 is -3.647 mm at the nominal beam energy of 5 GeV [Email
    # from Timothy Maxwell to Anton Lindahl on June 2 2014]
    return ( EBeamObject.ebeamEnergyBC2() / -364.7 + 1 ) * 5e3

def getEBeamCharge_nC(evt=None, verbose=False):
    if not _checkEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamCharge()

def getEBeamPkCurrentBC2_A(evt=None, verbose=False):
    if not _checkEvent(evt):
        return np.nan
    EBeamObject = getEBeamObject(verbose=verbose)
    if EBeamObject == None:
        return np.nan
    return EBeamObject.ebeamPkCurrBC2()

# e-beam setup
def getEBeamObject(evt=None, verbose=False):
    if not _checkEvent(evt, verbose):
        return None
    # Initialize the EBeam type
    if _EBeamType is None:
        if verbose:
            print 'Initializing the EBeam type.'
        _determineEBeamType(verbose=verbose)
    return _evt.get(_EBeamType, _EBeamSource)

def _determineEBeamType(evt=None, verbose=False):
    global _EBeamType
    if not _checkEvent(evt, verbose):
        return None
    if verbose:
        print 'Find the correct EBeam type.'
    for type in EBeamTypeList:
        if verbose:
            print 'Trying {};'.format(type),
        data = _evt.get(type, _EBeamSource)
        if data is not None:
            _EBeamType = type
            if verbose:
                print ' correct.'
            break
        elif verbose:
            print ' wrong one.'

    return _EBeamType

# fee
def getPulseEnergy_mJ(evt=None, verbose=False):
    if not _checkEvent(evt):
        return np.nan
    fee = getFeeObject(_evt, verbose)
    if fee is None:
        return np.array([np.nan for i in range(4)])
    return np.array( [ fee.f_11_ENRC(),
            fee.f_12_ENRC(),
            fee.f_21_ENRC(),
            fee.f_22_ENRC()] )

def getFeeObject(evt, verbose=False):
    if _feeType is None:
        _determineFeeType(evt, verbose)
    return evt.get(_feeType, _feeSource)

def _determineFeeType(evt, verbose=False):
    global _feeType
    if verbose:
        print 'Find the correct fee type.'
    for type in feeTypeList:
        if verbose:
            print 'Trying {};'.format(type),
        data = evt.get(type, _feeSource)
        if data is not None:
            _feeType = type
            if verbose:
                print ' correct.'
            break
        elif verbose:
            print ' wrong one.'


# Event id
def getIdObject(evt, verbose=False):
    return evt.get(psana.EventId)

def getEventFiducial(evt=None, verbose=False):
    if not _checkEvent(evt, verbose):
        return 0
    id = getIdObject(_evt, verbose)
    if id == None:
        return 0
    return id.fiducials()

def getEventTime(evt=None, verbose=False):
    if not _checkEvent(evt, verbose):
        return np.nan
    id = getIdObject(_evt, verbose)
    if id == None:
        return np.nan
    time =  id.time()
    return time[0] + time[1] * 1e-9


if __name__ == '__main__':
    print 'Connecting to data source.'
    ds = psana.DataSource('exp=amoc8114:run=24')
    print 'Geting an event.'
    evt = ds.events().next()
    print 'Set the event'
    setEvent(evt, verbose=True)
    print 'E at L3 is {} MeV'.format(getEBeamEnergyL3_MeV(verbose=True))
    print 'E at BC2 is {} MeV'.format(getEBeamEnergyBC2_MeV(verbose=True))
    print 'Q is {} nC'.format(getEBeamCharge_nC(verbose=True))
    print 'I is {} A'.format(getEBeamPkCurrentBC2_A(verbose=True))
    print 'fee is {} mJ'.format(getPulseEnergy_mJ(verbose=True))

